import os
import time

import colorlover as cl
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_deck
from dash.dependencies import Input, Output, State
from lyft_dataset_sdk.lyftdataset import LyftDataset, LyftDatasetExplorer
from lyft_dataset_sdk.utils.data_classes import Box, LidarPointCloud, RadarPointCloud
import numpy as np
import pandas as pd
from PIL import Image
import plotly.graph_objects as go
import plotly.express as px
import pydeck as pdk


def Header(name, app):
    title = html.H2(name, style={"margin-top": 5})
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 60}
    )
    link = html.A(logo, href="https://plotly.com/dash/")

    return dbc.Row([dbc.Col(title, md=8), dbc.Col(link, md=4)])


def unsnake(st):
    """BECAUSE_WE_DONT_READ_LIKE_THAT"""
    return st.replace("_", " ").title()


def build_deck(mode, pc_df, polygon_data):
    if mode == "first_person":
        view = pdk.View(type="FirstPersonView", controller=True)
        view_state = pdk.ViewState(latitude=0, longitude=0, bearing=-90, pitch=15)
        point_size = 10
    elif mode == "orbit":
        view = pdk.View(type="OrbitView", controller=True)
        view_state = pdk.ViewState(
            target=[0, 0, 1e-5],
            controller=True,
            zoom=23,
            rotation_orbit=-90,
            rotation_x=15,
        )
        point_size = 3

    else:
        view_state = pdk.ViewState(
            latitude=0,
            longitude=0,
            bearing=45,
            pitch=50,
            zoom=20,
            max_zoom=30,
            position=[0, 0, 1e-5],
        )
        view = pdk.View(type="MapView", controller=True)
        point_size = 1

    pc_layer = pdk.Layer(
        "PointCloudLayer",
        data=pc_df,
        get_position=["x", "y", "z"],
        get_color=[255, 255, 255],
        auto_highlight=True,
        pickable=False,
        point_size=point_size,
        coordinate_system=2,
        coordinate_origin=[0, 0],
    )

    box_layer = pdk.Layer(
        "PolygonLayer",
        data=polygon_data,
        stroked=True,
        pickable=True,
        filled=True,
        extruded=True,
        opacity=0.2,
        wireframe=True,
        line_width_min_pixels=1,
        get_polygon="polygon",
        get_fill_color="color",
        get_line_color=[255, 255, 255],
        get_line_width=0,
        coordinate_system=2,
        get_elevation="elevation",
    )

    tooltip = {"html": "<b>Label:</b> {name}"}

    r = pdk.Deck(
        [pc_layer, box_layer],
        initial_view_state=view_state,
        views=[view],
        tooltip=tooltip,
        map_provider=None,
    )

    return r


def compute_pointcloud_for_image(
    lv5,
    sample_token: str,
    dot_size: int = 2,
    pointsensor_channel: str = "LIDAR_TOP",
    camera_channel: str = "CAM_FRONT",
    out_path: str = None,
):
    """Scatter-plots a point-cloud on top of image.
    Args:
        sample_token: Sample token.
        dot_size: Scatter plot dot size.
        pointsensor_channel: RADAR or LIDAR channel name, e.g. 'LIDAR_TOP'.
        camera_channel: Camera channel name, e.g. 'CAM_FRONT'.
        out_path: Optional path to save the rendered figure to disk.
    Returns:
        tuple containing the points, array of colors and a pillow image
    """
    sample_record = lv5.get("sample", sample_token)

    # Here we just grab the front camera and the point sensor.
    pointsensor_token = sample_record["data"][pointsensor_channel]
    camera_token = sample_record["data"][camera_channel]

    points, coloring, im = lv5.explorer.map_pointcloud_to_image(
        pointsensor_token, camera_token
    )

    return points, coloring, im


def render_box_in_image(lv5, im, sample: str, camera_channel: str):
    camera_token = sample["data"][camera_channel]
    data_path, boxes, camera_intrinsic = lv5.get_sample_data(
        camera_token, flat_vehicle_coordinates=False
    )

    arr = np.array(im)

    for box in boxes:
        c = NAME2COLOR[box.name]
        box.render_cv2(arr, normalize=True, view=camera_intrinsic, colors=(c, c, c))

    new = Image.fromarray(arr)
    return new


def get_token_list(scene):
    token_list = [scene["first_sample_token"]]
    sample = lv5.get("sample", token_list[0])

    while sample["next"] != "":
        token_list.append(sample["next"])
        sample = lv5.get("sample", sample["next"])

    return token_list


def build_figure(lv5, sample, lidar, camera, overlay):
    points, coloring, im = compute_pointcloud_for_image(
        lv5, sample["token"], pointsensor_channel=lidar, camera_channel=camera
    )

    if "boxes" in overlay:
        im = render_box_in_image(lv5, im, sample, camera_channel=camera)

    fig = px.imshow(im, binary_format="jpeg", binary_compression_level=2)

    if "pointcloud" in overlay:
        fig.add_trace(
            go.Scattergl(
                x=points[0,],
                y=points[1,],
                mode="markers",
                opacity=0.4,
                marker_color=coloring,
                marker_size=3,
            )
        )

    fig.update_layout(
        margin=dict(l=10, r=10, t=0, b=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        hovermode=False,
    )
    fig.update_xaxes(showticklabels=False, showgrid=False, range=(0, im.size[0]))
    fig.update_yaxes(showticklabels=False, showgrid=False, range=(im.size[1], 0))

    return fig


# Variables
CAMERAS = [
    "CAM_FRONT",
    "CAM_BACK",
    "CAM_FRONT_ZOOMED",
    "CAM_FRONT_LEFT",
    "CAM_FRONT_RIGHT",
    "CAM_BACK_RIGHT",
    "CAM_BACK_LEFT",
]
LIDARS = ["LIDAR_TOP", "LIDAR_FRONT_RIGHT", "LIDAR_FRONT_LEFT"]

NAME2COLOR = dict(
    zip(
        ["bus", "car", "other_vehicle", "pedestrian", "truck"],
        cl.to_numeric(cl.scales["5"]["div"]["Spectral"]),
    )
)

# Create Lyft object
lv5 = LyftDataset(data_path="./data", json_path="./data/train_data", verbose=True)
# Load a single scene
scene = lv5.scene[0]
token_list = get_token_list(scene)
INITIAL_TOKEN = scene["first_sample_token"]


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.CYBORG])
server = app.server

controls = [
    dbc.FormGroup(
        [
            dbc.Label("Camera Position"),
            dbc.Select(
                id="camera",
                options=[
                    {"label": unsnake(s.replace("CAM_", "")), "value": s}
                    for s in CAMERAS
                ],
                value=CAMERAS[0],
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Image Overlay"),
            dbc.Checklist(
                id="overlay",
                value=[],
                options=[
                    {"label": x.title(), "value": x} for x in ["pointcloud", "boxes"]
                ],
                inline=True,
                switch=True,
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Frame"),
            html.Br(),
            dbc.Spinner(
                dbc.ButtonGroup(
                    [
                        dbc.Button(
                            "Prev", id="prev", n_clicks=0, color="primary", outline=True
                        ),
                        dbc.Button("Next", id="next", n_clicks=0, color="primary"),
                    ],
                    id="button-group",
                    style={"width": "100%"},
                ),
                spinner_style={"margin-top": 0, "margin-bottom": 0},
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Progression"),
            dbc.Spinner(
                dbc.Input(
                    id="progression", type="range", min=0, max=len(token_list), value=0
                ),
                spinner_style={"margin-top": 0, "margin-bottom": 0},
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Lidar Position"),
            dbc.Select(
                id="lidar",
                value=LIDARS[0],
                options=[
                    {"label": unsnake(s.replace("LIDAR_", "")), "value": s}
                    for s in LIDARS
                ],
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Lidar View Mode"),
            dbc.Select(
                id="view-mode",
                value="map",
                options=[
                    {"label": unsnake(x), "value": x}
                    for x in ["first_person", "orbit", "map"]
                ],
            ),
        ]
    ),
]

deck_card = dbc.Card(
    dash_deck.DeckGL(id="deck-pointcloud", tooltip={"html": "<b>Label:</b> {name}"}),
    body=True,
    style={"height": "calc(95vh - 215px)"},
)

app.layout = dbc.Container(
    [
        Header("Dash Lyft Perception", app),
        html.Br(),
        dbc.Card(dbc.Row([dbc.Col(c) for c in controls], form=True), body=True),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(dbc.Card(dcc.Graph(id="graph-camera"), body=True), md=5),
                dbc.Col(deck_card, md=7),
            ]
        ),
        dcc.Store(id="sample-token", data=INITIAL_TOKEN),
    ],
    fluid=True,
)


@app.callback(
    Output("progression", "value"),
    [Input("prev", "n_clicks"), Input("next", "n_clicks")],
    [State("progression", "value")],
)
def update_current_token(btn_prev, btn_next, curr_progress):
    ctx = dash.callback_context
    prop_id = ctx.triggered[0]["prop_id"]

    if "next" in prop_id:
        return min(int(curr_progress) + 1, len(token_list))
    elif "prev" in prop_id:
        return max(0, int(curr_progress) - 1)
    else:
        return dash.no_update


@app.callback(
    [
        Output("graph-camera", "figure"),
        Output("deck-pointcloud", "data"),
        Output("button-group", "children"),
        Output("progression", "type"),
    ],
    [
        Input("progression", "value"),
        Input("camera", "value"),
        Input("lidar", "value"),
        Input("overlay", "value"),
        Input("view-mode", "value"),
    ],
)
def update_graphs(progression, camera, lidar, overlay, view_mode):
    token = token_list[int(progression)]
    sample = lv5.get("sample", token)
    pointsensor_token = sample["data"][lidar]
    pointsensor = lv5.get("sample_data", pointsensor_token)
    pc = LidarPointCloud.from_file(lv5.data_path / pointsensor["filename"])
    pc_df = pd.DataFrame(pc.points.T, columns=["x", "y", "z", "intensity"])

    if lidar in ["LIDAR_FRONT_LEFT", "LIDAR_FRONT_RIGHT"]:
        pc_df.z = -pc_df.z + 1

    _, boxes, camera_intrinsic = lv5.get_sample_data(
        pointsensor["token"], flat_vehicle_coordinates=False
    )

    polygon_data = [
        {
            "name": box.name,
            "polygon": box.bottom_corners().T.tolist(),
            "width": box.wlh[0],
            "length": box.wlh[1],
            "elevation": box.wlh[2],
            "color": NAME2COLOR[box.name],
            "token": box.token,
            "distance": np.sqrt(np.square(boxes[0].center).sum()),
        }
        for box in boxes
    ]

    # Build figure and pydeck object
    fig = build_figure(lv5, sample, lidar, camera, overlay)
    r = build_deck(view_mode, pc_df, polygon_data)

    return fig, r.to_json(), dash.no_update, dash.no_update


if __name__ == "__main__":
    app.run_server(debug=True)
