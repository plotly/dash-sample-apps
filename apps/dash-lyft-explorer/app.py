
from dash import Dash, html, dcc, Input, Output, State, callback_context, no_update
import dash_bootstrap_components as dbc
import pandas as pd
import numpy as np
from lyft_dataset_sdk.utils.data_classes import LidarPointCloud

from constants import INITIAL_TOKEN, token_list, lv5, NAME2COLOR
from utils.components import Header, CONTROLS, DECK_CARD
from utils.model import build_figure, build_deck

app = Dash(
    __name__, 
    external_stylesheets=[dbc.themes.CYBORG],
    title="Lyft Interactive Dashboard",
)
server = app.server


app.layout = dbc.Container(
    [
        Header("Dash Lyft Perception", app),
        dbc.Card(dbc.Row([dbc.Col(c) for c in CONTROLS]), body=True),
        dbc.Row(
            [
                dbc.Col(dbc.Card(dcc.Graph(id="graph-camera"), body=True), md=5),
                dbc.Col(DECK_CARD, md=7),
            ],
            className="app-body"
        ),
        dcc.Store(id="sample-token", data=INITIAL_TOKEN),
    ],
    fluid=True,
)


@app.callback(
    Output("progression", "value"),
    Input("prev", "n_clicks"), 
    Input("next", "n_clicks"),
    State("progression", "value"),
)
def update_current_token(btn_prev, btn_next, curr_progress):
    ctx = callback_context
    prop_id = ctx.triggered[0]["prop_id"]

    if "next" in prop_id:
        return min(int(curr_progress) + 1, len(token_list))
    elif "prev" in prop_id:
        return max(0, int(curr_progress) - 1)
    else:
        return no_update


@app.callback(
    Output("graph-camera", "figure"),
    Output("deck-pointcloud", "data"),
    Output("button-group", "children"),
    Output("progression", "type"),
    Input("progression", "value"),
    Input("camera", "value"),
    Input("lidar", "value"),
    Input("overlay", "value"),
    Input("view-mode", "value"),
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

    return fig, r.to_json(), no_update, no_update


if __name__ == "__main__":
    app.run_server(debug=True)
