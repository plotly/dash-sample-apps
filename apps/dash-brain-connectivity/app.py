import dash
import json
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State, ALL
from dash_slicer import VolumeSlicer

import dash_bootstrap_components as dbc
import plotly.express as px
import nibabel as nib
import numpy as np
from nilearn import image
import pandas as pd
import dash_table
from plotly import graph_objects as go

from nibabel.affines import apply_affine
from matplotlib import pyplot as plt

from time import time
from datetime import datetime

t1 = time()
stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, external_stylesheets=stylesheets)
server = app.server

with open("assets/list_of_available_parcellations.json", "r") as f:
    atlases = json.load(f)
atlas_keys = sorted(atlases.keys())
atlas_names = [atlases[key] for key in atlas_keys]

template_img = nib.load("assets/MNI152_T1_1mm.nii.gz")


def get_atlas(atlas_key):
    with open(f"assets/{atlas_key}.json", "r") as f:
        information = json.load(f)
    atlas_table = pd.DataFrame(information["table"])
    atlas_description = information["description"]
    atlas_reference = information["reference"]
    atlas_connmat = np.load(f"assets/{atlas_key}.npy")
    tmp_img = nib.load(f"assets/{atlas_key}.nii.gz")
    atlas_img = image.resample_to_img(tmp_img, template_img, interpolation="nearest")
    parcellation = atlas_img.get_fdata().squeeze().astype(np.uint8)
    return atlas_table, atlas_description, atlas_reference, atlas_connmat, parcellation


def make_heatmap(df, colormap=None):
    axis_ticks = df.columns
    n_elements = len(axis_ticks)
    gap_width = 20 / n_elements
    if n_elements > 100:
        gap_width = 0
    layout = go.Layout(
        yaxis=dict(
            autorange="reversed",
            scaleanchor="x",
            constrain="domain",
            scaleratio=1,
            showgrid=True,
            showticklabels=False,
            tickvals=np.arange(-0.5, df.shape[0], 1),
        ),
        xaxis=dict(
            constrain="domain",
            tickvals=np.arange(-0.5, df.shape[0], 1),
            showgrid=True,
            showticklabels=False,
        ),
        margin=dict(l=0, r=0, b=0, t=0, pad=0),
    )
    heatmap = go.Heatmap(
        z=df,
        zmin=0,
        zmax=1,
        x=axis_ticks,
        y=axis_ticks,
        xgap=gap_width,
        ygap=gap_width,
        hoverongaps=False,
        autocolorscale=True if colormap is None else False,
        colorscale=colormap,
    )
    fig = go.Figure(heatmap, layout)
    fig.update_traces(showscale=True)

    return fig


# Load the initial data
(table, description, reference, conn_mat, parcellation) = get_atlas(atlas_keys[0])
region_labels = table["Region Name"]
vol = template_img.get_fdata().squeeze()

parcellation_colormap = px.colors.qualitative.Alphabet
connectivity_colormap = [
    [
        val,
        px.colors.label_rgb(
            plt.cm.plasma(val, alpha=1 if not val == 0 else 0, bytes=True)
        ),
    ]
    for val in np.linspace(0, 1, 100)
]
sagittal_slicer = VolumeSlicer(app, vol, axis=0, scene_id="brain", color="white")
coronal_slicer = VolumeSlicer(app, vol, axis=1, scene_id="brain", color="white")
axial_slicer = VolumeSlicer(app, vol, axis=2, scene_id="brain", color="white")

# Update all slicers in order
slicer_list = []
for slicer_idx, slicer in enumerate([sagittal_slicer, coronal_slicer, axial_slicer]):
    slicer.graph.figure.update_layout(dragmode=False, plot_bgcolor="rgb(0, 0, 0)")
    slicer.overlay_data.data = slicer.create_overlay_data(
        parcellation, parcellation_colormap
    )
    aux_slider = dcc.Slider(id=f"aux-slider-{slicer_idx}", max=slicer.nslices)
    slicer_list.append(
        html.Div(
            [
                slicer.graph,
                html.Div([slicer.slider, aux_slider], style={"display": "none"}),
                *slicer.stores,
            ]
        )
    )

setpos_store = dcc.Store(
    id={"context": "app", "scene": sagittal_slicer.scene_id, "name": "setpos"}
)

t2 = time()
print(t2 - t1)

# Define Navbar
# button_gh = "Hello"
button_gh = dbc.Button(
    "View Code on github",
    outline=True,
    color="primary",
    href="https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-image-annotation",
    id="gh-link",
    style={"text-transform": "none"},
)
navbar = dbc.Navbar(
    dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        html.A(
                            html.Img(
                                src=app.get_asset_url("dash-logo-new.png"),
                                height="30px",
                            ),
                            href="https://plot.ly",
                        )
                    ),
                    dbc.Col(dbc.NavbarBrand("Brain Functional Connecitivity App")),
                ],
                align="center",
            ),
            dbc.Row(
                dbc.Col(
                    [
                        dbc.NavbarToggler(id="navbar-toggler"),
                        dbc.Collapse(
                            dbc.Nav(
                                [dbc.NavItem(button_gh)],  # dbc.NavItem(button_howto),
                                className="ml-auto",
                                navbar=True,
                            ),
                            id="navbar-collapse",
                            navbar=True,
                        ),
                        # modal_overlay,
                    ]
                ),
                align="center",
            ),
        ],
        fluid=True,
    ),
    color="dark",
    dark=True,
    className="mb-5",
)

# Layout Slicer Positions
saggital_slicer_card = dbc.Card(
    [dbc.CardHeader("Saggital View"), dbc.CardBody(dbc.CardBody(slicer_list[0])),],
    color="dark",
    inverse=True,
)
coronal_slicer_card = dbc.Card(
    [dbc.CardHeader("Coronal View"), dbc.CardBody(slicer_list[1]),],
    color="dark",
    inverse=True,
)

axial_slicer_card = dbc.Card(
    [dbc.CardHeader("Axial View"), dbc.CardBody(dbc.CardBody(slicer_list[2])),],
    color="dark",
    inverse=True,
)

nav_table = html.Div(
    [
        dbc.Label("Brain viewer coordinates", html_for="nav-table-title"),
        dbc.Row(
            [dbc.Col("", width=1), dbc.Col("Voxel space"), dbc.Col("MNI space")],
            id="nav-table-title",
        ),
        dbc.Row(
            [
                dbc.Col("X", style={"text-align": "right"}, width=1),
                dbc.Col(
                    dcc.Input(
                        id="x-vox",
                        type="number",
                        placeholder="X value",
                        style={"width": "100%"},
                    )
                ),
                dbc.Col(
                    dcc.Input(
                        id="x-world",
                        type="number",
                        disabled=True,
                        placeholder="X world",
                        style={"width": "100%"},
                    )
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col("Y", style={"text-align": "right"}, width=1),
                dbc.Col(
                    dcc.Input(
                        id="y-vox",
                        type="number",
                        placeholder="Y value",
                        style={"width": "100%"},
                    )
                ),
                dbc.Col(
                    dcc.Input(
                        id="y-world",
                        type="number",
                        disabled=True,
                        placeholder="Y world",
                        style={"width": "100%"},
                    )
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col("Z", style={"text-align": "right"}, width=1),
                dbc.Col(
                    dcc.Input(
                        id="z-vox",
                        type="number",
                        placeholder="Z value",
                        style={"width": "100%"},
                    )
                ),
                dbc.Col(
                    dcc.Input(
                        id="z-world",
                        type="number",
                        disabled=True,
                        placeholder="Z world",
                        style={"width": "100%"},
                    )
                ),
            ]
        ),
    ],
    style={"width": "70%"},
)

connectivity_info = html.Div(
    [
        dbc.Row(
            [
                dbc.Col("The current region: ", width=3),
                dbc.Col(html.Pre(id="position-region-info"), width="auto",),
            ],
        ),
        dbc.Row(
            [
                dbc.Col("is connected to: ", width=3),
                dbc.Col(html.Pre(id="position-conn-info")),
            ],
        ),
        dbc.Row(
            [
                dbc.Col("with a strength of: ", width=3),
                dbc.Col(html.Pre(id="connection-conn-info"), width="auto",),
            ],
        ),
    ],
    id="connectivity-info",
)

atlas_drop = dcc.Dropdown(
    id="atlas-drop-menu",
    options=[
        {"label": atlases[col_name], "value": col_name} for col_name in atlas_keys
    ],
    value=atlas_keys[0],
)


viewer_card = dbc.Card(
    [
        dbc.CardBody(
            [
                dbc.Row(
                    dbc.Col(
                        dbc.CardDeck(
                            [
                                saggital_slicer_card,
                                coronal_slicer_card,
                                axial_slicer_card,
                            ]
                        )
                    ),
                    style={"margin-top": "1em"},
                ),
                dbc.Row(
                    dbc.Col(html.H3("Brain viewer context information")),
                    style={"margin-top": "1em"},
                ),
                dbc.Row(dbc.Col(html.Hr())),
                dbc.Row(
                    [
                        dbc.Col(nav_table, width=4),
                        dbc.Col(
                            [
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            [
                                                dbc.Label("Viewer mode: "),
                                                dbc.Badge(id="view-mode-alert"),
                                            ],
                                            width="auto",
                                        ),
                                    ],
                                ),
                                connectivity_info,
                            ],
                            width=8,
                        ),
                    ]
                ),
            ]
        ),
        dbc.CardFooter(
            dbc.Row(
                [
                    dbc.Col("Please select a brain parcellation: "),
                    dbc.Col(atlas_drop, width=8),
                ]
            )
        ),
    ],
    color="light",
    inverse=False,
)

table_card = dbc.Card(
    [
        dbc.CardHeader("Data Table"),
        dbc.CardBody(
            id="table_container",
            children=[
                html.Div(
                    dash_table.DataTable(
                        id="region_table",
                        columns=[{"name": col, "id": col} for col in table.columns],
                        data=table.to_dict("records"),
                        fixed_rows={"headers": True},
                        style_table={
                            "overflowX": "auto",
                            "maxHeight": "40vh",
                            "minWidth": "100%",
                        },
                        style_cell={
                            "textOverflow": "ellipsis",
                            "maxWidth": "300px",
                            "minWidth": "80px",
                            "overflow": "hidden",
                        },
                        style_cell_conditional=[
                            {"if": {"column_id": "Region Name"}, "width": "300px"},
                        ],
                        tooltip_data=[
                            {
                                column: {"value": str(value), "type": "markdown"}
                                for column, value in row.items()
                                if column == "Region Name"
                            }
                            for row in table.to_dict("records")
                        ],
                        tooltip_duration=None,
                    )
                )
            ],
        ),
        dbc.CardFooter(),
    ],
    color="light",
    inverse=False,
)

connectivity_card = dbc.Card(
    [
        dbc.CardHeader("Functional Connectivity"),
        dbc.CardBody([dcc.Graph(id="conn_mat"),]),
        dbc.CardFooter(
            [
                html.Div(
                    "Select the connectivity threshold and thresholding method below:"
                ),
                html.Div(dcc.Slider(id="val_slider")),
                dbc.RadioItems(
                    options=[
                        {"label": "absolute threshold", "value": "absolute"},
                        {"label": "percentage threshold", "value": "percentage"},
                    ],
                    value="percentage",
                    id="threshold-mode",
                ),
            ]
        ),
        dcc.Store(id="store_region_names"),
        dcc.Store(id="store_conn_mat"),
        dcc.Store(id="store_parcellation", data=parcellation),
        dcc.Store(id="store_conn_thr_mat"),
        dcc.Store(id="store_slicer_region"),
        dcc.Store(id="store_graph_position"),
        dcc.Store(id="store_graph_region"),
        dcc.Store(id="store_highlight_region"),
        setpos_store,
    ]
)

# Define the layout
app.layout = html.Div(
    [
        navbar,
        dbc.Container(
            [
                dbc.Row(dbc.Col(viewer_card)),
                dbc.Row(
                    dbc.Col(dbc.CardDeck([table_card, connectivity_card])),
                    style={"margin-top": "1em", "max-height": "min-content"},
                ),
            ],
            fluid=True,
        ),
    ]
)


# TODO move this to the helper functions?
def map_color(vector):
    # We will ignore all nonzero or nan values here because dash slicer adds a zero in front
    values = np.unique(vector[np.logical_not(np.isnan(vector))])
    mask = np.zeros(vector.shape, dtype=np.uint8)
    colormap = []
    for vid, value in enumerate(values):
        mask[vector == value] = vid + 1
        if value == 0:
            colormap.append(plt.cm.plasma(value, alpha=0, bytes=True))
        else:
            colormap.append(plt.cm.plasma(value, alpha=1, bytes=True))

    return mask, colormap


@app.callback(
    [
        Output("store_region_names", "data"),
        Output("store_conn_mat", "data"),
        Output("region_table", "data"),
        Output("region_table", "tooltip_data"),
        Output("store_parcellation", "data"),
    ],
    Input("atlas-drop-menu", "value"),
)
def update_atlas(atlas_selection):
    now = datetime.now()
    print("cb 0", now.minute, now.second)
    (table, description, reference, conn_mat, parcellation) = get_atlas(atlas_selection)
    region_labels = table["Region Name"]
    region_tooltips = [
        {
            column: {"value": str(value), "type": "markdown"}
            for column, value in row.items()
            if column == "Region Name"
        }
        for row in table.to_dict("records")
    ]
    return (
        region_labels,
        conn_mat,
        table.to_dict("records"),
        region_tooltips,
        parcellation,
    )


@app.callback(
    [Output("view-mode-alert", "children"), Output("view-mode-alert", "color")],
    Input("store_graph_region", "data"),
)
def announce_viewmode(graph_data):
    now = datetime.now()
    print("cb 1", now.minute, now.second)
    view_mode = "Brain atlas view mode"
    color = "primary"
    if graph_data is not None:
        view_mode = "Functional Connectivity view mode"
        color = "success"
    return view_mode, color


@app.callback(
    [
        Output("position-region-info", "children"),
        Output("position-conn-info", "children"),
        Output("connection-conn-info", "children"),
    ],
    [
        Input("store_graph_region", "data"),
        Input("store_slicer_region", "data"),
        Input("store_region_names", "data"),
    ],
    [State("store_conn_mat", "data")],
)
def position_info(seed_region, brain_region, region_names, raw_conn_mat):
    now = datetime.now()
    print("cb 2", now.minute, now.second)
    conn_val = "Select a seed region to display its connectivity here"
    conn_mat = np.array(raw_conn_mat)
    if seed_region is not None:
        conn_val = "Point the cursor at a brain region to display its connectivity with the seed region."
        if brain_region is not None:
            seed_id = region_names.index(seed_region)
            brain_id = region_names.index(brain_region)
            conn_val = f"r = {conn_mat[seed_id, brain_id]:.2f}"
    if seed_region is None:
        seed_region = "Select a seed region in the table or heatmap below."
    if brain_region is None:
        brain_region = (
            "Not inside a region. Point the cursor at a colored brain region above."
        )

    return brain_region, seed_region, conn_val


# Update the slicer app overlay so it shows the connectivity information
@app.callback(
    [
        Output(axial_slicer.overlay_data.id, "data"),
        Output(sagittal_slicer.overlay_data.id, "data"),
        Output(coronal_slicer.overlay_data.id, "data"),
    ],
    Input("store_graph_region", "data"),
    Input("store_conn_thr_mat", "data"),
    Input("store_parcellation", "data"),
    [State("region_table", "data"), State("threshold-mode", "value"),],
)
def update_brain_overlay(
    graph_region, conn_mat_data, parcellation_data, table_data, threshold,
):
    now = datetime.now()
    print("cb 3", now.minute, now.second)
    conn_mat = np.array(conn_mat_data)
    parcellation = np.array(parcellation_data).astype(np.uint8)
    overlay = None
    cm = None
    if graph_region is None:
        overlay = parcellation
        cm = parcellation_colormap
    else:
        table = pd.DataFrame(table_data)
        region_id = np.where(table["Region Name"] == graph_region)[0][0]
        # Prepend a zero for the background
        conn_vec = np.insert(
            np.array(conn_mat[region_id, :], dtype=np.float), 0, np.nan
        )
        # Masking out all negative values (there shouldn't be many).
        conn_vec[conn_vec < 0] = np.nan
        mask, cm = map_color(conn_vec)
        overlay = mask[parcellation]
    axial_overlay = axial_slicer.create_overlay_data(overlay, cm)
    sagittal_overlay = sagittal_slicer.create_overlay_data(overlay, cm)
    coronal_overlay = coronal_slicer.create_overlay_data(overlay, cm)
    return axial_overlay, sagittal_overlay, coronal_overlay


# Update the nav table values to show the current slicer positions
# Note that we are listening to the "drag_value" attributes here to avoid circularity of the callbacks
@app.callback(
    [Output("x-vox", "value"), Output("y-vox", "value"), Output("z-vox", "value")],
    [
        Input(sagittal_slicer.slider.id, "drag_value"),
        Input(coronal_slicer.slider.id, "drag_value"),
        Input(axial_slicer.slider.id, "drag_value"),
        Input("store_graph_position", "data"),
    ],
)
def update_navtable(x_pos, y_pos, z_pos, requested_position):
    now = datetime.now()
    print("cb 4", now.minute, now.second)
    ctx = dash.callback_context
    if ctx.triggered[0]["prop_id"] == "store_graph_position.data":
        x_pos, y_pos, z_pos = np.floor(
            apply_affine(np.linalg.inv(template_img.affine), requested_position)
        ).astype(int)
    return x_pos, y_pos, z_pos


# Listen for changes to the nav table values and update the slicer position accordingly
@app.callback(
    [
        Output(setpos_store.id, "data"),
        Output("x-world", "value"),
        Output("y-world", "value"),
        Output("z-world", "value"),
        Output("store_slicer_region", "data"),
    ],
    [
        Input("x-vox", "value"),
        Input("y-vox", "value"),
        Input("z-vox", "value"),
        Input("store_parcellation", "data"),
    ],
    State("region_table", "data"),
)
def write_table_values_to_slicer(x_vox, y_vox, z_vox, parcellation_data, table_data):
    # Find the overlay value at the current slicer position
    now = datetime.now()
    print("cb 5", now.minute, now.second)
    parcellation = np.array(parcellation_data)
    region_number = parcellation[x_vox, y_vox, z_vox]
    # Compute the world positions
    x_world, y_world, z_world = apply_affine(template_img.affine, (x_vox, y_vox, z_vox))
    if not region_number == 0:
        table = pd.DataFrame(table_data)
        region_name = table.query("`Region ID` == @region_number")[
            "Region Name"
        ].values[0]
    else:
        region_name = None
    return (
        (z_vox, y_vox, x_vox),
        x_world,
        y_world,
        z_world,
        region_name,
    )


# React to clicks on graphs and update the corresponding region
@app.callback(
    [Output("store_graph_region", "data"), Output("store_graph_position", "data")],
    [Input("conn_mat", "clickData"), Input("region_table", "active_cell")],
    [State("region_table", "data"), State("store_graph_region", "data")],
)
def listen_to_graph_clicks(
    heatmap_clickdata, active_table_cell, table_data, previous_region
):
    now = datetime.now()
    print("cb 6", now.minute, now.second)
    ctx = dash.callback_context
    table = pd.DataFrame(table_data)
    region_name = None
    if ctx.triggered[0]["prop_id"] == "conn_mat.clickData":
        region_name = heatmap_clickdata["points"][0]["y"]
    elif ctx.triggered[0]["prop_id"] == "region_table.active_cell":
        region_name = table.iloc[active_table_cell["row"]]["Region Name"]
    if region_name == previous_region:
        return None, dash.no_update
    else:
        # Find the corresponding position
        region_position = tuple(
            table.query("`Region Name` == @region_name")[["X", "Y", "Z"]].values[0]
        )
        return region_name, region_position


# Set the threshold slider according to the threshold mode
@app.callback(
    [
        Output("val_slider", "min"),
        Output("val_slider", "max"),
        Output("val_slider", "step"),
        Output("val_slider", "marks"),
        Output("val_slider", "value"),
    ],
    Input("threshold-mode", "value"),
)
def select_thresholding(val):
    now = datetime.now()
    print("cb 7", now.minute, now.second)
    min_val = 0
    max_val = 1
    step = 0.1
    ticks = np.arange(min_val, max_val + step, step)

    if val == "absolute":
        # TODO: check if we still need this complex rounding for float precision
        marks = {
            key: val
            for key, val in [
                (int(i), f"{int(i)}") if i.is_integer() else (i, f"{i:.1f}")
                for i in map(lambda x: np.round(x, 1), ticks)
            ]
        }

        value = 0.1
        return min_val, max_val, step, marks, value
    elif val == "percentage":
        marks = {
            key: val
            for key, val in [
                (int(i), f"{int(i * 100)}%")
                if i.is_integer()
                else (i, f"{int(i * 100)}%")
                for i in ticks
            ]
        }
        value = 0.5
        return min_val, max_val, step, marks, value
    else:
        raise Exception(f"received unknown threshold mode {val}")


@app.callback(
    Output("store_highlight_region", "data"),
    [Input("store_graph_region", "data"), Input("store_slicer_region", "data")],
)
def update_highlight_region(graph_region, slicer_region):
    now = datetime.now()
    print("cb 8", now.minute, now.second)
    ctx = dash.callback_context
    if ctx.triggered[0]["prop_id"] == "store_graph_region.data":
        return graph_region
    elif ctx.triggered[0]["prop_id"] == "store_slicer_region.data":
        return slicer_region
    else:
        raise Exception("This is not supposed to happen")


@app.callback(
    Output("store_conn_thr_mat", "data"),
    [Input("val_slider", "value"), Input("store_conn_mat", "data")],
    [State("threshold-mode", "value")],
)
def update_connmat(thr_value, raw_conn_mat, thr_mode):
    now = datetime.now()
    print("cb 9", now.minute, now.second)
    conn_mat = np.array(raw_conn_mat)
    if thr_mode == "percentage":
        conn_vec = np.abs(conn_mat[np.tril_indices(conn_mat.shape[0], -1)])
        thr_value = np.percentile(conn_vec, thr_value * 100)

    conn_mat[np.abs(conn_mat) < thr_value] = np.nan
    return conn_mat


@app.callback(
    Output("conn_mat", "figure"),
    [
        Input("store_graph_region", "data"),
        Input("store_slicer_region", "data"),
        Input("store_conn_thr_mat", "data"),
    ],
    State("store_region_names", "data"),
)
def connectivity_heatmap(graph_region, slicer_region, conn_mat, region_names):
    now = datetime.now()
    print("cb 10", now.minute, now.second)
    conn_mat = np.array(conn_mat, dtype=np.float)
    conn_mat[conn_mat < 0] = 0
    conn_mat_df = pd.DataFrame(conn_mat, columns=region_names, index=region_names)
    # Redraw the figure
    fig = make_heatmap(conn_mat_df, colormap=connectivity_colormap)

    if graph_region is not None:
        # We have a selected region and are in connectivity drawing mode
        fig = add_region_shape(fig, graph_region, line=dict(color="#c3e6cb", width=4))

    if slicer_region is not None:
        fig = add_region_shape(
            fig,
            slicer_region,
            line=dict(color="#b8daff", width=4),
            orientation="v" if graph_region is not None else "h",
        )
    return fig


def add_region_shape(fig, region_id, line=None, orientation="h"):
    # We just use the position along the y axis here but because the matrix is symmetrical, there is no difference
    draw_pos = np.where(np.array(fig.data[0]["y"]) == region_id)[0][0]
    fig.add_shape(
        type="rect",
        xref="x",
        yref="y",
        x0=-0.5 if orientation == "h" else draw_pos - 0.5,
        y0=draw_pos - 0.5 if orientation == "h" else -0.5,
        x1=len(fig.data[0]["x"]) - 0.5 if orientation == "h" else draw_pos + 0.5,
        y1=draw_pos + 0.5 if orientation == "h" else len(fig.data[0]["y"]) - 0.5,
        line=line if line is not None else dict(color="LightSeaGreen", width=3),
    )
    return fig


# Highlight the table region
@app.callback(
    Output("region_table", "style_data_conditional"),
    [Input("store_graph_region", "data"), Input("store_slicer_region", "data")],
    prevent_initial_call=True,
)
def higlight_row(graph_region, slicer_region):
    """
    When a region is selected, highlight the corresponding table row.
    """
    now = datetime.now()
    print("cb 11", now.minute, now.second)
    style = [
        {
            "if": {"filter_query": "{Region Name} = '%s'" % graph_region},
            "backgroundColor": "#d4edda",
        },
        {
            "if": {"filter_query": "{Region Name} = '%s'" % slicer_region},
            "backgroundColor": "#cce5ff",
            "color": "#004085",
            "fontWeight": "bold",
        },
    ]
    return style


if __name__ == "__main__":
    app.run_server(debug=True, dev_tools_props_check=False)
