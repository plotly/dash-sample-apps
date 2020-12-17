import dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State, ALL
from dash_slicer import VolumeSlicer

import dash_bootstrap_components as dbc
import plotly.express as px
import nibabel as nib
import numpy as np
from nilearn import datasets, image
import pandas as pd
import dash_table
import dash_bio as dashbio
#from plotly import graph_objects as go

from nibabel.affines import apply_affine
import pathlib

import sys

sys.path.append("./")
import helper_functions as hf

stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, external_stylesheets=stylesheets)
server = app.server

APP_PATH = pathlib.Path(__file__).parent.resolve()
nilearn_datapath = str(APP_PATH / "nilearn_datasets")
conn_mat = np.load("./assets/yeo_17_connectivity_matrix.npy")
table = pd.read_csv("./assets/yeo_17_region_table.csv")
region_labels = table["Region Name"]

img = nib.load("./assets/MNI152_T1_1mm_brain.nii.gz")
yeo = datasets.fetch_atlas_yeo_2011(data_dir=nilearn_datapath)
parcellation_img = nib.load(yeo["thick_17"])
parc_img_res = image.resample_to_img(parcellation_img, img, interpolation="nearest")
parcellation = parc_img_res.get_fdata().squeeze().astype(np.uint8)
data = img.get_fdata().squeeze()
vol = data
origin = img.affine[:3, 3]
origin = [-90, -126, -72]  # first position has inverse sign
origin = [0, 0, 0]
spacing = np.diag(img.affine[:3, :3])
spacing = [1, 1, 1]

colormap = px.colors.qualitative.Alphabet
saggital_slicer = VolumeSlicer(
    app, vol, spacing=spacing, origin=origin, axis=0, scene_id="brain"
)
coronal_slicer = VolumeSlicer(
    app, vol, spacing=spacing, origin=origin, axis=1, scene_id="brain"
)
axial_slicer = VolumeSlicer(
    app, vol, spacing=spacing, origin=origin, axis=2, scene_id="brain"
)

# Update all slicers in order
# TODO figure out if we need to explicitly define these in order for the Cards to be equal height
slicer_list = []
for slicer_idx, slicer in enumerate([saggital_slicer, coronal_slicer, axial_slicer]):
    slicer.graph.figure.update_layout(dragmode=False, plot_bgcolor="rgb(0, 0, 0)")
    slicer.overlay_data.data = slicer.create_overlay_data(parcellation, colormap)
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
    id={"context": "app", "scene": saggital_slicer.scene_id, "name": "setpos"}
)

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
        dbc.Row([dbc.Col(""), dbc.Col("Voxel Position")]),
        dbc.Row(
            [
                dbc.Col("X"),
                dbc.Col(dcc.Input(id="x-vox", type="number", placeholder="X value", style={"width": "100%"})),
                dbc.Col(dcc.Input(id="x-world", type="number", placeholder="X world", style={"width": "100%"})),
            ]
        ),
        dbc.Row(
            [
                dbc.Col("Y"),
                dbc.Col(dcc.Input(id="y-vox", type="number", placeholder="Y value", style={"width": "100%"})),
                dbc.Col(dcc.Input(id="y-world", type="number", placeholder="Y world", style={"width": "100%"})),
            ]
        ),
        dbc.Row(
            [
                dbc.Col("Z"),
                dbc.Col(dcc.Input(id="z-vox", type="number", placeholder="Z value", style={"width": "100%"})),
                dbc.Col(dcc.Input(id="z-world", type="number", placeholder="Z world", style={"width": "100%"})),
            ]
        ),
    ],
    style={"width": "40%"},
)


controls_card = dbc.Card(
    [
        dbc.CardHeader("Controls"),
        dbc.CardBody(
            dbc.Row(
                [
                    dbc.Col(nav_table),
                    dbc.Col(
                        [
                            html.Pre(id="write_selected_position"),
                            html.Pre(id="write_slicer_region"),
                            html.Pre(id="write_graph_region"),
                            html.Pre(id="write_highlight_region"),
                        ]
                    ),
                ]
            )
        ),
        dbc.CardFooter(
            [
                html.Div(dcc.Slider(id="val_slider"), style={"width": "50%"}),
                html.Div(id="slider-val", children=""),
                dcc.RadioItems(
                    options=[
                        {"label": "Absolute threshold", "value": "absolute"},
                        {"label": "Percentage threshold", "value": "percentage"},
                    ],
                    value="percentage",
                    id="threshold-mode",
                ),
                html.Div(id="dumpster"),
            ]
        ),
    ],
    color="light",
    inverse=False,
    style={"width": "80%"},
)

table_card = dbc.Card(
    [
        dbc.CardHeader("Data Table"),
        dbc.CardBody(
            id="table_container",
            children=dash_table.DataTable(
                id="region_table",
                columns=[{"name": col, "id": col} for col in table.columns],
                data=table.to_dict("records"),
                style_table={"overflowY": "scroll"},
            ),
        ),
    ],
    color="light",
    inverse=False,
)

heatmap_tab = dbc.Card(dbc.CardBody([dcc.Graph(id="conn_mat")]))

circos_tab = dbc.Card(
    [dbc.CardBody(id="circos_graph"), dbc.CardFooter(html.Div(id="circos_footer"))]
)

connectivity_tabs = dbc.Tabs(
    [
        dbc.Tab(heatmap_tab, label="Heatmap Connectivity"),
        dbc.Tab(circos_tab, label="Circos Connectivity"),
    ]
)
connectivity_card = dbc.Card(
    [
        dbc.CardHeader("Functional Connectivity"),
        dbc.CardBody(connectivity_tabs),
        dcc.Store(id="store_region_names", data=region_labels),
        dcc.Store(id="store_conn_mat", data=conn_mat),
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
                dbc.Row(
                    [
                        dbc.Col(saggital_slicer_card),
                        dbc.Col(coronal_slicer_card),
                        dbc.Col(axial_slicer_card),
                    ],
                    style={"margin-top": "1em"},
                ),
                dbc.Row(dbc.Col(controls_card), style={"margin-top": "1em"}),
                dbc.Row(
                    [dbc.Col(table_card), dbc.Col(connectivity_card)],
                    style={"margin-top": "1em"},
                ),
            ],
            fluid=True,
        ),
    ]
)


# Update the nav table values to show the current slicer positions
# Note that we are listening to the "drag_value" attributes here to avoid circularity of the callbacks
@app.callback(
    [Output("x-vox", "value"), Output("y-vox", "value"), Output("z-vox", "value")],
    [
        Input(saggital_slicer.slider.id, "drag_value"),
        Input(coronal_slicer.slider.id, "drag_value"),
        Input(axial_slicer.slider.id, "drag_value"),
        Input("store_graph_position", "data"),
    ],
)
def update_navtable(x_pos, y_pos, z_pos, requested_position):
    ctx = dash.callback_context
    if ctx.triggered[0]["prop_id"] == "store_graph_position.data":
        x_pos, y_pos, z_pos = np.floor(apply_affine(np.linalg.inv(img.affine), requested_position)).astype(int)
    return x_pos, y_pos, z_pos


# Listen for changes to the nav table values and update the slicer position accordingly
@app.callback(
    [
        Output(setpos_store.id, "data"),
        Output("x-world", "value"),
        Output("y-world", "value"),
        Output("z-world", "value"),
        Output("write_selected_position", "children"),
        Output("store_slicer_region", "data"),
    ],
    [Input("x-vox", "value"), Input("y-vox", "value"), Input("z-vox", "value"),],
    State("region_table", "data"),
)
def write_table_values_to_slicer(x_vox, y_vox, z_vox, table_data):
    # Find the overlay value at the current slicer position
    region_number = parcellation[x_vox, y_vox, z_vox]
    # Compute the world positions
    x_world, y_world, z_world = apply_affine(img.affine, (x_vox, y_vox, z_vox))
    if not region_number == 0:
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
        f"Current voxel positions\nx: {x_vox}, y: {y_vox}, z: {z_vox}",
        region_name,
    )


# TODO remove this
@app.callback(
    Output("write_slicer_region", "children"), Input("store_slicer_region", "data")
)
def temp_writer1(string):
    return f"The current slicer region is: {string}"


@app.callback(
    Output("write_graph_region", "children"), Input("store_graph_region", "data")
)
def temp_writer2(string):
    return f"The current graph region is: {string}"


@app.callback(
    Output("write_highlight_region", "children"),
    Input("store_highlight_region", "data"),
)
def temp_writer3(string):
    return f"The current highlight region is: {string}"


# React to clicks on graphs and update the corresponding region
@app.callback(
    [Output("store_graph_region", "data"), Output("store_graph_position", "data")],
    [Input("conn_mat", "clickData"), Input("region_table", "active_cell")],
    [State("region_table", "data"), State("store_graph_region", "data")],
)
def listen_to_graph_clicks(
    heatmap_clickdata, active_table_cell, table_data, previous_region
):
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
    State("store_conn_mat", "data"),
)
def select_thresholding(val, conn_mat):
    min_val = 0
    max_val = 1
    step = 0.05
    ticks = np.arange(min_val, max_val + step, step)

    if val == "absolute":
        # TODO: check if we still need this complex rounding for float precision
        marks = {
            key: val
            for key, val in [
                (int(i), f"{int(i)}") if i.is_integer() else (i, f"{i:.2f}")
                for i in map(lambda x: np.round(x, 2), ticks)
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
    ctx = dash.callback_context
    if ctx.triggered[0]["prop_id"] == "store_graph_region.data":
        return graph_region
    elif ctx.triggered[0]["prop_id"] == "store_slicer_region.data":
        return slicer_region
    else:
        # TODO remove this, 'tis silly
        raise Exception("This is not supposed to happen")


@app.callback(
    Output("conn_mat", "figure"),
    [Input("val_slider", "value"), Input("store_highlight_region", "data")],
    [
        State("store_conn_mat", "data"),
        State("store_region_names", "data"),
        State("threshold-mode", "value"),
    ],
)
def connectivity_heatmap(
    threshold, current_region, conn_mat, region_names, thr_mode
):
    conn_mat = np.array(conn_mat)
    conn_mat_df = hf.thr_conn_mat(conn_mat, threshold, region_names, thr_mode)
    # Redraw the figure
    fig = hf.make_heatmap(conn_mat_df)

    if current_region is not None:
        fig = hf.add_region_shape(fig, current_region)
    return fig


# Draw circos here
@app.callback(
    Output("circos_graph", "children"),
    [Input("val_slider", "value"), Input("store_highlight_region", "data")],
    [
        State("store_conn_mat", "data"),
        State("store_region_names", "data"),
        State("threshold-mode", "value"),
    ],
)
def connectivity_circos(threshold, current_region, conn_mat, region_names, thr_mode):
    conn_mat = np.array(conn_mat)
    conn_mat_df = hf.thr_conn_mat(conn_mat, threshold, region_names, thr_mode)

    layout, highlight, connection = hf.make_circos(
        conn_mat_df, focus_region=current_region
    )

    svg_size = 500
    inner = 200
    outer = 220
    circos = dashbio.Circos(
        size=svg_size,
        layout=layout,
        selectEvent={"0": "hover"},
        config={"innerRadius": inner, "outerRadius": outer},
        tracks=[
            {
                "type": "HIGHLIGHT",
                "data": highlight,
                "config": {
                    "tooltipContent": {"name": "name"},
                    "color": {"name": "color"},
                    "innerRadius": inner,
                    "outerRadius": outer,
                },
            },
            {
                "type": "CHORDS",
                "data": connection,
                "config": {
                    "opacity": 0.5,
                    "tooltipContent": {"name": "name"},
                    "color": {"name": "color"},
                },
            },
        ],
    )
    return circos


# Highlight the table region
@app.callback(
    Output("region_table", "style_data_conditional"),
    Input("store_highlight_region", "data"),
    prevent_initial_call=True,
)
def higlight_row(region_name):
    """
    When a region is selected, highlight the corresponding table row.
    """
    if region_name is not None:
        return [
            {
                "if": {"filter_query": "{Region Name} = '%s'" % region_name},
                "backgroundColor": "#3D9970",
                "color": "white",
            }
        ]
    return None


if __name__ == "__main__":
    app.run_server(debug=True, dev_tools_props_check=False)
