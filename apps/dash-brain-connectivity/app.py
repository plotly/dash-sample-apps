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
from plotly import graph_objects as go

from nibabel.affines import apply_affine
import pathlib
import plotly

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

colormap = px.colors.qualitative.Alphabet
saggital_slicer = VolumeSlicer(app, vol, axis=0, scene_id="brain")
coronal_slicer = VolumeSlicer(app, vol, axis=1, scene_id="brain")
axial_slicer = VolumeSlicer(app, vol, axis=2, scene_id="brain")

# Update all slicers in order
for slicer in [saggital_slicer, coronal_slicer, axial_slicer]:
    slicer.graph.figure.update_layout(dragmode=False, plot_bgcolor="rgb(0, 0, 0)")
    slicer.overlay_data.data = slicer.create_overlay_data(parcellation, colormap)

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
    [
        dbc.CardHeader("Saggital View"),
        dbc.CardBody(
            [
                saggital_slicer.graph,
                html.Div(saggital_slicer.slider, style={"display": "none"}),
                *saggital_slicer.stores,
            ]
        ),
    ],
    color="dark",
    inverse=True,
)

coronal_slicer_card = dbc.Card(
    [
        dbc.CardHeader("Coronal View"),
        dbc.CardBody(
            [
                coronal_slicer.graph,
                html.Div(coronal_slicer.slider, style={"display": "none"}),
                *coronal_slicer.stores,
            ]
        ),
    ],
    color="dark",
    inverse=True,
)

axial_slicer_card = dbc.Card(
    [
        dbc.CardHeader("Axial View"),
        dbc.CardBody(
            [
                axial_slicer.graph,
                html.Div(axial_slicer.slider, style={"display": "none"}),
                *axial_slicer.stores,
            ]
        ),
    ],
    color="dark",
    inverse=True,
)

controls_card = dbc.Card(
    [
        dbc.CardHeader("Controls"),
        dbc.CardBody([html.Div(id="markdown")]),
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

heatmap_tab = dbc.Card(dbc.CardBody([dcc.Graph(id="conn_mat"),]))

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
        dcc.Store(id="store_selected_region"),
        dcc.Store(id="store_selected_brain_location"),
        dcc.Store(id="store_conn_mat", data=conn_mat),
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


# Update the currently selected brain region
@app.callback(
    [
        Output("store_selected_region", "data"),
        Output("store_selected_brain_location", "data"),
    ],
    [
        Input("conn_mat", "clickData"),
        Input("region_table", "active_cell"),
        Input(saggital_slicer.graph.id, "clickData"),
        Input(coronal_slicer.graph.id, "clickData"),
        Input(axial_slicer.graph.id, "clickData")
    ],
    [
        State({"scene": axial_slicer.scene_id, "context": ALL, "name": "state"}, "data"),
        State("store_selected_region", "data"),
        State("region_table", "data"),
    ],
)
def update_selected_region(
    click_data,
    current_table_cell,
    saggital_trigger,
    coronal_trigger,
    axial_trigger,
    brain_index,
    current_region,
    table_data,
):
    # Understand what has happened
    ctx = dash.callback_context
    clicked_region = None
    selected_brain_location = dash.no_update
    table = pd.DataFrame(table_data)
    # See if the heatmap was clicked
    if ctx.triggered[0]["prop_id"] == "conn_mat.clickData":
        last_region_selector = None
        clicked_region = click_data["points"][0]["y"]
        corresponding_row = table.query("`Region Name` == @clicked_region")
        selected_brain_location = tuple(corresponding_row[["X", "Y", "Z"]].values[0])

    # See if the table row was clicked
    elif ctx.triggered[0]["prop_id"] == "region_table.active_cell":
        last_region_selector = None
        table_row = table.iloc[current_table_cell["row"]]
        clicked_region = table_row["Region Name"]
        selected_brain_location = tuple(table_row[["X", "Y", "Z"]])

    elif "slicer" in ctx.triggered[0]["prop_id"]:
        if brain_index is None or all(map(lambda x: x is None, brain_index)):
            return None, selected_brain_location

        x_pos = brain_index[0]["index"]
        y_pos = brain_index[1]["index"]
        z_pos = brain_index[2]["index"]
        clicked_region_number = parcellation[x_pos, y_pos, z_pos]

        # TODO: do something smarter with the fact that we have selected background
        if clicked_region_number == 0:
            return None, selected_brain_location

        corresponding_row = table.query("`Region ID` == @clicked_region_number")
        clicked_region = corresponding_row["Region Name"].values[0]
        # selected_brain_location = tuple(corresponding_row[["X", "Y", "Z"]].values[0])

    if clicked_region == current_region:
        return None, selected_brain_location
    else:
        return clicked_region, selected_brain_location


@app.callback(
    Output(setpos_store.id, "data"),
    Input("store_selected_brain_location", "data"),
    prevent_inital_call=True,
)
def focus_brain_region(selected_brain_index):
    # Extract the MNI coordinates
    if selected_brain_index is None:
        return dash.no_update  # , dash.no_update, dash.no_update
    x, y, z = selected_brain_index
    # Convert them to voxel coordinate
    x_vox, y_vox, z_vox = np.floor(
        apply_affine(np.linalg.inv(parc_img_res.affine), (x, y, z))
    ).astype(int)
    return x_vox, y_vox, z_vox


# Draw the connectivity matrix given the current threshold
# also: redraw it quickly if only a click event has occured
@app.callback(
    [Output("conn_mat", "figure"), Output("markdown", "children"),],
    [Input("val_slider", "value"), Input("store_selected_region", "data")],
    [
        State("conn_mat", "figure"),
        State("store_conn_mat", "data"),
        State("store_region_names", "data"),
        State("threshold-mode", "value"),
    ],
)
def connectivity_heatmap(
    threshold, current_region, fig, conn_mat, region_names, thr_mode
):
    conn_mat = np.array(conn_mat)
    conn_mat_df = hf.thr_conn_mat(conn_mat, threshold, region_names, thr_mode)

    # Understand what has happened
    ctx = dash.callback_context
    if not ctx.triggered:
        # Nothing has happened yet, this is the setup
        # Draw the figure and move on
        fig = hf.make_heatmap(conn_mat_df)
        return fig, "nothing happened yet"

    # If we start because the threshold has changed
    elif ctx.triggered[0]["prop_id"] == "store_selected_region.data":
        fig = go.Figure(fig)
        fig.layout.shapes = []
        # The selected region has been updated by some event
        # We either draw a new shape or we draw only the figure
        if current_region is None:
            return fig, f"removed {current_region}"
        fig = hf.add_region_shape(fig, current_region)
        return fig, f"{current_region}"

    elif ctx.triggered[0]["prop_id"] == "val_slider.value":
        # The threshold has changed
        # Redraw the figure
        fig = hf.make_heatmap(conn_mat_df)
        if current_region is not None:
            # We also want to draw the current region shape on the heatmap
            fig = hf.add_region_shape(fig, current_region)
        return fig, f"{threshold} and {current_region}"

    else:
        raise Exception(f"You shouldn't have come here: {ctx.triggered[0]['prop_id']}")


# Draw circos here
@app.callback(
    Output("circos_graph", "children"),
    [Input("val_slider", "value"), Input("store_selected_region", "data")],
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
    Input("store_selected_region", "data"),
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
    app.run_server(debug=True, dev_tools_props_check=False, port=8052)
