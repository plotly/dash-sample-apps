import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_table
import plotly.express as px
import re
import time
from skimage import io

DEBUG = True

NUM_ATYPES = 15
DEFAULT_FIG_MODE = "layout"
annotation_colormap = px.colors.qualitative.Light24
annotation_types = [
    "tree",
    "building",
    "sky",
    "road",
    "sidewalk",
    "car",
    "pedestrian",
    "cyclist",
    "stop sign",
    "parking sign",
    "traffic light",
    "lamp post",
    "star",  # e.g., sun or moon as to not confuse them with artificial lighting
]
DEFAULT_ATYPE = annotation_types[0]

# prepare bijective type<->color mapping
typ_col_pairs = [
    (t, annotation_colormap[n % len(annotation_colormap)])
    for n, t in enumerate(annotation_types)
]
# types to colors
color_dict = {}
# colors to types
type_dict = {}
for typ, col in typ_col_pairs:
    color_dict[typ] = col
    type_dict[col] = typ

options = list(color_dict.keys())
columns = ["Type", "X0", "Y0", "X1", "Y1"]
# Open the readme for use in the context info
with open("README.md", "r") as f:
    readme = f.readlines()


def debug_print(*args):
    if DEBUG:
        print(*args)


def coord_to_tab_column(coord):
    return coord.upper()


def time_passed(start=0):
    return round(time.mktime(time.localtime())) - start


def format_float(f):
    return "%.2f" % (float(f),)


def shape_to_table_row(sh):
    return {
        "Type": type_dict[sh["line"]["color"]],
        "X0": format_float(sh["x0"]),
        "Y0": format_float(sh["y0"]),
        "X1": format_float(sh["x1"]),
        "Y1": format_float(sh["y1"]),
    }


def default_table_row():
    return {
        "Type": DEFAULT_ATYPE,
        "X0": format_float(10),
        "Y0": format_float(10),
        "X1": format_float(20),
        "Y1": format_float(20),
    }


def table_row_to_shape(tr):
    return {
        "editable": True,
        "xref": "x",
        "yref": "y",
        "layer": "above",
        "opacity": 1,
        "line": {"color": color_dict[tr["Type"]], "width": 4, "dash": "solid"},
        "fillcolor": "rgba(0, 0, 0, 0)",
        "fillrule": "evenodd",
        "type": "rect",
        "x0": tr["X0"],
        "y0": tr["Y0"],
        "x1": tr["X1"],
        "y1": tr["Y1"],
    }


def shape_cmp(s0, s1):
    """ Compare two shapes """
    return (
        (s0["x0"] == s1["x0"])
        and (s0["x1"] == s1["x1"])
        and (s0["y0"] == s1["y0"])
        and (s0["y1"] == s1["y1"])
        and (s0["line"]["color"] == s1["line"]["color"])
    )


def shape_in(se):
    """ check if a shape is in list (done this way to use custom compare) """
    return lambda s: any(shape_cmp(s, s_) for s_ in se)


def index_of_shape(shapes, shape):
    for i, shapes_item in enumerate(shapes):
        if shape_cmp(shapes_item, shape):
            return i
    raise ValueError  # not found


def annotations_table_shape_resize(annotations_table_data, fig_data):
    """
    Extract the shape that was resized (its index) and store the resized
    coordinates.
    """
    debug_print("fig_data", fig_data)
    debug_print("table_data", annotations_table_data)
    for key, val in fig_data.items():
        shape_nb, coord = key.split(".")
        # shape_nb is for example 'shapes[2].x0': this extracts the number
        shape_nb = shape_nb.split(".")[0].split("[")[-1].split("]")[0]
        # this should correspond to the same row in the data table
        # we have to format the float here because this is exactly the entry in
        # the table
        annotations_table_data[int(shape_nb)][
            coord_to_tab_column(coord)
        ] = format_float(fig_data[key])
        # (no need to compute a time stamp, that is done for any change in the
        # table values, so will be done later)
    return annotations_table_data


def shape_data_remove_timestamp(shape):
    """
    go.Figure complains if we include the 'timestamp' key when updating the
    figure
    """
    new_shape = dict()
    for k in shape.keys() - set(["timestamp"]):
        new_shape[k] = shape[k]
    return new_shape


external_stylesheets = [
    dbc.themes.BOOTSTRAP,
    "assets/style.css",
    "assets/app_bounding_box_style.css",
]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

filelist = [
    "assets/driving.jpg",
    "assets/professional-transport-autos-bridge-traffic-road-rush-hour.jpg",
    "assets/rocket.jpg",
]

server = app.server

fig = px.imshow(io.imread(filelist[0]), binary_backend="jpg")
fig.update_layout(
    newshape_line_color=color_dict[DEFAULT_ATYPE],
    margin=dict(l=0, r=0, b=0, t=0, pad=4),
    dragmode="drawrect",
)
app.layout = html.Div(
    id="main",
    children=[
        # Banner display
        html.Div(
            id="banner",
            children=[
                dbc.Row(
                    [
                        dbc.Col(
                            html.Img(
                                id="logo", src=app.get_asset_url("dash-logo-new.png")
                            ),
                            width=2,
                            align="center",
                        ),
                        dbc.Col(
                            html.H1("Bounding Box Classification App", id="title"),
                            width=4,
                            align="center",
                        ),
                        dbc.Col(
                            [
                                dbc.Button(
                                    "Readme",
                                    id="readme-open",
                                    color="secondary",
                                    style={"margin": "5px"},
                                ),
                                dbc.Modal(
                                    [
                                        dbc.ModalBody(
                                            html.Div(
                                                [dcc.Markdown(readme, id="readme-md")]
                                            )
                                        ),
                                        dbc.ModalFooter(
                                            dbc.Button(
                                                "Close",
                                                id="readme-close",
                                                className="ml-auto",
                                            )
                                        ),
                                    ],
                                    id="modal",
                                    size="md",
                                    style={"font-size": "small"},
                                ),
                                dbc.Button(
                                    "View Code on github",
                                    outline=True,
                                    color="primary",
                                    href="https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-image-annotation",
                                    id="gh-link",
                                ),
                            ],
                            width=2,
                            align="center",
                        ),
                    ]
                ),
            ],
            className="twelve columns",
        ),
        # Main body
        html.Div(
            id="app-container",
            children=[
                # Graph
                dcc.Graph(
                    id="graph",
                    figure=fig,
                    config={"modeBarButtonsToAdd": ["drawrect", "eraseshape"]},
                ),
            ],
            className="seven columns",
        ),
        # Sidebar
        html.Div(
            id="sidebar",
            children=[
                html.H2("Annotations", id="firsth2"),
                html.Div(
                    id="table-container",
                    children=[
                        # Timestamp table
                        dash_table.DataTable(
                            id="timestamp-table",
                            columns=[{"name": "Timestamp", "id": "Timestamp"}],
                            style_data={"height": 40},
                        ),
                        # Data table
                        dash_table.DataTable(
                            id="annotations-table",
                            columns=[
                                dict(
                                    name=n,
                                    id=n,
                                    presentation=(
                                        "dropdown" if n == "Type" else "input"
                                    ),
                                )
                                for n in columns
                            ],
                            editable=True,
                            style_data={"height": 40},
                            dropdown={
                                "Type": {
                                    "options": [
                                        {"label": o, "value": o}
                                        for o in annotation_types
                                    ],
                                    "clearable": False,
                                }
                            },
                            style_cell_conditional=[
                                {"if": {"column_id": "Type"}, "textAlign": "left"}
                            ],
                        ),
                    ],
                ),
                dcc.Store(id="graph-copy", data=fig),
                dcc.Store(
                    id="annotations-store",
                    data=dict(
                        **{filename: {"shapes": []} for filename in filelist},
                        **{"starttime": time_passed()}
                    ),
                ),
                dcc.Store(id="image_files", data={"files": filelist, "current": 0}),
                html.H2("Type of annotation"),
                dcc.Dropdown(
                    id="annotation-type-dropdown",
                    options=[{"label": t, "value": t} for t in annotation_types],
                    value=DEFAULT_ATYPE,
                    clearable=False,
                ),
                html.H2("Choose image"),
                html.Button("Previous", id="previous", className="button"),
                html.Button("Next", id="next", className="button"),
                # We use this pattern because we want to be able to download the
                # annotations by clicking on a button
                html.A(
                    id="download",
                    download="annotations.json",
                    # make invisble, we just want it to click on it
                    style={"display": "none"},
                ),
                html.Button(
                    "Download annotations", id="download-button", className="button"
                ),
                html.Div(id="dummy", style={"display": "none"}),
            ],
            className="five columns",
        ),
    ],
    className="twelve columns",
)


@app.callback(
    [Output("annotations-table", "data"), Output("image_files", "data")],
    [
        Input("previous", "n_clicks"),
        Input("next", "n_clicks"),
        Input("graph", "relayoutData"),
    ],
    [
        State("annotations-table", "data"),
        State("image_files", "data"),
        State("annotations-store", "data"),
        State("annotation-type-dropdown", "value"),
    ],
)
def modify_table_entries(
    previous_n_clicks,
    next_n_clicks,
    graph_relayoutData,
    annotations_table_data,
    image_files_data,
    annotations_store_data,
    annotation_type,
):
    cbcontext = [p["prop_id"] for p in dash.callback_context.triggered][0]
    if cbcontext == "graph.relayoutData":
        debug_print("graph_relayoutData:", graph_relayoutData)
        debug_print("annotations_table_data before:", annotations_table_data)
        if "shapes" in graph_relayoutData.keys():
            # this means all the shapes have been passed to this function via
            # graph_relayoutData, so we store them
            annotations_table_data = [
                shape_to_table_row(sh) for sh in graph_relayoutData["shapes"]
            ]
        elif re.match("shapes\[[0-9]+\].x0", list(graph_relayoutData.keys())[0]):
            # this means a shape was updated (e.g., by clicking and dragging its
            # vertices), so we just update the specific shape
            annotations_table_data = annotations_table_shape_resize(
                annotations_table_data, graph_relayoutData
            )
        if annotations_table_data is None:
            return dash.no_update
        else:
            debug_print("annotations_table_data after:", annotations_table_data)
            return (annotations_table_data, image_files_data)
    image_index_change = 0
    if cbcontext == "previous.n_clicks":
        image_index_change = -1
    if cbcontext == "next.n_clicks":
        image_index_change = 1
    image_files_data["current"] += image_index_change
    image_files_data["current"] %= len(image_files_data["files"])
    if image_index_change != 0:
        # image changed, update annotations_table_data with new data
        annotations_table_data = []
        filename = image_files_data["files"][image_files_data["current"]]
        debug_print(annotations_store_data[filename])
        for sh in annotations_store_data[filename]["shapes"]:
            annotations_table_data.append(shape_to_table_row(sh))
        return (annotations_table_data, image_files_data)
    else:
        return dash.no_update


@app.callback(
    [
        Output("graph", "figure"),
        Output("annotations-store", "data"),
        Output("timestamp-table", "data"),
    ],
    [Input("annotations-table", "data"), Input("annotation-type-dropdown", "value")],
    [State("image_files", "data"), State("annotations-store", "data")],
)
def send_figure_to_graph(
    annotations_table_data, annotation_type, image_files_data, annotations_store
):
    if annotations_table_data is not None:
        filename = image_files_data["files"][image_files_data["current"]]
        # convert table rows to those understood by fig.update_layout
        fig_shapes = [table_row_to_shape(sh) for sh in annotations_table_data]
        debug_print("fig_shapes:", fig_shapes)
        debug_print(
            "annotations_store[%s]['shapes']:" % (filename,),
            annotations_store[filename]["shapes"],
        )
        # find the shapes that are new
        new_shapes_i = []
        old_shapes_i = []
        for i, sh in enumerate(fig_shapes):
            if not shape_in(annotations_store[filename]["shapes"])(sh):
                new_shapes_i.append(i)
            else:
                old_shapes_i.append(i)
        # add timestamps to the new shapes
        for i in new_shapes_i:
            fig_shapes[i]["timestamp"] = time_passed(annotations_store["starttime"])
        # find the old shapes and look up their timestamps
        for i in old_shapes_i:
            old_shape_i = index_of_shape(
                annotations_store[filename]["shapes"], fig_shapes[i]
            )
            fig_shapes[i]["timestamp"] = annotations_store[filename]["shapes"][
                old_shape_i
            ]["timestamp"]
        shapes = fig_shapes
        debug_print("shapes:", shapes)
        fig = px.imshow(io.imread(filename), binary_backend="jpg")
        fig.update_layout(
            shapes=[shape_data_remove_timestamp(sh) for sh in shapes],
            # reduce space between image and graph edges
            newshape_line_color=color_dict[annotation_type],
            margin=dict(l=0, r=0, b=0, t=0, pad=4),
            dragmode="drawrect",
        )
        annotations_store[filename]["shapes"] = shapes
        return (fig, annotations_store, [{"Timestamp": s["timestamp"]} for s in shapes])
    return dash.no_update


@app.callback(
    Output("modal", "is_open"),
    [Input("readme-open", "n_clicks"), Input("readme-close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


# set the download url to the contents of the annotations-store (so they can be
# downloaded from the browser's memory)
app.clientside_callback(
    """
function(the_store_data) {
    let s = JSON.stringify(the_store_data);
    let b = new Blob([s],{type: 'text/plain'});
    let url = URL.createObjectURL(b);
    return url;
}
""",
    Output("download", "href"),
    [Input("annotations-store", "data")],
)

# click on download link via button
app.clientside_callback(
    """
function(download_button_n_clicks)
{
    let download_a=document.getElementById("download");
    download_a.click();
    return '';
}
""",
    Output("dummy", "children"),
    [Input("download-button", "n_clicks")],
)


if __name__ == "__main__":
    app.run_server()
