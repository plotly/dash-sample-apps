import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import dash_table
from dash_table.Format import Format
import plotly.express as px
import plotly.graph_objs as go

import numpy as np
from skimage import io, filters, measure
import pandas as pd

import PIL
from skimage import color, img_as_ubyte
from plotly import colors
from textwrap import dedent


external_stylesheets = [dbc.themes.BOOTSTRAP]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

# Read the image that will be segmented
filename = (
    "https://upload.wikimedia.org/wikipedia/commons/a/ac/Monocyte_no_vacuoles.JPG"
)
img = io.imread(filename, as_gray=True)[:660:2, :800:2]
labels = measure.label(img < filters.threshold_otsu(img))
props = measure.regionprops(labels, img)

# Define table columns
list_columns = [
    "label",
    "area",
    "perimeter",
    "eccentricity",
    "euler_number",
    "mean_intensity",
]
columns = [
    {"name": label_name, "id": label_name}
    if precision is None
    else {
        "name": label_name,
        "id": label_name,
        "type": "numeric",
        "format": Format(precision=precision),
    }
    for label_name, precision in zip(list_columns, (None, None, 4, 4, None, 3))
]
table = pd.DataFrame(
    [[getattr(prop, col) for col in list_columns] for prop in props],
    columns=list_columns,
)


def image_with_contour(img, labels, mode="lines", shape=None):
    """
    Figure with contour plot of labels superimposed on background image.

    Parameters
    ----------

    img : URL, dataURI or ndarray
        Background image. If a numpy array, it is transformed into a PIL
        Image object.
    labels : 2D ndarray
        Contours are the isolines of labels.
    shape: tuple, optional
        Shape of the arrays, to be provided if ``img`` is not a numpy array.
    """
    try:
        sh_y, sh_x = shape if shape is not None else img.shape
    except AttributeError:
        print(
            """the shape of the image must be provided with the
                 ``shape`` parameter if ``img`` is not a numpy array"""
        )
    if type(img) == np.ndarray:
        img = img_as_ubyte(color.gray2rgb(img))
        img = PIL.Image.fromarray(img)
    labels = labels.astype(np.float)
    custom_viridis = colors.PLOTLY_SCALES["Viridis"]
    custom_viridis.insert(0, [0, "#FFFFFF"])
    custom_viridis[1][0] = 1.0e-4
    # Contour plot of segmentation
    print("mode is", mode)
    opacity = 0.4 if mode is None else 1

    fig = px.imshow(img, binary_backend="jpg")
    fig.add_contour(
        z=labels,
        contours=dict(start=0, end=labels.max() + 1, size=1, coloring=mode),
        line=dict(width=1),
        showscale=False,
        colorscale=custom_viridis,
        opacity=opacity,
    )
    # Remove axis ticks and labels
    fig.update_xaxes(showticklabels=False).update_yaxes(showticklabels=False)

    return fig


# Define Modal
howto = """
    Hover over objects to highlight their properties in the table,
    select cell in table to highlight object in image, or
    filter objects in the table to display a subset of objects.

    Learn more about [DataTable filtering syntax](https://dash.plot.ly/datatable/filtering)
    for selecting ranges of properties.
    """
modal = dbc.Modal(
    [
        dbc.ModalHeader("How to use this app"),
        dbc.ModalBody(
            dbc.Row(
                dbc.Col(
                    [
                        dcc.Markdown(howto, id="howto-md"),
                        html.Img(
                            id="help",
                            src="assets/properties.gif",
                            width="80%",
                            style={
                                "border": "2px solid black",
                                "display": "block",
                                "margin-left": "auto",
                                "margin-right": "auto",
                            },
                        ),
                    ]
                )
            )
        ),
        dbc.ModalFooter(dbc.Button("Close", id="howto-close", className="howto-bn",)),
    ],
    id="modal",
    size="lg",
    style={"font-size": "small"},
)

# Buttons
button_gh = dbc.Button(
    "Learn more",
    id="howto-open",
    outline=True,
    color="secondary",
    # Turn off lowercase transformation for class .button in stylesheet
    style={"textTransform": "none"},
)

button_howto = dbc.Button(
    "View Code on github",
    outline=True,
    color="primary",
    href="https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-image-annotation",
    id="gh-link",
    style={"text-transform": "none"},
)

# Define Header Layout
header = dbc.Navbar(
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
                            href="https://plotly.com/dash/",
                        )
                    ),
                    dbc.Col(dbc.NavbarBrand("Object Properties App")),
                    modal,
                ],
                align="center",
            ),
            dbc.Row(
                dbc.Col(
                    [
                        dbc.NavbarToggler(id="navbar-toggler"),
                        dbc.Collapse(
                            dbc.Nav(
                                [dbc.NavItem(button_howto), dbc.NavItem(button_gh)],
                                className="ml-auto",
                                navbar=True,
                            ),
                            id="navbar-collapse",
                            navbar=True,
                        ),
                    ]
                ),
                align="center",
            ),
        ],
        fluid=True,
    ),
    color="dark",
    dark=True,
)

# Define Cards
image_card = dbc.Card(
    [
        dbc.CardHeader(html.H2("Explore object properties")),
        dbc.CardBody(
            dbc.Row(
                dbc.Col(
                    dcc.Graph(
                        id="graph", figure=image_with_contour(img, labels, mode=None),
                    )
                )
            )
        ),
    ]
)

table_card = dbc.Card(
    [
        dbc.CardHeader(html.H2("Data Table")),
        dbc.CardBody(
            dbc.Row(
                dbc.Col(
                    [
                        dash_table.DataTable(
                            id="table-line",
                            columns=columns,
                            data=table.to_dict("records"),
                            filter_action="native",
                            row_deletable=True,
                            style_table={"overflowY": "scroll"},
                            fixed_rows={"headers": False, "data": 0},
                            style_cell={"width": "85px"},
                        ),
                        dcc.Store(id="cache", data=labels),
                        html.Div(id="row", hidden=True, children=None),
                    ]
                )
            )
        ),
    ]
)

app.layout = html.Div(
    [
        header,
        dbc.Container(
            [dbc.Row([dbc.Col(image_card, md=5), dbc.Col(table_card, md=7)])],
            fluid=True,
        ),
    ]
)


@app.callback(
    Output("table-line", "style_data_conditional"), [Input("graph", "hoverData")]
)
def higlight_row(string):
    """
    When hovering hover label, highlight corresponding row in table,
    using label column.
    """
    if not dash.callback_context.triggered:
        # Nothing has happened yet
        return []
    index = string["points"][0]["z"]
    return [
        {
            "if": {"filter_query": "{label} eq %d" % index},
            "backgroundColor": "#3D9970",
            "color": "white",
        }
    ]


@app.callback(
    [Output("graph", "figure"), Output("cache", "data"), Output("row", "children")],
    [
        Input("table-line", "derived_virtual_indices"),
        Input("table-line", "active_cell"),
        Input("table-line", "data"),
    ],
    [State("cache", "data"), State("row", "children")],
)
def highlight_filter(indices, cell_index, data, current_labels, previous_row):
    """
    Updates figure and labels array when a selection is made in the table.

    When a cell is selected (active_cell), highlight this particular label
    with a red outline.

    When the set of filtered labels changes, or when a row is deleted. 
    """
    if cell_index and cell_index["row"] != previous_row:
        current_labels = np.asanyarray(current_labels)
        label = indices[cell_index["row"]] + 1
        mask = (labels == label).astype(np.float)
        fig = image_with_contour(img, current_labels, mode=None)
        # Add the outline of the selected label as a contour to the figure
        fig.add_contour(
            z=mask,
            contours=dict(coloring="lines"),
            showscale=False,
            line=dict(width=6),
            colorscale="YlOrRd",
            opacity=0.8,
            hoverinfo="skip",
        )
        return [fig, current_labels, cell_index["row"]]
    if not dash.callback_context.triggered:
        # Nothing has happened yet
        new_labels = labels
    else:
        filtered_labels = np.array(
            pd.DataFrame(data).lookup(np.array(indices), ["label",] * len(indices))
        )
        mask = np.in1d(labels.ravel(), filtered_labels).reshape(labels.shape)
        new_labels = np.copy(labels)
        new_labels *= mask
    fig = image_with_contour(img, new_labels, mode=None)
    return [fig, new_labels, previous_row]


@app.callback(
    Output("modal", "is_open"),
    [Input("howto-open", "n_clicks"), Input("howto-close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


# we use a callback to toggle the collapse on small screens
@app.callback(
    Output("navbar-collapse", "is_open"),
    [Input("navbar-toggler", "n_clicks")],
    [State("navbar-collapse", "is_open")],
)
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


if __name__ == "__main__":
    app.run_server(debug=True)
