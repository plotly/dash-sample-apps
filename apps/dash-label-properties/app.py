import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
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

    fig = px.imshow(img, binary_backend="jpg",)
    fig.add_contour(
        z=labels,
        contours=dict(start=0, end=labels.max() + 1, size=1, coloring=mode),
        line=dict(width=1),
        showscale=False,
        colorscale=custom_viridis,
        opacity=opacity,
    )
    # Remove axis ticks and labels
    fig.update_layout(
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks="",
            showticklabels=False,
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks="",
            showticklabels=False,
        ),
    )

    return fig


# Image to segment
filename = (
    "https://upload.wikimedia.org/wikipedia/commons/a/ac/Monocyte_no_vacuoles.JPG"
)
img = io.imread(filename, as_gray=True)[:660:2, :800:2]
labels = measure.label(img < filters.threshold_otsu(img))

height, width = img.shape
canvas_width = 600
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
columns = [{"name": i, "id": i} for i in list_columns]
columns[2]["format"] = Format(precision=4)
columns[2]["type"] = "numeric"
columns[3]["format"] = Format(precision=4)
columns[3]["type"] = "numeric"
columns[5]["format"] = Format(precision=3)
columns[5]["type"] = "numeric"

data = pd.DataFrame(
    [[getattr(prop, col) for col in list_columns] for prop in props],
    columns=list_columns,
)

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(
    [
        html.Div(
            [
                html.Div(
                    [
                        html.H4("Explore objects properties"),
                        dcc.Graph(
                            id="graph",
                            figure=image_with_contour(img, labels, mode=None),
                        ),
                    ],
                    className="six columns",
                ),
                html.Div(
                    [
                        html.Img(
                            src="https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png",
                            width="30px",
                        ),
                        html.A(
                            id="gh-link",
                            children=["View on GitHub"],
                            href="http://github.com/plotly/canvas-portal/"
                            "blob/master/apps/object-properties/app.py",
                            style={
                                "color": "black",
                                "border": "solid 1px black",
                                "float": "left",
                            },
                        ),
                        dash_table.DataTable(
                            id="table-line",
                            columns=columns,
                            data=data.to_dict("records"),
                            filter_action="native",
                            row_deletable=True,
                            style_table={"overflowY": "scroll"},
                            fixed_rows={"headers": False, "data": 0},
                            style_cell={"width": "85px"},
                        ),
                        dcc.Store(id="cache", data=labels),
                        html.Div(id="row", hidden=True, children=None),
                    ],
                    className="six columns",
                ),
            ],
            className="row",
        ),
        html.H4("How to use this app (see below)"),
        dcc.Markdown(
            dedent(
                """
    Hover over objects to highlight their properties in the table,
    select cell in table to highlight object in image, or
    filter objects in the table to display a subset of objects.

    Learn more about [DataTable filtering syntax](https://dash.plot.ly/datatable/filtering)
    for selecting ranges of properties.
    """
            )
        ),
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


if __name__ == "__main__":
    app.run_server(debug=True)
