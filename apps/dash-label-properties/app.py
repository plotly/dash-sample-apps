import dash
import dash_table
from dash_table.Format import Format
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import plotly.express as px

import numpy as np
from skimage import io, filters, measure, color, img_as_ubyte
import PIL
import pandas as pd
import matplotlib as mpl


# Set up the app
external_stylesheets = [dbc.themes.BOOTSTRAP, "assets/object_properties_style.css"]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

# Read the image that will be segmented
filename = (
    "https://upload.wikimedia.org/wikipedia/commons/a/ac/Monocyte_no_vacuoles.JPG"
)
img = io.imread(filename, as_gray=True)[:660:2, :800:2]
label_array = measure.label(img < filters.threshold_otsu(img))
current_labels = np.unique(label_array)[np.nonzero(np.unique(label_array))]
# Compute and store properties of the labeled image
prop_names = [
    "label",
    "area",
    "perimeter",
    "eccentricity",
    "euler_number",
    "mean_intensity",
]
prop_table = measure.regionprops_table(
    label_array, intensity_image=img, properties=prop_names
)
table = pd.DataFrame(prop_table)
# Format the Table columns
columns = [
    {"name": label_name, "id": label_name, "selectable": True}
    if precision is None
    else {
        "name": label_name,
        "id": label_name,
        "type": "numeric",
        "format": Format(precision=precision),
        "selectable": True,
    }
    for label_name, precision in zip(prop_names, (None, None, 4, 4, None, 3))
]
# Select the columns that are selected when the app starts
initial_columns = ["label", "area"]

img = img_as_ubyte(color.gray2rgb(img))
img = PIL.Image.fromarray(img)


def image_with_contour(img, active_labels, data_table, active_columns, color_column):
    """
    Returns a greyscale image that is segmented and superimposed with contour traces of
    the segmented regions, color coded by values from a data table.

    Parameters
    ----------
    img : PIL Image object.
    active_labels : list
        the currently visible labels in the datatable
    data_table : pandas.DataFrame
        the currently visible entries of the datatable
    active_columns: list
        the currently selected columns of the datatable
    color_column: str
        name of the datatable column that is used to define the colorscale of the overlay
    """

    # First we get the values from the selected datatable column and use them to define a colormap
    values = np.array(table[color_column].values)
    norm = mpl.colors.Normalize(vmin=values.min(), vmax=values.max())
    cmap = mpl.cm.get_cmap("viridis")

    # Now we convert our background image to a greyscale bytestring that is very small and can be transferred very
    # efficiently over the network. We do not want any hover-information for this image, so we disable it
    fig = px.imshow(img, binary_string=True, binary_backend="jpg",)
    fig.update_traces(hoverinfo="skip", hovertemplate=None)

    # For each region that is visible in the datatable, we compute and draw the filled contour, color it based on
    # the color_column value of this region, and add it to the figure
    # here is an small tutorial of this: https://scikit-image.org/docs/dev/auto_examples/segmentation/plot_regionprops.html#sphx-glr-auto-examples-segmentation-plot-regionprops-py
    for rid, row in data_table.iterrows():
        label = row.label
        value = row[color_column]
        contour = measure.find_contours(
            label_array == label, 0.5, fully_connected="high"
        )[0]
        y, x = contour.T
        # We add the values of the selected datatable columns to the hover information of the current region
        hoverinfo = (
            "<br>".join(
                [
                    # All numbers are passed as floats. If there are no decimals, cast to int for visibility
                    f"{prop_name}: {f'{int(prop_val):d}' if prop_val.is_integer() else f'{prop_val:.3f}'}"
                    if np.issubdtype(type(prop_val), "float")
                    else f"{prop_name}: {prop_val}"
                    for prop_name, prop_val in row[active_columns].iteritems()
                ]
            )
            # remove the trace name. See e.g. https://plotly.com/python/reference/#scatter-hovertemplate
            + " <extra></extra>"
        )
        fig.add_scatter(
            x=x,
            y=y,
            name=label,
            opacity=0.8,
            mode="lines",
            line=dict(color=mpl.colors.rgb2hex(cmap(norm(value))),),
            fill="toself",
            customdata=[label] * len(x),
            showlegend=False,
            hovertemplate=hoverinfo,
            hoveron="points+fills",
        )

    # Finally, because we color our contour traces one by one, we need to manually add a colorscale to explain the
    # mapping of our color_column values to the colormap. This also gets added to the figure
    fig.add_scatter(
        # We only care about the colorscale here, so the x and y values can be empty
        x=[None],
        y=[None],
        mode="markers",
        showlegend=False,
        marker=dict(
            colorscale=[mpl.colors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, 50)],
            showscale=True,
            # The cmin and cmax values here are arbitrary, we just set them to put our value ticks in the right place
            cmin=-5,
            cmax=5,
            colorbar=dict(
                tickvals=[-5, 5],
                ticktext=[f"{np.min(values[values!=0]):.2f}", f"{np.max(values):.2f}",],
                # We want our colorbar to scale with the image when it is resized, so we set them to
                # be a fraction of the total image container
                lenmode="fraction",
                len=0.6,
                thicknessmode="fraction",
                thickness=0.05,
                outlinewidth=1,
                # And finally we give the colorbar a title so the user may know what value the colormap is based on
                title=dict(text=f"<b>{color_column.capitalize()}</b>"),
            ),
        ),
        hoverinfo="none",
    )

    # Remove axis ticks and labels and have the image fill the container
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0, pad=0), template="simple_white")
    fig.update_xaxes(visible=False).update_yaxes(visible=False)

    return fig


# Define Modal
with open("assets/modal.md", "r") as f:
    howto_md = f.read()

modal_overlay = dbc.Modal(
    [
        dbc.ModalBody(html.Div([dcc.Markdown(howto_md)], id="howto-md")),
        dbc.ModalFooter(dbc.Button("Close", id="howto-close", className="howto-bn")),
    ],
    id="modal",
    size="lg",
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
    href="https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-label-properties",
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
                    modal_overlay,
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

# Color selector dropdown
color_drop = dcc.Dropdown(
    id="color-drop-menu",
    options=[
        {"label": col_name.capitalize(), "value": col_name}
        for col_name in table.columns
    ],
    value="label",
)

# Define Cards

image_card = dbc.Card(
    [
        dbc.CardHeader(html.H2("Explore object properties")),
        dbc.CardBody(
            dbc.Row(
                dbc.Col(
                    dcc.Graph(
                        id="graph",
                        figure=image_with_contour(
                            img,
                            current_labels,
                            table,
                            initial_columns,
                            color_column="area",
                        ),
                    ),
                )
            )
        ),
        dbc.CardFooter(
            dbc.Row(
                [
                    dbc.Col(
                        "Use the dropdown menu to select which variable to base the colorscale on:"
                    ),
                    dbc.Col(color_drop),
                    dbc.Toast(
                        [
                            html.P(
                                "In order to use all functions of this app, please select a value "
                                "to compute the colorscale on",
                                className="mb-0",
                            )
                        ],
                        id="auto-toast",
                        header="No colorscale value selected",
                        icon="danger",
                        style={
                            "position": "fixed",
                            "top": 66,
                            "left": 10,
                            "width": 350,
                        },
                    ),
                ],
                align="center",
            ),
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
                            tooltip_header={
                                col: "Select columns with the checkbox to include them in the hover info of the image."
                                for col in table.columns
                            },
                            style_header={
                                "textDecoration": "underline",
                                "textDecorationStyle": "dotted",
                            },
                            tooltip_delay=0,
                            tooltip_duration=None,
                            filter_action="native",
                            row_deletable=True,
                            column_selectable="multi",
                            selected_columns=initial_columns,
                            style_table={"overflowY": "scroll"},
                            fixed_rows={"headers": False, "data": 0},
                            style_cell={"width": "85px"},
                        ),
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
    Output("table-line", "style_data_conditional"),
    [Input("graph", "hoverData")],
    prevent_initial_call=True,
)
def higlight_row(string):
    """
    When hovering hover label, highlight corresponding row in table,
    using label column.
    """
    index = string["points"][0]["customdata"]
    return [
        {
            "if": {"filter_query": "{label} eq %d" % index},
            "backgroundColor": "#3D9970",
            "color": "white",
        }
    ]


@app.callback(
    [
        Output("graph", "figure"),
        Output("row", "children"),
        Output("auto-toast", "is_open"),
    ],
    [
        Input("table-line", "derived_virtual_indices"),
        Input("table-line", "active_cell"),
        Input("table-line", "data"),
        Input("table-line", "selected_columns"),
        Input("color-drop-menu", "value"),
    ],
    [State("row", "children")],
    prevent_initial_call=True,
)
def highlight_filter(
    indices, cell_index, data, active_columns, color_column, previous_row
):
    """
    Updates figure and labels array when a selection is made in the table.

    When a cell is selected (active_cell), highlight this particular label
    with a red outline.

    When the set of filtered labels changes, or when a row is deleted.
    """
    # If no color column is selected, open a popup to ask the user to select one.
    if color_column is None:
        return [dash.no_update, dash.no_update, True]

    _table = pd.DataFrame(data)
    filtered_labels = _table.loc[indices, "label"].values
    filtered_table = _table.query("label in @filtered_labels")
    fig = image_with_contour(
        img, filtered_labels, filtered_table, active_columns, color_column
    )

    if cell_index and cell_index["row"] != previous_row:
        label = filtered_labels[cell_index["row"]]
        mask = (label_array == label).astype(np.float)
        # Add the outline of the selected label as a contour to the figure
        fig.add_contour(
            z=mask,
            contours=dict(coloring="lines"),
            showscale=False,
            line=dict(width=8),
            colorscale="YlOrRd",
            opacity=0.8,
            hoverinfo="skip",
        )
        return [fig, cell_index["row"], False]

    return [fig, previous_row, False]


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
    app.run_server(debug=False)
