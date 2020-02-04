import os
import json
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_colorscales as dcs
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from mni import create_mesh_data, default_colorscale


app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

server = app.server

GITHUB_LINK = os.environ.get(
    "GITHUB_LINK",
    "https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-brain-viewer",
)

default_colorscale_index = [ea[1] for ea in default_colorscale]

axis_template = {
    "showbackground": True,
    "backgroundcolor": "#141414",
    "gridcolor": "rgb(255, 255, 255)",
    "zerolinecolor": "rgb(255, 255, 255)",
}

plot_layout = {
    "title": "",
    "margin": {"t": 0, "b": 0, "l": 0, "r": 0},
    "font": {"size": 12, "color": "white"},
    "showlegend": False,
    "plot_bgcolor": "#141414",
    "paper_bgcolor": "#141414",
    "scene": {
        "xaxis": axis_template,
        "yaxis": axis_template,
        "zaxis": axis_template,
        "aspectratio": {"x": 1, "y": 1.2, "z": 1},
        "camera": {"eye": {"x": 1.25, "y": 1.25, "z": 1.25}},
        "annotations": [],
    },
}

app.layout = html.Div(
    [
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.Img(
                                            src=app.get_asset_url("dash-logo.png")
                                        ),
                                        html.H4("MRI Reconstruction"),
                                    ],
                                    className="header__title",
                                ),
                                html.Div(
                                    [
                                        html.P(
                                            "Click on the brain to add an annotation. Drag the black corners of the graph to rotate."
                                        )
                                    ],
                                    className="header__info pb-20",
                                ),
                                html.Div(
                                    [
                                        html.A(
                                            "View on GitHub",
                                            href=GITHUB_LINK,
                                            target="_blank",
                                        )
                                    ],
                                    className="header__button",
                                ),
                            ],
                            className="header pb-20",
                        ),
                        html.Div(
                            [
                                dcc.Graph(
                                    id="brain-graph",
                                    figure={
                                        "data": create_mesh_data("human_atlas"),
                                        "layout": plot_layout,
                                    },
                                    config={"editable": True, "scrollZoom": False},
                                )
                            ],
                            className="graph__container",
                        ),
                    ],
                    className="container",
                )
            ],
            className="two-thirds column app__left__section",
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.P(
                                    "Click colorscale to change", className="subheader"
                                ),
                                dcs.DashColorscales(
                                    id="colorscale-picker",
                                    colorscale=default_colorscale_index,
                                ),
                            ]
                        )
                    ],
                    className="colorscale pb-20",
                ),
                html.Div(
                    [
                        html.P("Select option", className="subheader"),
                        dcc.RadioItems(
                            options=[
                                {"label": "Brain Atlas", "value": "human_atlas"},
                                {"label": "Cortical Thickness", "value": "human"},
                                {"label": "Mouse Brain", "value": "mouse"},
                            ],
                            value="human_atlas",
                            id="radio-options",
                            labelClassName="label__option",
                            inputClassName="input__option",
                        ),
                    ],
                    className="pb-20",
                ),
                html.Div(
                    [
                        html.Span("Click data", className="subheader"),
                        html.Span("  |  "),
                        html.Span(
                            "Click on points in the graph.", className="small-text"
                        ),
                        dcc.Loading(
                            html.Pre(id="click-data", className="info__container"),
                            type="dot",
                        ),
                    ],
                    className="pb-20",
                ),
                html.Div(
                    [
                        html.Span("Relayout data", className="subheader"),
                        html.Span("  |  "),
                        html.Span(
                            "Drag the graph corners to rotate it.",
                            className="small-text",
                        ),
                        dcc.Loading(
                            html.Pre(id="relayout-data", className="info__container"),
                            type="dot",
                        ),
                    ],
                    className="pb-20",
                ),
                html.Div(
                    [
                        html.P(
                            [
                                "Dash/Python code on ",
                                html.A(
                                    children="GitHub.",
                                    target="_blank",
                                    href=GITHUB_LINK,
                                    className="red-ish",
                                ),
                            ]
                        ),
                        html.P(
                            [
                                "Brain data from Mcgill's ACE Lab ",
                                html.A(
                                    children="Surface Viewer.",
                                    target="_blank",
                                    href="https://brainbrowser.cbrain.mcgill.ca/surface-viewer#ct",
                                    className="red-ish",
                                ),
                            ]
                        ),
                    ]
                ),
            ],
            className="one-third column app__right__section",
        ),
        dcc.Store(id="annotation_storage"),
    ]
)


def add_marker(x, y, z):
    """ Create a plotly marker dict. """

    return {
        "x": [x],
        "y": [y],
        "z": [z],
        "mode": "markers",
        "marker": {"size": 25, "line": {"width": 3}},
        "name": "Marker",
        "type": "scatter3d",
        "text": ["Click point to remove annotation"],
    }


def add_annotation(x, y, z):
    """ Create plotly annotation dict. """

    return {
        "x": x,
        "y": y,
        "z": z,
        "font": {"color": "black"},
        "bgcolor": "white",
        "borderpad": 5,
        "bordercolor": "black",
        "borderwidth": 1,
        "captureevents": True,
        "ay": -100,
        "arrowcolor": "white",
        "arrowwidth": 2,
        "arrowhead": 0,
        "text": "Click here to annotate<br>(Click point to remove)",
    }


def marker_in_points(points, marker):
    """
    Checks if the marker is in the list of points.
    
    :params points: a list of dict that contains x, y, z
    :params marker: a dict that contains x, y, z
    :returns: index of the matching marker in list
    """

    for index, point in enumerate(points):
        if (
            point["x"] == marker["x"]
            and point["y"] == marker["y"]
            and point["z"] == marker["z"]
        ):
            return index
    return None


@app.callback(
    Output("brain-graph", "figure"),
    [
        Input("brain-graph", "clickData"),
        Input("radio-options", "value"),
        Input("colorscale-picker", "colorscale"),
    ],
    [State("brain-graph", "figure"), State("annotation_storage", "data")],
)
def brain_graph_handler(click_data, val, colorscale, figure, current_anno):
    """ Listener on colorscale, option picker, and graph on click to update the graph. """

    # new option select
    if figure["data"][0]["name"] != val:
        figure["data"] = create_mesh_data(val)
        figure["layout"] = plot_layout
        cs = [[i / (len(colorscale) - 1), rgb] for i, rgb in enumerate(colorscale)]
        figure["data"][0]["colorscale"] = cs
        return figure

    # modify graph markers
    if click_data is not None and "points" in click_data:

        y_value = click_data["points"][0]["y"]
        x_value = click_data["points"][0]["x"]
        z_value = click_data["points"][0]["z"]

        marker = add_marker(x_value, y_value, z_value)
        point_index = marker_in_points(figure["data"], marker)

        # delete graph markers
        if len(figure["data"]) > 1 and point_index is not None:

            figure["data"].pop(point_index)
            anno_index_offset = 2 if val == "mouse" else 1
            try:
                figure["layout"]["scene"]["annotations"].pop(
                    point_index - anno_index_offset
                )
            except Exception as error:
                print(error)
                pass

        # append graph markers
        else:

            # iterate through the store annotations and save it into figure data
            if current_anno is not None:
                for index, annotations in enumerate(
                    figure["layout"]["scene"]["annotations"]
                ):
                    for key in current_anno.keys():
                        if str(index) in key:
                            figure["layout"]["scene"]["annotations"][index][
                                "text"
                            ] = current_anno[key]

            figure["data"].append(marker)
            figure["layout"]["scene"]["annotations"].append(
                add_annotation(x_value, y_value, z_value)
            )

    cs = [[i / (len(colorscale) - 1), rgb] for i, rgb in enumerate(colorscale)]
    figure["data"][0]["colorscale"] = cs

    return figure


@app.callback(Output("click-data", "children"), [Input("brain-graph", "clickData")])
def display_click_data(click_data):
    return json.dumps(click_data, indent=4)


@app.callback(
    Output("relayout-data", "children"), [Input("brain-graph", "relayoutData")]
)
def display_relayout_data(relayout_data):
    return json.dumps(relayout_data, indent=4)


@app.callback(
    Output("annotation_storage", "data"),
    [Input("brain-graph", "relayoutData")],
    [State("annotation_storage", "data")],
)
def save_annotations(relayout_data, current_data):
    """ Update the annotations in the dcc store. """

    if relayout_data is None:
        raise PreventUpdate

    if current_data is None:
        return {}

    for key in relayout_data.keys():

        # to determine if the relayout has to do with annotations
        if "scene.annotations" in key:
            current_data[key] = relayout_data[key]

    return current_data


if __name__ == "__main__":
    app.run_server(debug=True)
