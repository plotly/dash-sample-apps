import json
import dash
from dash import html, dcc, Input, Output, State, callback

import dash_colorscales as dcs
from mni import create_mesh_data

from utils.figures import brain_graph_handler, save_annotations
from utils.components import header, brain_graph, color_picker


app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

app.title = "Brain Surface Viewer"

server = app.server

app.layout = html.Div(
    [
        html.Div(
            [
                html.Div(
                    [
                        header(
                            app,
                            "#141414",
                            "MRI Reconstruction",
                            subheader="Click on the brain to add an annotation. Drag the black corners of the graph to rotate.",
                        ),
                        brain_graph("brain-graph"),
                    ],
                    className="container",
                ),
            ],
            className="two-thirds column app__left__section",
        ),
        color_picker("colorscale-picker"),
        dcc.Store(id="annotation_storage"),
    ]
)


@callback(
    Output("brain-graph", "figure"),
    Input("brain-graph", "clickData"),
    Input("radio-options", "value"),
    Input("colorscale-picker", "colorscale"),
    State("brain-graph", "figure"),
    State("annotation_storage", "data"),
)
def return_brain_graph_handler(click_data, val, colorscale, figure, current_anno):
    return brain_graph_handler(click_data, val, colorscale, figure, current_anno)


@callback(Output("click-data", "children"), Input("brain-graph", "clickData"))
def display_click_data(click_data):
    return json.dumps(click_data, indent=4)


@callback(Output("relayout-data", "children"), Input("brain-graph", "relayoutData"))
def display_relayout_data(relayout_data):
    return json.dumps(relayout_data, indent=4)


@callback(
    Output("annotation_storage", "data"),
    Input("brain-graph", "relayoutData"),
    State("annotation_storage", "data"),
)
def return_save_annotations(relayout_data, current_data):
    return save_annotations(relayout_data, current_data)


if __name__ == "__main__":
    app.run_server(debug=True)
