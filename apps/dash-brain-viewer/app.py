from dash import Dash, dcc, Input, Output, State, callback
import dash_bootstrap_components as dbc
import json

import utils.figures as figs
from utils.components import header, brain_graph, control_and_output


app = Dash(__name__, title = "Brain Surface Viewer", external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

app.layout = dbc.Container(
    [
        dbc.Row([
            dbc.Col([
                header(
                    app,
                    header_color="#F4F6F8",
                    header="MRI Reconstruction",
                    subheader="Click on the brain to add an annotation. Drag the black corners of the graph to rotate.",
                ),
                brain_graph("brain-graph"),
            ], width=8
            ),
            dbc.Col([
                dbc.Card(control_and_output(), className="right-card"),
                
            ], width=4)
        ]),
        
        dcc.Store(id="annotation_storage"),
    ],
    fluid=True
)



@callback(
    Output("brain-graph", "figure"),
    Output("click-data", "children"), 
    Input("brain-graph", "clickData"),
    Input("radio-options", "value"),
    Input("colorscale-picker", "colorscale"),
    State("brain-graph", "figure"),
    State("annotation_storage", "data"),
)
def return_brain_graph_handler_display_data(click_data, val, colorscale, figure, current_anno):
    fig = figs.brain_graph_handler(click_data, val, colorscale, figure, current_anno)
    click_data = json.dumps(click_data, indent=4)
    return fig, click_data

@callback(
    Output("relayout-data", "children"), 
    Output("annotation_storage", "data"),
    Input("brain-graph", "relayoutData"),
    State("annotation_storage", "data"),)
def display_relayout_data_update_storage(relayout_data, current_data):
    return json.dumps(relayout_data, indent=4), figs.save_annotations(relayout_data, current_data)


if __name__ == "__main__":
    app.run_server(debug=True)