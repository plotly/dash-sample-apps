import random

import dash
import dash_vtk
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State

random.seed(42)

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

vtk_view = dash_vtk.View(
    id="geometry-view",
    children=[
        dash_vtk.GeometryRepresentation(
            [
                dash_vtk.Algorithm(
                    id="vtk-algorithm",
                    vtkClass="vtkConeSource",
                    state={"capping": False, "resolution": 60},
                )
            ]
        )
    ],
)

controls = dbc.Card(
    [
        dbc.CardHeader("Controls"),
        dbc.CardBody(
            [
                html.P("Resolution:"),
                dcc.Slider(
                    id="slider-resolution",
                    min=10,
                    max=100,
                    step=1,
                    value=60,
                    marks={10: "10", 100: "100"},
                ),
                html.Br(),
                dbc.Checklist(
                    options=[{"label": "Capping", "value": "capping"}],
                    value=[],
                    id="capping-checklist",
                    switch=True,
                ),
            ]
        ),
    ]
)

app.layout = dbc.Container(
    fluid=True,
    style={"margin-top": "15px", "height": "calc(100vh - 30px)"},
    children=[
        dbc.Row(
            [
                dbc.Col(width=4, children=controls),
                dbc.Col(
                    width=8,
                    children=[
                        html.Div(vtk_view, style={"height": "100%", "width": "100%"},)
                    ],
                ),
            ],
            style={"height": "100%"},
        ),
    ],
)


@app.callback(
    [Output("vtk-algorithm", "state"), Output("geometry-view", "triggerResetCamera")],
    [Input("slider-resolution", "value"), Input("capping-checklist", "value")],
)
def update_cone(slider_val, checked_values):
    new_state = {"resolution": slider_val, "capping": "capping" in checked_values}
    return new_state, random.random()


if __name__ == "__main__":
    app.run_server(debug=True)
