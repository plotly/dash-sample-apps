import os
import random


import dash
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc

from dash.dependencies import Input, Output, State

import dash_vtk
from dash_vtk.utils import to_mesh_state, preset_as_options

import vtk

random.seed(42)

# -----------------------------------------------------------------------------
# VTK Pipeline
# -----------------------------------------------------------------------------


class Viz:
    def __init__(self, data_directory):
        self.color_range = [0, 1]
        bike_filename = os.path.join(data_directory, "bike.vtp")
        tunnel_filename = os.path.join(data_directory, "tunnel.vtu")

        # Seeds settings
        self.resolution = 10
        self.point1 = [-0.4, 0, 0.05]
        self.point2 = [-0.4, 0, 1.5]

        # VTK Pipeline setup
        bikeReader = vtk.vtkXMLPolyDataReader()
        bikeReader.SetFileName(bike_filename)
        bikeReader.Update()
        self.bike_mesh = to_mesh_state(bikeReader.GetOutput())

        tunnelReader = vtk.vtkXMLUnstructuredGridReader()
        tunnelReader.SetFileName(tunnel_filename)
        tunnelReader.Update()

        self.lineSeed = vtk.vtkLineSource()
        self.lineSeed.SetPoint1(*self.point1)
        self.lineSeed.SetPoint2(*self.point2)
        self.lineSeed.SetResolution(self.resolution)

        streamTracer = vtk.vtkStreamTracer()
        streamTracer.SetInputConnection(tunnelReader.GetOutputPort())
        streamTracer.SetSourceConnection(self.lineSeed.GetOutputPort())
        streamTracer.SetIntegrationDirectionToForward()
        streamTracer.SetIntegratorTypeToRungeKutta45()
        streamTracer.SetMaximumPropagation(3)
        streamTracer.SetIntegrationStepUnit(2)
        streamTracer.SetInitialIntegrationStep(0.2)
        streamTracer.SetMinimumIntegrationStep(0.01)
        streamTracer.SetMaximumIntegrationStep(0.5)
        streamTracer.SetMaximumError(0.000001)
        streamTracer.SetMaximumNumberOfSteps(2000)
        streamTracer.SetTerminalSpeed(0.00000000001)

        self.tubeFilter = vtk.vtkTubeFilter()
        self.tubeFilter.SetInputConnection(streamTracer.GetOutputPort())
        self.tubeFilter.SetRadius(0.01)
        self.tubeFilter.SetNumberOfSides(6)
        self.tubeFilter.CappingOn()
        self.tubeFilter.Update()

    def updateSeedPoints(self, p1_y, p2_y, resolution):
        self.point1[1] = p1_y
        self.point2[1] = p2_y
        self.resolution = resolution

        self.lineSeed.SetPoint1(*self.point1)
        self.lineSeed.SetPoint2(*self.point2)
        self.lineSeed.SetResolution(resolution)

    def getTubesMesh(self, color_by_field_name):
        self.tubeFilter.Update()
        ds = self.tubeFilter.GetOutput()
        mesh_state = to_mesh_state(ds, color_by_field_name)
        self.color_range = mesh_state["field"]["dataRange"]
        return mesh_state

    def getBikeMesh(self):
        return self.bike_mesh

    def getColorRange(self):
        return self.color_range

    def getSeedState(self):
        return {
            "point1": self.point1,
            "point2": self.point2,
            "resolution": self.resolution,
        }


# -----------------------------------------------------------------------------
# GUI setup
# -----------------------------------------------------------------------------

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server
viz = Viz(os.path.join(os.path.dirname(__file__), "data"))

# -----------------------------------------------------------------------------
# 3D Viz
# -----------------------------------------------------------------------------

vtk_view = dash_vtk.View(
    id="vtk-view",
    children=[
        dash_vtk.GeometryRepresentation(
            id="bike-rep",
            children=[dash_vtk.Mesh(id="bike", state=viz.getBikeMesh(),)],
        ),
        dash_vtk.GeometryRepresentation(
            id="tubes-rep",
            colorMapPreset="erdc_rainbow_bright",
            colorDataRange=viz.getColorRange(),
            children=[dash_vtk.Mesh(id="tubes-mesh", state=viz.getTubesMesh("p"),)],
        ),
        dash_vtk.GeometryRepresentation(
            id="seed-rep",
            property={"color": [0.8, 0, 0], "representation": 0, "pointSize": 8,},
            children=[
                dash_vtk.Algorithm(
                    id="seed-line", vtkClass="vtkLineSource", state=viz.getSeedState(),
                )
            ],
        ),
    ],
)

# -----------------------------------------------------------------------------
# Control UI
# -----------------------------------------------------------------------------

controls = [
    dbc.Card(
        [
            dbc.CardHeader("Seeds"),
            dbc.CardBody(
                [
                    html.P("Seed line:"),
                    dcc.Slider(
                        id="point-1",
                        min=-1,
                        max=1,
                        step=0.01,
                        value=0,
                        marks={-1: "-1", 1: "+1"},
                    ),
                    dcc.Slider(
                        id="point-2",
                        min=-1,
                        max=1,
                        step=0.01,
                        value=0,
                        marks={-1: "-1", 1: "+1"},
                    ),
                    html.Br(),
                    html.P("Line resolution:"),
                    dcc.Slider(
                        id="seed-resolution",
                        min=5,
                        max=50,
                        step=1,
                        value=10,
                        marks={5: "5", 50: "50"},
                    ),
                ]
            ),
        ]
    ),
    html.Br(),
    dbc.Card(
        [
            dbc.CardHeader("Color By"),
            dbc.CardBody(
                [
                    html.P("Field name"),
                    dcc.Dropdown(
                        id="color-by",
                        options=[
                            {"label": "p", "value": "p"},
                            {"label": "Rotation", "value": "Rotation"},
                            {"label": "U", "value": "U"},
                            {"label": "Vorticity", "value": "Vorticity"},
                            {"label": "k", "value": "k"},
                        ],
                        value="p",
                    ),
                    html.Br(),
                    html.P("Color Preset"),
                    dcc.Dropdown(
                        id="preset",
                        options=preset_as_options,
                        value="erdc_rainbow_bright",
                    ),
                ]
            ),
        ]
    ),
]

# -----------------------------------------------------------------------------
# App UI
# -----------------------------------------------------------------------------

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
                        html.Div(vtk_view, style={"height": "100%", "width": "100%"})
                    ],
                ),
            ],
            style={"height": "100%"},
        ),
    ],
)

# -----------------------------------------------------------------------------
# Handle controls
# -----------------------------------------------------------------------------


@app.callback(
    [
        Output("seed-line", "state"),
        Output("tubes-mesh", "state"),
        Output("tubes-rep", "colorDataRange"),
        Output("tubes-rep", "colorMapPreset"),
        Output("vtk-view", "triggerRender"),
    ],
    [
        Input("point-1", "value"),
        Input("point-2", "value"),
        Input("seed-resolution", "value"),
        Input("color-by", "value"),
        Input("preset", "value"),
    ],
)
def update_seeds(y1, y2, resolution, colorByField, presetName):
    viz.updateSeedPoints(y1, y2, resolution)
    return [
        viz.getSeedState(),
        viz.getTubesMesh(colorByField),
        viz.getColorRange(),
        presetName,
        random.random(),  # trigger a render
    ]


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    app.run_server(debug=True)
