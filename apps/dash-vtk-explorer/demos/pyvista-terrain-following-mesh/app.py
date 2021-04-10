import dash
import dash_vtk
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State

import random
import json
import numpy as np
import pyvista as pv
from pyvista import examples
from vtk.util.numpy_support import vtk_to_numpy

from dash_vtk.utils import presets

random.seed(42)


def toDropOption(name):
    return {"label": name, "value": name}


# Get point cloud data from PyVista
uniformGrid = examples.download_crater_topo()
subset = uniformGrid.extract_subset((500, 900, 400, 800, 0, 0), (5, 5, 1))


def updateWarp(factor=1):
    terrain = subset.warp_by_scalar(factor=factor)
    polydata = terrain.extract_geometry()
    points = polydata.points.ravel()
    polys = vtk_to_numpy(polydata.GetPolys().GetData())
    elevation = polydata["scalar1of1"]
    min_elevation = np.amin(elevation)
    max_elevation = np.amax(elevation)
    return [points, polys, elevation, [min_elevation, max_elevation]]


points, polys, elevation, color_range = updateWarp(1)

# Setup VTK rendering of PointCloud
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

vtk_view = dash_vtk.View(
    id="vtk-view",
    pickingModes=["hover"],
    children=[
        dash_vtk.GeometryRepresentation(
            id="vtk-representation",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata",
                    points=points,
                    polys=polys,
                    children=[
                        dash_vtk.PointData(
                            [
                                dash_vtk.DataArray(
                                    id="vtk-array",
                                    registration="setScalars",
                                    name="elevation",
                                    values=elevation,
                                )
                            ]
                        )
                    ],
                )
            ],
            colorMapPreset="erdc_blue2green_muted",
            colorDataRange=color_range,
            property={"edgeVisibility": True},
            showCubeAxes=True,
            cubeAxesStyle={"axisLabels": ["", "", "Altitude"]},
        ),
        dash_vtk.GeometryRepresentation(
            id="pick-rep",
            actor={"visibility": False},
            children=[
                dash_vtk.Algorithm(
                    id="pick-sphere", vtkClass="vtkSphereSource", state={"radius": 100},
                )
            ],
        ),
    ],
)

app.layout = dbc.Container(
    fluid=True,
    style={"height": "100vh"},
    children=[
        dbc.Row(
            [
                dbc.Col(
                    children=dcc.Slider(
                        id="scale-factor",
                        min=0.1,
                        max=5,
                        step=0.1,
                        value=1,
                        marks={0.1: "0.1", 5: "5"},
                    )
                ),
                dbc.Col(
                    children=dcc.Dropdown(
                        id="dropdown-preset",
                        options=list(map(toDropOption, presets)),
                        value="erdc_rainbow_bright",
                    ),
                ),
                dbc.Col(
                    children=dcc.Checklist(
                        id="toggle-cube-axes",
                        options=[{"label": " Show axis grid", "value": "grid"},],
                        value=[],
                        labelStyle={"display": "inline-block"},
                    ),
                ),
            ],
            style={"height": "12%", "align-items": "center"},
        ),
        html.Div(
            html.Div(vtk_view, style={"height": "100%", "width": "100%"}),
            style={"height": "88%"},
        ),
        html.Pre(
            id="tooltip",
            style={
                "position": "absolute",
                "bottom": "25px",
                "left": "25px",
                "zIndex": 1,
                "color": "white",
            },
        ),
    ],
)


@app.callback(
    [
        Output("vtk-representation", "showCubeAxes"),
        Output("vtk-representation", "colorMapPreset"),
        Output("vtk-representation", "colorDataRange"),
        Output("vtk-polydata", "points"),
        Output("vtk-polydata", "polys"),
        Output("vtk-array", "values"),
        Output("vtk-view", "triggerResetCamera"),
    ],
    [
        Input("dropdown-preset", "value"),
        Input("scale-factor", "value"),
        Input("toggle-cube-axes", "value"),
    ],
)
def updatePresetName(name, scale_factor, cubeAxes):
    points, polys, elevation, color_range = updateWarp(scale_factor)
    return [
        "grid" in cubeAxes,
        name,
        color_range,
        points,
        polys,
        elevation,
        random.random(),
    ]


@app.callback(
    [
        Output("tooltip", "children"),
        Output("pick-sphere", "state"),
        Output("pick-rep", "actor"),
    ],
    [Input("vtk-view", "clickInfo"), Input("vtk-view", "hoverInfo"),],
)
def onInfo(clickData, hoverData):
    info = hoverData if hoverData else clickData
    if info:
        if (
            "representationId" in info
            and info["representationId"] == "vtk-representation"
        ):
            return (
                [json.dumps(info, indent=2)],
                {"center": info["worldPosition"]},
                {"visibility": True},
            )
        return dash.no_update, dash.no_update, dash.no_update
    return [""], {}, {"visibility": False}


if __name__ == "__main__":
    app.run_server(debug=True)
