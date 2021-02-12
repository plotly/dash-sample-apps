import dash
import dash_vtk
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State

import random
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
            property={
                "edgeVisibility": True,
            },
        )
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
            ],
            style={"height": "12%", "align-items": "center"},
        ),
        html.Div(
            html.Div(vtk_view, style={"height": "100%", "width": "100%"}),
            style={"height": "88%"},
        ),
    ],
)


@app.callback(
    [
        Output("vtk-representation", "colorMapPreset"),
        Output("vtk-representation", "colorDataRange"),
        Output("vtk-polydata", "points"),
        Output("vtk-polydata", "polys"),
        Output("vtk-array", "values"),
        Output("vtk-view", "triggerResetCamera"),
    ],
    [Input("dropdown-preset", "value"), Input("scale-factor", "value")],
)
def updatePresetName(name, scale_factor):
    points, polys, elevation, color_range = updateWarp(scale_factor)
    return [name, color_range, points, polys, elevation, random.random()]


if __name__ == "__main__":
    app.run_server(debug=True)
