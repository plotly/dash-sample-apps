import os
import glob
import random

import dash
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc

from dash.dependencies import Input, Output, State

import dash_vtk
from dash_vtk.utils import to_mesh_state, preset_as_options

import vtk

DATA_PATH = "data"

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------


def _load_vtp(filepath, fieldname=None, point_arrays=[], cell_arrays=[]):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filepath)
    reader.Update()
    return to_mesh_state(reader.GetOutput(), fieldname, point_arrays, cell_arrays)


# -----------------------------------------------------------------------------
# GUI setup
# -----------------------------------------------------------------------------

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True,
)
server = app.server

# -----------------------------------------------------------------------------
# Populate scene
# -----------------------------------------------------------------------------

# vehicle geometry
vehicle_vtk = []
vehicle_mesh_ids = []
vehicle_meshes = []

for filename in glob.glob(os.path.join(DATA_PATH, "vehicle") + "/*.vtp"):
    mesh = _load_vtp(filename, point_arrays=["U", "p"])
    part_name = filename.split("/")[-1].replace(".vtp", "")
    child = dash_vtk.GeometryRepresentation(
        id=f"{part_name}-rep",
        colorMapPreset="erdc_rainbow_bright",
        colorDataRange=[0, 100],
        actor={"visibility": 1},
        mapper={"scalarVisibility": False},
        children=[dash_vtk.Mesh(id=f"{part_name}-mesh", state=mesh)],
        # children=[dash_vtk.Mesh(id=f"{part_name}-mesh")],
    )
    vehicle_vtk.append(child)

    vehicle_mesh_ids.append(f"{part_name}-mesh")
    vehicle_meshes.append(mesh)

# isosurfaces
isosurfs_vtk = []
isomesh_ids = []
isosurfs_meshes = []

for filename in glob.glob(os.path.join(DATA_PATH, "isosurfaces") + "/*.vtp"):
    mesh = _load_vtp(filename)

    surf_name = filename.split("/")[-1].replace(".vtp", "")
    child = dash_vtk.GeometryRepresentation(
        id=f"{surf_name}-rep",
        property={"color": [1, 0, 0]},
        actor={"visibility": 0},
        children=[dash_vtk.Mesh(id=f"{surf_name}-mesh", state=mesh)],
        # children=[dash_vtk.Mesh(id=f"{surf_name}-mesh")],
    )

    isosurfs_vtk.append(child)

    isomesh_ids.append(f"{surf_name}-mesh")
    isosurfs_meshes.append(mesh)


# -----------------------------------------------------------------------------
# 3D Viz
# -----------------------------------------------------------------------------

# vtk_view = dash_vtk.View(id="vtk-view", children=vehicle_vtk + isosurfs_vtk)
vtk_view = dash_vtk.View(id="vtk-view")

# -----------------------------------------------------------------------------
# Control UI
# -----------------------------------------------------------------------------

controls = [
    dbc.Card(
        [
            dbc.CardHeader("Geometry"),
            dbc.CardBody(
                [
                    dcc.Checklist(
                        id="geometry",
                        options=[
                            {"label": " body", "value": "body"},
                            {"label": " drivetrain", "value": "drive-train"},
                            {"label": " front-wing", "value": "front-wing"},
                            {"label": " rear-wing", "value": "rear-wing"},
                        ],
                        labelStyle={"display": "block"},
                        value=["body", "drive-train", "front-wing", "rear-wing"],
                    )
                ]
            ),
        ]
    ),
    html.Br(),
    dbc.Card(
        [
            dbc.CardHeader("Surface Coloring"),
            dbc.CardBody(
                [
                    dcc.Dropdown(
                        id="surfcolor",
                        options=[
                            {"label": "solid", "value": "solid"},
                            {"label": "U", "value": "U"},
                            {"label": "p", "value": "p"},
                        ],
                        value="solid",
                    )
                ]
            ),
        ]
    ),
    html.Br(),
    dbc.Card(
        [
            dbc.CardHeader("Isosurfaces"),
            dbc.CardBody(
                [
                    dcc.Checklist(
                        id="isosurfaces",
                        options=[{"label": " Cp", "value": "cp"}],
                        labelStyle={"display": "block"},
                        value=[],
                    )
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
    children=[
        html.H2("Vehicle Geometry with OpenFOAM"),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(width=4, children=controls),
                dbc.Col(
                    width=8,
                    children=[
                        html.Div(
                            dbc.Spinner(
                                html.Div(
                                    id="vtk-view-container",
                                    style={
                                        "height": "calc(100vh - 230px)",
                                        "width": "100%",
                                    },
                                ),
                                color="light",
                            ),
                            style={"background-color": "#334c66"},
                        )
                    ],
                ),
            ],
            style={"margin-top": "15px", "height": "calc(100vh - 230px)"},
        ),
    ],
)

# -----------------------------------------------------------------------------
# This Handle controls
# -----------------------------------------------------------------------------

COLOR_RANGES = {"solid": [0, 1], "U": [0, 100], "p": [-4464, 1700]}


@app.callback(Output("vtk-view-container", "children"), [Input("geometry", "value")])
def initial_loading(geometry):
    triggered = dash.callback_context.triggered
    if triggered:
        return dash.no_update

    return dash_vtk.View(id="vtk-view", children=vehicle_vtk + isosurfs_vtk)


@app.callback(
    [Output("vtk-view", "triggerRender")]
    + [Output(item.id, "mapper") for item in vehicle_vtk]
    + [Output(item.id, "actor") for item in vehicle_vtk]
    + [Output(item.id, "colorDataRange") for item in vehicle_vtk]
    + [Output("cp-rep", "actor")],
    [
        Input("geometry", "value"),
        Input("isosurfaces", "value"),
        Input("surfcolor", "value"),
    ],
    prevent_initial_call=True,
)
def update_scene(geometry, isosurfaces, surfcolor):
    triggered = dash.callback_context.triggered

    # update geometry visibility
    geo_viz, iso_viz = [], []
    if triggered and "geometry" in triggered[0]["prop_id"]:
        geo_viz = [
            {"visibility": 1}
            if part.id.replace("-rep", "").split("_")[0] in triggered[0]["value"]
            else {"visibility": 0}
            for part in vehicle_vtk
        ]
    else:
        geo_viz = [dash.no_update for item in vehicle_vtk]

    # update isosurface visibility
    if triggered and "isosurfaces" in triggered[0]["prop_id"]:
        if "cp" in triggered[0]["value"]:
            iso_viz = {"visibility": 1}
        else:
            iso_viz = {"visibility": 0}
    else:
        iso_viz = dash.no_update

    # update surface coloring
    if triggered and "surfcolor" in triggered[0]["prop_id"]:
        color_range = COLOR_RANGES[triggered[0]["value"]]
        mapper = {
            "colorByArrayName": triggered[0]["value"],
            "scalarMode": 3,
            "interpolateScalarsBeforeMapping": True,
            "scalarVisibility": True,
        }
        if triggered[0]["value"] == "solid":
            mapper = {"scalarVisibility": False}

        surf_state = [mapper for item in vehicle_vtk]
        color_ranges = [color_range for item in vehicle_vtk]
    else:
        surf_state = [dash.no_update] * len(vehicle_vtk)
        color_ranges = [dash.no_update] * len(vehicle_vtk)

    return [random.random()] + surf_state + geo_viz + color_ranges + [iso_viz]


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    app.run_server(debug=True)
