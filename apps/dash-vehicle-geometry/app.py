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


def _load_vtp(filepath, fieldname=None):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filepath)
    reader.Update()
    if fieldname == None:
        return to_mesh_state(reader.GetOutput())
    else:
        return to_mesh_state(reader.GetOutput(), fieldname)


# -----------------------------------------------------------------------------
# GUI setup
# -----------------------------------------------------------------------------

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

# -----------------------------------------------------------------------------
# Populate scene
# -----------------------------------------------------------------------------

# vehicle geometry
vehicle_vtk = []
for filename in glob.glob(os.path.join(DATA_PATH, "vehicle") + "/*.vtp"):
    mesh = _load_vtp(filename)
    part_name = filename.split("/")[-1].replace(".vtp", "")
    child = dash_vtk.GeometryRepresentation(
        id=f"{part_name}-rep",
        colorMapPreset="erdc_rainbow_bright",
        colorDataRange=[0, 100],
        property={"opacity": 1},
        children=[dash_vtk.Mesh(id=f"{part_name}-mesh", state=mesh,)],
    )
    vehicle_vtk.append(child)

# isosurfaces
isosurfs_vtk = []
for filename in glob.glob(os.path.join(DATA_PATH, "isosurfaces") + "/*.vtp"):
    mesh = _load_vtp(filename)

    surf_name = filename.split("/")[-1].replace(".vtp", "")
    child = dash_vtk.GeometryRepresentation(
        id=f"{surf_name}-rep",
        property={"opacity": 0, "color": [1, 0, 0]},
        children=[dash_vtk.Mesh(id=f"{surf_name}-mesh", state=mesh,)],
    )

    isosurfs_vtk.append(child)

# -----------------------------------------------------------------------------
# 3D Viz
# -----------------------------------------------------------------------------

vtk_view = dash_vtk.View(id="vtk-view", children=vehicle_vtk + isosurfs_vtk)

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
                    ),
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
                        ],
                        value="solid",
                    ),
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
    children=[
        html.H2("Vehicle Geometry with OpenFOAM"),
        html.Hr(),
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
            style={"margin-top": "15px", "height": "calc(100vh - 230px)"},
        ),
    ],
)

# -----------------------------------------------------------------------------
# Handle controls
# -----------------------------------------------------------------------------


@app.callback(
    [Output("vtk-view", "triggerRender")]
    + [Output(item.id.replace("rep", "mesh"), "state") for item in vehicle_vtk]
    + [Output(item.id, "property") for item in vehicle_vtk]
    + [Output("cp-rep", "property")],
    [
        Input("geometry", "value"),
        Input("isosurfaces", "value"),
        Input("surfcolor", "value"),
    ],
)
def update_scene(geometry, isosurfaces, surfcolor):
    triggered = dash.callback_context.triggered

    # update geometry visibility
    geo_viz, iso_viz = [], []
    if triggered and "geometry" in triggered[0]["prop_id"]:
        geo_viz = [
            {"opacity": 1}
            if part.id.replace("-rep", "").split("_")[0] in triggered[0]["value"]
            else {"opacity": 0}
            for part in vehicle_vtk
        ]
    else:
        geo_viz = [dash.no_update for item in vehicle_vtk]

    # update isosurface visibility
    if triggered and "isosurfaces" in triggered[0]["prop_id"]:
        if "cp" in triggered[0]["value"]:
            iso_viz = {"opacity": 1}
        else:
            iso_viz = {"opacity": 0}
    else:
        iso_viz = dash.no_update

    # update surface coloring
    if triggered and "surfcolor" in triggered[0]["prop_id"]:
        surf_state = []

        for filename in glob.glob(os.path.join(DATA_PATH, "vehicle") + "/*.vtp"):
            if triggered[0]["value"] == "solid":
                mesh = _load_vtp(filename)
            else:
                mesh = _load_vtp(filename, triggered[0]["value"])
            surf_state.append(mesh)
    else:
        surf_state = [dash.no_update for item in vehicle_vtk]

    return [random.random()] + surf_state + geo_viz + [iso_viz]


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    app.run_server(debug=True)
