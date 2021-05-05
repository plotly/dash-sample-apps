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


def cache_mesh(filepath, di, fieldname=None, point_arrays=[], cell_arrays=[]):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filepath)
    reader.Update()

    # Cache mesh for future lookup
    part_name = filepath.split("/")[-1].replace(".vtp", "")
    di[part_name] = reader.GetOutput()


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

for filename in glob.glob(os.path.join(DATA_PATH, "vehicle-vect") + "/*.vtp"):
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

# time to cache meshes
vtk_datasets = {}
for filename in glob.glob(os.path.join(DATA_PATH, "vehicle-vect") + "/*.vtp"):
    cache_mesh(filename, vtk_datasets)

for filename in glob.glob(os.path.join(DATA_PATH, "isosurfaces") + "/*.vtp"):
    cache_mesh(filename, vtk_datasets)

# -----------------------------------------------------------------------------
# 3D Viz
# -----------------------------------------------------------------------------

# vtk_view = dash_vtk.View(id="vtk-view", children=vehicle_vtk + isosurfs_vtk)
# vtk_view = dash_vtk.View(id="vtk-view", pickingModes=['hover'])

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
                            children=[
                                html.Div(
                                    dbc.Spinner(color="light"),
                                    style={
                                        "background-color": "#334c66",
                                        "height": "calc(100vh - 230px)",
                                        "width": "100%",
                                        "text-align": "center",
                                        "padding-top": "calc(50vh - 105px)",
                                    },
                                ),
                            ],
                            id="vtk-view-container",
                            style={"height": "calc(100vh - 230px)", "width": "100%",},
                        ),
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

    cone_pointer = dash_vtk.GeometryRepresentation(
        property={"color": [1, 0, 0]},
        children=[dash_vtk.Algorithm(id="pointer", vtkClass="vtkConeSource")],
    )

    tooltip = html.Pre(
        id="tooltip",
        style={
            "position": "absolute",
            "bottom": "25px",
            "left": "25px",
            "zIndex": 1,
            "color": "white",
        },
    )

    return dash_vtk.View(
        id="vtk-view",
        children=vehicle_vtk + isosurfs_vtk + [cone_pointer, tooltip],
        pickingModes=["hover"],
    )


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


SCALE_P = 0.0001
SCALE_U = 0.01


@app.callback(
    [Output("tooltip", "children"), Output("pointer", "state"),],
    [Input("vtk-view", "hoverInfo"),],
)
def probe_data(info):
    cone_state = {"resolution": 12}
    if info:
        if "representationId" not in info:
            return dash.no_update, dash.no_update
        ds_name = info["representationId"].replace("-rep", "")
        mesh = vtk_datasets[ds_name]
        if mesh:
            xyx = info["worldPosition"]
            idx = mesh.FindPoint(xyx)
            if idx > -1:
                cone_state["center"] = mesh.GetPoints().GetPoint(idx)
                messages = []
                pd = mesh.GetPointData()
                size = pd.GetNumberOfArrays()
                for i in range(size):
                    array = pd.GetArray(i)
                    name = array.GetName()
                    nb_comp = array.GetNumberOfComponents()
                    value = array.GetValue(idx)
                    value_str = f"{array.GetValue(idx):.2f}"
                    norm_str = ""
                    if nb_comp == 3:
                        value = array.GetTuple3(idx)
                        norm = (value[0] ** 2 + value[1] ** 2 + value[2] ** 2) ** 0.5
                        norm_str = f" norm({norm:.2f})"
                        value_str = ", ".join([f"{v:.2f}" for v in value])

                        cone_state["height"] = SCALE_U * norm
                        cone_state["direction"] = [v / norm for v in value]

                    if name == "p":
                        cone_state["radius"] = array.GetValue(idx) * SCALE_P

                    messages.append(f"{name}: {value_str} {norm_str}")

        if "height" in cone_state:
            new_center = [v for v in cone_state["center"]]
            for i in range(3):
                new_center[i] -= 0.5 * cone_state["height"] * cone_state["direction"][i]
            cone_state["center"] = new_center

        return ["\n".join(messages)], cone_state
    return [""], cone_state


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    app.run_server(debug=True)
