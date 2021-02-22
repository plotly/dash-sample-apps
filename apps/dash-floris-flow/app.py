# Simpler gist here: https://gist.github.com/xhlulu/773fe238773ea69c8bc2b26560ab67d7
import random

import dash
import dash_bootstrap_components as dbc
import dash_vtk
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.express as px
import floris.tools as wfct
import numpy as np


def Header(name, app):
    title = html.H1(name)
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 60}
    )
    link = html.A(logo, href="https://plotly.com/dash/")

    return dbc.Row([dbc.Col(title, md=8), dbc.Col(link, md=4)], align="center")


def CustomSlider(data: np.array, label: str, id: str):
    n_unique = np.unique(data).shape[0]

    return dbc.FormGroup(
        [
            dbc.Label(label),
            dcc.Slider(
                min=0,
                max=n_unique,
                value=n_unique / 2,
                step=1,
                id=id,
            ),
        ]
    )


def build_view_child(
    dimensions, spacing, origin, field, enabled, i, j, k, window, level
):
    slice_prop = {"colorWindow": window, "colorLevel": level}
    child = [
        dash_vtk.ShareDataSet(
            dash_vtk.ImageData(
                dimensions=dimensions,
                spacing=spacing,
                origin=origin,
                children=dash_vtk.PointData(
                    dash_vtk.DataArray(registration="setScalars", values=field)
                ),
            ),
        ),
    ]

    if "Volume" in enabled:
        child.append(
            dash_vtk.VolumeRepresentation(
                [dash_vtk.VolumeController(), dash_vtk.ShareDataSet()],
            )
        )
    if "i" in enabled:
        child.append(
            dash_vtk.SliceRepresentation(
                iSlice=int(round(i)),
                property=slice_prop,
                children=dash_vtk.ShareDataSet(),
            )
        )

    if "j" in enabled:
        child.append(
            dash_vtk.SliceRepresentation(
                jSlice=int(round(j)),
                property=slice_prop,
                children=dash_vtk.ShareDataSet(),
            )
        )

    if "k" in enabled:
        child.append(
            dash_vtk.SliceRepresentation(
                kSlice=int(round(k)),
                property=slice_prop,
                children=dash_vtk.ShareDataSet(),
            )
        )

    return child


# Initialize the FLORIS interface fi
fi = wfct.floris_interface.FlorisInterface("./data/example_input.json")
fd = fi.get_flow_data()

# compute dimensions, origins and spacing based on flow data
origin = [axis.mean().round().astype(int) for axis in [fd.x, fd.y, fd.z]]
ranges = np.array([axis.ptp().round().astype(int) for axis in [fd.x, fd.y, fd.z]])
dimensions = np.array([np.unique(axis).shape[0] for axis in [fd.x, fd.y, fd.z]])
x, y, z = dimensions
spacing = np.round(ranges / dimensions).astype(int)

# Create the volumes
views = dict(u=fd.u, v=fd.v, w=fd.w)


controls = [
    dbc.FormGroup(
        [
            dbc.Label("Dimension"),
            dbc.RadioItems(
                options=[{"label": x, "value": x} for x in ["u", "v", "w"]],
                value="u",
                id="radio-dimension",
                inline=True,
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Enabled"),
            dbc.Checklist(
                options=[
                    {"label": x, "value": x.replace("Slice ", "")}
                    for x in ["Volume", "Slice i", "Slice j", "Slice k"]
                ],
                value=["Volume", "i"],
                id="child-enabled",
                inline=True,
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Color Window"),
            dcc.Slider(
                min=0.01,
                max=1,
                value=0.5,
                step=0.01,
                id="color-window",
                tooltip={},
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Color Level"),
            dcc.Slider(
                min=0.01,
                max=1,
                value=0.5,
                step=0.01,
                id="color-level",
                tooltip={},
            ),
        ]
    ),
    CustomSlider(data=fd.x, id="slider-slice-i", label="Slice i"),
    CustomSlider(data=fd.y, id="slider-slice-j", label="Slice j"),
    CustomSlider(data=fd.z, id="slider-slice-k", label="Slice k"),
]

vtk_view = dash_vtk.View(id="vtk-view")


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server  # gunicorn needs this for deployment

app.layout = dbc.Container(
    fluid=True,
    style={"height": "100vh"},
    children=[
        Header("Flow Visualization with FLORIS and VTK", app),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    width=4,
                    children=dbc.Card(
                        [
                            dbc.CardHeader("Controls"),
                            dbc.CardBody(controls),
                        ]
                    ),
                ),
                dbc.Col(
                    width=8,
                    children=dbc.Card(
                        [
                            dbc.CardHeader("Flow Visualization"),
                            dbc.CardBody(vtk_view, style={"height": "100%"}),
                        ],
                        style={"height": "80vh"},
                    ),
                ),
            ],
        ),
    ],
)


@app.callback(
    Output("vtk-view", "children"),
    Output("vtk-view", "triggerRender"),
    Input("radio-dimension", "value"),
    Input("child-enabled", "value"),
    Input("color-window", "value"),
    Input("color-level", "value"),
    Input("slider-slice-i", "value"),
    Input("slider-slice-j", "value"),
    Input("slider-slice-k", "value"),
)
def update_flow_viz(
    selected_dim, enabled, window_coef, level_coef, slice_i, slice_j, slice_k
):
    field = views[selected_dim]
    window = (field.max() - field.min()) * window_coef
    level = (field.max() + field.min()) * level_coef

    new_view_child = build_view_child(
        dimensions,
        spacing,
        origin,
        window=window,
        level=level,
        field=field,
        enabled=enabled,
        i=slice_i,
        j=slice_j,
        k=slice_k,
    )
    return new_view_child, random.random()


if __name__ == "__main__":
    app.run_server(debug=True)