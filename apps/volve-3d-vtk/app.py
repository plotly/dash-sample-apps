#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import pyvista as pv
import xarray as xr
from jupyter_dash import JupyterDash
import dash
import dash_html_components as html
import vtk
import dash_vtk
from dash_vtk.utils import to_volume_state
from pyvista import examples
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import os
from dash_vtk.utils import presets
import segyio
import random
import dash_core_components as dcc
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import PVGeo


# # Data Import
# ## Wells

# In[9]:


# import of the wells and store them in a dict
surveys = {}
for filename in os.listdir(r"data/Wells_seismic"):
    surveys[filename[:-4]] = pd.read_csv(r"data/Wells_seismic/{}".format(filename))


# In[10]:


# generation of numpy arrays of xyz/inline-xline-z for each well
points_dict = {}
for points in surveys:
    points_dict["{}".format(points)] = np.array(
        list(zip(surveys[points].inline, surveys[points].xline, surveys[points].TVD))
    )


# In[11]:


# generation of lines from the numpy arrays by using PVGeo for each well
lines_dict = {}
for lines in points_dict:
    poly = PVGeo.points_to_poly_data(points_dict[lines])
    lines_dict["{}".format(lines)] = PVGeo.filters.AddCellConnToPoints().apply(poly)


# In[12]:


# from the previous lines, a dictionary of vtk points and vtk lines was generated for each well
points_well = {}
lines_well = {}
for lines in lines_dict:
    points_well[lines] = lines_dict[lines].points.ravel()
    lines_well[lines] = lines_dict[lines].lines.ravel()


# ## Horizon

# In[13]:


# import of the top of the reservoir horizon
reservoir_hor = pd.read_csv(
    r"data/Horizons/Hugin-reservoir.dat",
    names=["x", "y", "inline", "xline", "z"],
    sep="\t",
)

# to generate the vtk surface of the horizon, we are going to divide the generation of the polys and the points,
# the polys are being generated from the x, y, z data and the points from inline, xline and z, this in order to get correct the scale of the horizon.


# ### Polys

# In[14]:


# convertion of the x,y,z to numpy for polys
hor_array_polys = reservoir_hor.loc[:, ["x", "y", "z"]].to_numpy()

# use of pyvista to generate a polydata object
hor_cloud_polys = pv.PolyData(hor_array_polys)

# assign the depth coordinate as the Depth value
hor_cloud_polys["Depth"] = hor_array_polys[:, -1]

# generation of a surface from the polydata by using ptvista delaunay
surf_polys = hor_cloud_polys.delaunay_2d()


# In[15]:


# from the polydata object, we get the vtk polys for the 3d visualization
polydata_polys = surf_polys.extract_geometry()
polys_hor = vtk_to_numpy(polydata_polys.GetPolys().GetData())

# extract the depth values of the horizon and range of depth for the 3d visualization
depth = polydata_polys["Depth"]
min_depth = np.amin(depth)
max_depth = np.amax(depth)
color_range = [min_depth, max_depth]


# ### Points

# In[16]:


# convertion of the inline, xline,z to numpy for points
hor_array_points = reservoir_hor.loc[:, ["inline", "xline", "z"]].to_numpy()

# use of pyvista to generate a polydata object
hor_cloud_points = pv.PolyData(hor_array_points)

# generation of a surface from the polydata by using ptvista delaunay
surf_points = hor_cloud_points.delaunay_2d()


# In[17]:


# from the polydata object, we get the vtk points for the 3d visualization
polydata_points = surf_points.extract_geometry()
points_hor = polydata_points.points.ravel()


# ## Seismic

# In[18]:


# import of the seismic in numpy format
# seis = np.load('data/Seismic/final_cube.npy')
# seis.shape


# In[20]:


# z = seis[:,:,0:375:25]
# np.save(file='data/Seismic/final_cube_z', arr=z[:,:,:])
z = np.load("data/Seismic/final_cube_z.npy")

z_mesh = pv.wrap(z)

z_slice = (
    dash_vtk.ImageData(
        dimensions=[231, 392, 15],
        origin=[10030, 2145, 0],
        spacing=[1, 1, 250],
        children=[
            dash_vtk.PointData(
                [
                    dash_vtk.DataArray(
                        registration="setScalars", values=z_mesh["values"] * 10000,
                    )
                ]
            )
        ],
    ),
)


# In[21]:


# x = seis[0:220:20,:,:]
# np.save(file='data/Seismic/final_cube_in', arr=z[:,:,:])

x = np.load("data/Seismic/final_cube_in.npy")

x_mesh = pv.wrap(x)

x_slice = (
    dash_vtk.ImageData(
        dimensions=[11, 392, 376],
        origin=[10030, 2145, 0],
        spacing=[20, 1, 10],
        children=[
            dash_vtk.PointData(
                [
                    dash_vtk.DataArray(
                        registration="setScalars", values=x_mesh["values"] * 10000,
                    )
                ]
            )
        ],
    ),
)


# In[22]:


# y= seis[:,0:380:20,:]
# np.save(file='data/Seismic/final_cube_xl', arr=z[:,:,:])

y = np.load("data/Seismic/final_cube_xl.npy")

y_mesh = pv.wrap(y)

y_slice = (
    dash_vtk.ImageData(
        dimensions=[231, 19, 376],
        origin=[10030, 2145, 0],
        spacing=[1, 20, 10],
        children=[
            dash_vtk.PointData(
                [
                    dash_vtk.DataArray(
                        registration="setScalars", values=y_mesh["values"] * 10000,
                    )
                ]
            )
        ],
    ),
)


# ## Grid

# In[23]:


grid_points = np.array(
    [
        100300,
        21450,
        0,
        100300,
        25360,
        0,
        102600,
        21450,
        0,
        102600,
        25360,
        0,
        100300,
        21450,
        3750,
        100300,
        25360,
        3750,
        102600,
        21450,
        3750,
        102600,
        25360,
        3750,
    ]
)


# # VTK Visualization

# ## View

# In[312]:


slice_view = dash_vtk.View(
    id="slice-view",
    cameraPosition=[1, 0, 0],
    cameraViewUp=[0, 0, -1],
    cameraParallelProjection=False,
    background=[0.1137, 0.2078, 0.3411],
    children=[
        dash_vtk.SliceRepresentation(
            id="slice-repr-in",
            iSlice=5,
            property={"colorWindow": 2000, "colorLevel": 0},
            actor={"scale": (10, 10, 1)},
            colorMapPreset="Grayscale",
            children=x_slice,
        ),
        dash_vtk.SliceRepresentation(
            id="slice-repr-xl",
            jSlice=8,
            actor={"scale": (10, 10, 1)},
            property={"colorWindow": 2000, "colorLevel": 0},
            colorMapPreset="Grayscale",
            children=y_slice,
        ),
        dash_vtk.SliceRepresentation(
            id="slice-repr-z",
            kSlice=2,
            actor={"scale": (10, 10, 1)},
            property={"colorWindow": 2000, "colorLevel": 0},
            colorMapPreset="Grayscale",
            children=z_slice,
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-15-9_19A",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata1",
                    points=points_well["15-9_19A"],
                    lines=lines_well["15-9_19A"],
                ),
            ],
            property={"edgeVisibility": False, "lineWidth": 100, "color": (255, 0, 0)},
            actor={"scale": (10, 10, 1)},
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-15-9_19BT2",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata2",
                    points=points_well["15-9_19BT2"],
                    lines=lines_well["15-9_19BT2"],
                )
            ],
            property={"edgeVisibility": False, "lineWidth": 100, "color": (255, 0, 0)},
            actor={"scale": (10, 10, 1)},
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-15-9_19SR",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata3",
                    points=points_well["15-9_19SR"],
                    lines=lines_well["15-9_19SR"],
                )
            ],
            property={"edgeVisibility": False, "lineWidth": 100, "color": (255, 0, 0)},
            actor={"scale": (10, 10, 1)},
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-15-9_F11B",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata4",
                    points=points_well["15-9_F11B"],
                    lines=lines_well["15-9_F11B"],
                )
            ],
            property={"edgeVisibility": False, "lineWidth": 100, "color": (255, 0, 0)},
            actor={"scale": (10, 10, 1)},
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-15-9_F12",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata5",
                    points=points_well["15-9_F12"],
                    lines=lines_well["15-9_F12"],
                )
            ],
            property={"edgeVisibility": False, "lineWidth": 100, "color": (255, 0, 0)},
            actor={"scale": (10, 10, 1)},
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-15-9_F14",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata6",
                    points=points_well["15-9_F14"],
                    lines=lines_well["15-9_F14"],
                )
            ],
            property={"edgeVisibility": False, "lineWidth": 100, "color": (255, 0, 0)},
            actor={"scale": (10, 10, 1)},
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-15-9_F15C",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata7",
                    points=points_well["15-9_F15C"],
                    lines=lines_well["15-9_F15C"],
                )
            ],
            property={"edgeVisibility": False, "lineWidth": 100, "color": (255, 0, 0)},
            actor={"scale": (10, 10, 1)},
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-15-9_F1C",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata8",
                    points=points_well["15-9_F1C"],
                    lines=lines_well["15-9_F1C"],
                )
            ],
            property={"edgeVisibility": False, "lineWidth": 100, "color": (255, 0, 0)},
            actor={"scale": (10, 10, 1)},
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-horizon",
            children=[
                dash_vtk.PolyData(
                    id="vtk-polydata9",
                    points=points_hor,
                    polys=polys_hor,
                    children=[
                        dash_vtk.PointData(
                            [
                                dash_vtk.DataArray(
                                    id="vtk-array",
                                    registration="setScalars",
                                    name="Depth",
                                    values=depth,
                                )
                            ]
                        )
                    ],
                )
            ],
            colorMapPreset="erdc_rainbow_bright",
            colorDataRange=color_range,
            property={"edgeVisibility": False},
            actor={"scale": (10, 10, 1)},
        ),
        dash_vtk.GeometryRepresentation(
            id="vtk-grid",
            children=[dash_vtk.PolyData(id="vtk-polydata10", points=grid_points)],
            property={"edgeVisibility": False},
            #             actor={'scale':(10,10,1)},
            showCubeAxes=True,
            cubeAxesStyle={"axisLabels": ["Inline", "Xline", "Depth (m)"]},
        ),
    ],
)


# In[313]:


controls = [
    dbc.FormGroup(
        [
            dbc.Label("Display Wells:"),
            dbc.RadioItems(
                options=[
                    {"label": "Yes", "value": "Yes"},
                    {"label": "No", "value": "No"},
                ],
                value="Yes",
                id="wells_id",
                inline=True,
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Display Reservoir Horizon:"),
            dbc.RadioItems(
                options=[
                    {"label": "Yes", "value": "Yes"},
                    {"label": "No", "value": "No"},
                ],
                value="Yes",
                id="horizon_id",
                inline=True,
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Display seismic:"),
            dbc.Checklist(
                options=[
                    {"label": "Inline", "value": "inline"},
                    {"label": "Xline", "value": "xline"},
                    {"label": "Depth Slice", "value": "z"},
                ],
                value=["inline", "xline", "z"],
                id="enabled",
                inline=True,
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Seismic controls"),
            html.Br(),
            html.Div(id="inline_val", style={"margin-right": 0, "margin-left": 0,}),
            dcc.Slider(id="slider-in", min=0, max=11, value=5),
            html.Div(id="xline_val", style={"margin-right": 0, "margin-left": 0,}),
            dcc.Slider(id="slider-xl", min=0, max=19, value=9),
            html.Div(id="z_val", style={"margin-right": 0, "margin-left": 0,}),
            dcc.Slider(id="slider-z", min=0, max=15, value=2),
        ],
    ),
]


# In[314]:


app = dash.Dash(__name__)

server = app.server

app.layout = dbc.Container(
    fluid=True,
    children=[
        html.Div(
            [
                html.H3(
                    "Geoscience Data Visualizer",
                    style={
                        "fontSize": "2vw",
                        "lineHeight": 1.3,
                        "letterSpacing": "-1px",
                        "marginBottom": "0px",
                        "textAlign": "center",
                        "marginTop": "0px",
                        "fontFamily": "sans-serif",
                    },
                ),
                html.H5(
                    "Volve Dataset",
                    style={
                        "fontSize": "1.5vw",
                        "lineHeight": 1.5,
                        "letterSpacing": "-0.5px",
                        "marginBottom": "0px",
                        "textAlign": "center",
                        "marginTop": "0px",
                        "fontFamily": "sans-serif",
                    },
                ),
                dcc.Markdown(
                    "Volve is a Norwegian Oil Field in the North Sea, that in 2018, Equinor and the Volve license partners decided to releasea a big dataset containing all kinds of data, from seismic to production related to the field, all of these under an Open License, that you can check [here.](https://www.equinor.com/content/dam/statoil/documents/what-we-do/Equinor-HRS-Terms-and-conditions-for-licence-to-data-Volve.pdf)",
                    style={
                        "fontSize": ".7vw",
                        "marginBottom": "1vh",
                        "textAlign": "center",
                        "marginTop": "0px",
                        "fontFamily": "sans-serif",
                    },
                ),
            ]
        ),
        html.Div(
            [
                html.Div(style={"margin-bottom": "9px"}),
                dbc.Row(
                    [
                        dbc.Col(
                            dbc.Card(
                                [dbc.CardBody(controls),], style={"width": "15vw"}
                            ),
                            style={"margin-left": "5vw"},
                            width={"size": 2},
                        ),
                        dbc.Col(
                            dbc.Card(
                                [dbc.CardBody(slice_view),],
                                style={"height": "85vh", "width": "73vw"},
                            ),
                            width={"size": 9},
                        ),
                    ]
                ),
            ],
            style={"margin-bottom": "0px"},
        ),
    ],
)


# In[315]:


@app.callback(
    [
        Output("inline_val", "children"),
        Output("xline_val", "children"),
        Output("z_val", "children"),
    ],
    [
        Input("slider-in", "value"),
        Input("slider-xl", "value"),
        Input("slider-z", "value"),
    ],
)
def inl_value(in_val, xl_val, z_val):
    inl_dict = {}
    inl = 100300
    for i in range(0, 12):
        inl_dict[i] = inl
        inl += 200
    in_display = inl_dict[in_val]

    xl_dict = {}
    xl = 21450
    for i in range(0, 20):
        xl_dict[i] = xl
        xl += 200
    xl_display = xl_dict[xl_val]

    z_dict = {}
    z = 0
    for i in range(0, 16):
        z_dict[i] = z
        z += 250
    z_display = z_dict[z_val]

    return (
        f"Inlines: {in_display}",
        f"Xlines: {xl_display}",
        f"Depth slice: {z_display} m",
    )


@app.callback(
    [
        Output("slice-view", "triggerRender"),
        Output("slice-repr-in", "iSlice"),
        Output("slice-repr-xl", "jSlice"),
        Output("slice-repr-z", "kSlice"),
        Output("vtk-15-9_19A", "actor"),
        Output("vtk-15-9_19BT2", "actor"),
        Output("vtk-15-9_19SR", "actor"),
        Output("vtk-15-9_F11B", "actor"),
        Output("vtk-15-9_F12", "actor"),
        Output("vtk-15-9_F14", "actor"),
        Output("vtk-15-9_F15C", "actor"),
        Output("vtk-15-9_F1C", "actor"),
        Output("vtk-horizon", "actor"),
        Output("slice-repr-in", "actor"),
        Output("slice-repr-xl", "actor"),
        Output("slice-repr-z", "actor"),
    ],
    [
        Input("slider-in", "value"),
        Input("slider-xl", "value"),
        Input("slider-z", "value"),
        Input("wells_id", "value"),
        Input("horizon_id", "value"),
        Input("enabled", "value"),
    ],
)
def update_slice_property(i, j, k, well, hor, seis):
    render_call = random.random()
    if well == "No":
        act = {"scale": (10, 10, 1), "visibility": 0}
        act2 = {"scale": (10, 10, 1), "visibility": 0}
        act3 = {"scale": (10, 10, 1), "visibility": 0}
        act4 = {"scale": (10, 10, 1), "visibility": 0}
        act5 = {"scale": (10, 10, 1), "visibility": 0}
        act6 = {"scale": (10, 10, 1), "visibility": 0}
        act7 = {"scale": (10, 10, 1), "visibility": 0}
        act8 = {"scale": (10, 10, 1), "visibility": 0}
    elif well == "Yes":
        act = {"scale": (10, 10, 1), "visibility": 1}
        act2 = {"scale": (10, 10, 1), "visibility": 1}
        act3 = {"scale": (10, 10, 1), "visibility": 1}
        act4 = {"scale": (10, 10, 1), "visibility": 1}
        act5 = {"scale": (10, 10, 1), "visibility": 1}
        act6 = {"scale": (10, 10, 1), "visibility": 1}
        act7 = {"scale": (10, 10, 1), "visibility": 1}
        act8 = {"scale": (10, 10, 1), "visibility": 1}

    if hor == "No":
        act9 = {"scale": (10, 10, 1), "visibility": 0}
    elif hor == "Yes":
        act9 = {"scale": (10, 10, 1), "visibility": 1}

    if "inline" in seis:
        act10 = {"scale": (10, 10, 1), "visibility": 1}
    elif "inline" not in seis:
        act10 = {"scale": (10, 10, 1), "visibility": 0}

    if "xline" in seis:
        act11 = {"scale": (10, 10, 1), "visibility": 1}
    elif "xline" not in seis:
        act11 = {"scale": (10, 10, 1), "visibility": 0}

    if "z" in seis:
        act12 = {"scale": (10, 10, 1), "visibility": 1}
    elif "z" not in seis:
        act12 = {"scale": (10, 10, 1), "visibility": 0}

    return (
        render_call,
        i,
        j,
        k,
        act,
        act2,
        act3,
        act4,
        act5,
        act6,
        act7,
        act8,
        act9,
        act10,
        act11,
        act12,
    )


if __name__ == "__main__":
    app.run_server(debug=True)


# In[ ]:
