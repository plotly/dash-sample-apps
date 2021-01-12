from time import time

import numpy as np
from nilearn import image
from skimage import draw, filters, exposure, measure
from scipy import ndimage

import plotly.graph_objects as go
import plotly.express as px

import dash
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
from dash_slicer import VolumeSlicer


external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, update_title=None, external_stylesheets=external_stylesheets)
server = app.server


t1 = time()

# ------------- I/O and data massaging ---------------------------------------------------

img = image.load_img("assets/radiopaedia_org_covid-19-pneumonia-7_85703_0-dcm.nii")
mat = img.affine
img = img.get_data()
img = np.copy(np.moveaxis(img, -1, 0))[:, ::-1]

spacing = abs(mat[2, 2]), abs(mat[1, 1]), abs(mat[0, 0])

# Create smoothed image and histogram
med_img = filters.median(img, selem=np.ones((1, 3, 3), dtype=np.bool))
hi = exposure.histogram(med_img)

# Create mesh
verts, faces, _, _ = measure.marching_cubes(med_img, 200, step_size=5)
x, y, z = verts.T
i, j, k = faces.T
fig_mesh = go.Figure()
fig_mesh.add_trace(go.Mesh3d(x=z, y=y, z=x, opacity=0.2, i=k, j=j, k=i))

# Create slicers
slicer1 = VolumeSlicer(app, img, axis=0, spacing=spacing, thumbnail=False)
slicer1.graph.figure.update_layout(
    dragmode="drawclosedpath", newshape_line_color="cyan", plot_bgcolor="rgb(0, 0, 0)"
)
slicer1.graph.config.update(
    modeBarButtonsToAdd=["drawclosedpath", "eraseshape",]
)

slicer2 = VolumeSlicer(app, img, axis=1, spacing=spacing, thumbnail=False)
slicer2.graph.figure.update_layout(
    dragmode="drawrect", newshape_line_color="cyan", plot_bgcolor="rgb(0, 0, 0)"
)
slicer2.graph.config.update(
    modeBarButtonsToAdd=["drawrect", "eraseshape",]
)


def path_to_coords(path):
    """From SVG path to numpy array of coordinates, each row being a (row, col) point"""
    indices_str = [
        el.replace("M", "").replace("Z", "").split(",") for el in path.split("L")
    ]
    return np.array(indices_str, dtype=float)


def largest_connected_component(mask):
    labels, _ = ndimage.label(mask)
    sizes = np.bincount(labels.ravel())[1:]
    return labels == (np.argmax(sizes) + 1)


t2 = time()
print("initial calculations", t2 - t1)

# ------------- Define App Layout ---------------------------------------------------
axial_card = dbc.Card(
    [
        dbc.CardHeader("Axial view of the lung"),
        dbc.CardBody([slicer1.graph, slicer1.slider, *slicer1.stores]),
        dbc.CardFooter(
            [
                html.H6(
                    [
                        "Step 1: Draw a rough outline that encompasses all ground glass occlusions across ",
                        html.Span(
                            "all axial slices",
                            id="tooltip-target-1",
                            className="tooltip-target",
                        ),
                        ".",
                    ]
                ),
                dbc.Tooltip(
                    "Use the slider to scroll vertically through the image and look for the ground glass occlusions.",
                    target="tooltip-target-1",
                ),
            ]
        ),
    ]
)

saggital_card = dbc.Card(
    [
        dbc.CardHeader("Sagittal view of the lung"),
        dbc.CardBody([slicer2.graph, slicer2.slider, *slicer2.stores]),
        dbc.CardFooter(
            [
                html.H6(
                    [
                        "Step 2:\n\nDraw a rectangle to determine the ",
                        html.Span(
                            "min and max height ",
                            id="tooltip-target-2",
                            className="tooltip-target",
                        ),
                        "of the occlusion.",
                    ]
                ),
                dbc.Tooltip(
                    "Only the min and max height of the rectangle are used, the width is ignored",
                    target="tooltip-target-2",
                ),
            ]
        ),
    ]
)

histogram_card = dbc.Card(
    [
        dbc.CardHeader("Histogram of intensity values"),
        dbc.CardBody(
            [
                dcc.Graph(
                    id="graph-histogram",
                    figure=px.bar(
                        x=hi[1],
                        y=hi[0],
                        labels={"x": "intensity", "y": "count"},
                        template="plotly_white",
                    ),
                    config={
                        "modeBarButtonsToAdd": [
                            "drawline",
                            "drawclosedpath",
                            "drawrect",
                            "eraseshape",
                        ]
                    },
                ),
            ]
        ),
        dbc.CardFooter(
            [
                dbc.Toast(
                    [
                        html.P(
                            "Before you can select value ranges in this histogram, you need to define a region"
                            " of interest in the slicer views above (step 1 and 2)!",
                            className="mb-0",
                        )
                    ],
                    id="roi-warning",
                    header="Please select a volume of interest first",
                    icon="danger",
                    is_open=True,
                    dismissable=False,
                ),
                "Step 3: Select a range of values to segment the occlusion. Hover on slices to find the typical "
                "values of the occlusion.",
            ]
        ),
    ]
)

mesh_card = dbc.Card(
    [
        dbc.CardHeader("3D mesh representation of the image data and annotation"),
        dbc.CardBody([dcc.Graph(id="graph-helper", figure=fig_mesh)]),
    ]
)

# Define Modal
with open("assets/modal.md", "r") as f:
    howto_md = f.read()

modal_overlay = dbc.Modal(
    [
        dbc.ModalBody(html.Div([dcc.Markdown(howto_md)], id="howto-md")),
        dbc.ModalFooter(dbc.Button("Close", id="howto-close", className="howto-bn")),
    ],
    id="modal",
    size="lg",
)

# Buttons
button_gh = dbc.Button(
    "Learn more",
    id="howto-open",
    outline=True,
    color="secondary",
    # Turn off lowercase transformation for class .button in stylesheet
    style={"textTransform": "none"},
)

button_howto = dbc.Button(
    "View Code on github",
    outline=True,
    color="primary",
    href="https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-covid-xray",
    id="gh-link",
    style={"text-transform": "none"},
)

nav_bar = dbc.Navbar(
    dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Row(
                            [
                                dbc.Col(
                                    html.A(
                                        html.Img(
                                            src=app.get_asset_url("dash-logo-new.png"),
                                            height="30px",
                                        ),
                                        href="https://plotly.com/dash/",
                                    ),
                                    style={"width": "min-content"},
                                ),
                                dbc.Col(
                                    html.Div(
                                        [
                                            html.H3("Covid X-Ray app"),
                                            html.P(
                                                "Exploration and annotation of CT images"
                                            ),
                                        ],
                                        id="app_title",
                                    )
                                ),
                            ],
                            align="center",
                            style={"display": "inline-flex"},
                        )
                    ),
                    dbc.Col(
                        [
                            dbc.NavbarToggler(id="navbar-toggler"),
                            dbc.Collapse(
                                dbc.Nav(
                                    [dbc.NavItem(button_howto), dbc.NavItem(button_gh)],
                                    className="ml-auto",
                                    navbar=True,
                                ),
                                id="navbar-collapse",
                                navbar=True,
                            ),
                        ]
                    ),
                    modal_overlay,
                ],
                align="center",
                style={"width": "100%"},
            ),
        ],
        fluid=True,
    ),
    color="dark",
    dark=True,
)


app.layout = html.Div(
    [
        nav_bar,
        dbc.Container(
            [
                dbc.Row([dbc.Col(axial_card), dbc.Col(saggital_card)]),
                dbc.Row([dbc.Col(histogram_card), dbc.Col(mesh_card),]),
            ],
            fluid=True,
        ),
        dcc.Store(id="annotations", data={}),
        dcc.Store(id="occlusion-surface", data={}),
    ],
)

t3 = time()
print("layout definition", t3 - t2)


# ------------- Define App Interactivity ---------------------------------------------------
@app.callback(
    [Output("graph-histogram", "figure"), Output("roi-warning", "is_open")],
    [Input("annotations", "data")],
)
def update_histo(annotations):
    if (
        annotations is None
        or annotations.get("x") is None
        or annotations.get("z") is None
    ):
        return dash.no_update, dash.no_update
    # Horizontal mask for the xy plane (z-axis)
    path = path_to_coords(annotations["z"]["path"])
    rr, cc = draw.polygon(path[:, 1] / spacing[1], path[:, 0] / spacing[2])
    if len(rr) == 0 or len(cc) == 0:
        return dash.no_update, dash.no_update
    mask = np.zeros(img.shape[1:])
    mask[rr, cc] = 1
    mask = ndimage.binary_fill_holes(mask)
    # top and bottom, the top is a lower number than the bottom because y values
    # increase moving down the figure
    top, bottom = sorted([int(annotations["x"][c] / spacing[0]) for c in ["y0", "y1"]])
    intensities = med_img[top:bottom, mask].ravel()
    if len(intensities) == 0:
        return dash.no_update, dash.no_update
    hi = exposure.histogram(intensities)
    fig = px.bar(
        x=hi[1],
        y=hi[0],
        # Histogram
        labels={"x": "intensity", "y": "count"},
    )
    fig.update_layout(dragmode="select", title_font=dict(size=20, color="blue"))
    return fig, False


@app.callback(
    [
        Output("occlusion-surface", "data"),
        Output(slicer1.overlay_data.id, "data"),
        Output(slicer2.overlay_data.id, "data"),
    ],
    [Input("graph-histogram", "selectedData"), Input("annotations", "data")],
)
def update_segmentation_slices(selected, annotations):
    ctx = dash.callback_context
    # When shape annotations are changed, reset segmentation visualization
    if (
        ctx.triggered[0]["prop_id"] == "annotations.data"
        or annotations is None
        or annotations.get("x") is None
        or annotations.get("z") is None
    ):
        mask = np.zeros_like(med_img)
        overlay1 = slicer1.create_overlay_data(mask)
        overlay2 = slicer2.create_overlay_data(mask)
        return go.Mesh3d(), overlay1, overlay2
    elif selected is not None and "range" in selected:
        if len(selected["points"]) == 0:
            return dash.no_update
        v_min, v_max = selected["range"]["x"]
        t_start = time()
        # Horizontal mask
        path = path_to_coords(annotations["z"]["path"])
        rr, cc = draw.polygon(path[:, 1] / spacing[1], path[:, 0] / spacing[2])
        mask = np.zeros(img.shape[1:])
        mask[rr, cc] = 1
        mask = ndimage.binary_fill_holes(mask)
        # top and bottom, the top is a lower number than the bottom because y values
        # increase moving down the figure
        top, bottom = sorted(
            [int(annotations["x"][c] / spacing[0]) for c in ["y0", "y1"]]
        )
        img_mask = np.logical_and(med_img > v_min, med_img <= v_max)
        img_mask[:top] = False
        img_mask[bottom:] = False
        img_mask[top:bottom, np.logical_not(mask)] = False
        img_mask = largest_connected_component(img_mask)
        # img_mask_color = mask_to_color(img_mask)
        t_end = time()
        print("build the mask", t_end - t_start)
        t_start = time()
        # Update 3d viz
        verts, faces, _, _ = measure.marching_cubes(
            filters.median(img_mask, selem=np.ones((1, 7, 7))), 0.5, step_size=3
        )
        t_end = time()
        print("marching cubes", t_end - t_start)
        x, y, z = verts.T
        i, j, k = faces.T
        trace = go.Mesh3d(x=z, y=y, z=x, color="red", opacity=0.8, i=k, j=j, k=i)
        overlay1 = slicer1.create_overlay_data(img_mask)
        overlay2 = slicer2.create_overlay_data(img_mask)
        # todo: do we need an output to trigger an update?
        return trace, overlay1, overlay2
    else:
        return (dash.no_update,) * 3


@app.callback(
    Output("annotations", "data"),
    [Input(slicer1.graph.id, "relayoutData"), Input(slicer2.graph.id, "relayoutData"),],
    [State("annotations", "data")],
)
def update_annotations(relayout1, relayout2, annotations):
    if relayout1 is not None and "shapes" in relayout1:
        if len(relayout1["shapes"]) >= 1:
            shape = relayout1["shapes"][-1]
            annotations["z"] = shape
        else:
            annotations.pop("z", None)
    elif relayout1 is not None and "shapes[2].path" in relayout1:
        annotations["z"]["path"] = relayout1["shapes[2].path"]

    if relayout2 is not None and "shapes" in relayout2:
        if len(relayout2["shapes"]) >= 1:
            shape = relayout2["shapes"][-1]
            annotations["x"] = shape
        else:
            annotations.pop("x", None)
    elif relayout2 is not None and (
        "shapes[2].y0" in relayout2 or "shapes[2].y1" in relayout2
    ):
        annotations["x"]["y0"] = relayout2["shapes[2].y0"]
        annotations["x"]["y1"] = relayout2["shapes[2].y1"]
    return annotations


app.clientside_callback(
    """
function(surf, fig){
        let fig_ = {...fig};
        fig_.data[1] = surf;
        return fig_;
    }
""",
    output=Output("graph-helper", "figure"),
    inputs=[Input("occlusion-surface", "data"),],
    state=[State("graph-helper", "figure"),],
)


@app.callback(
    Output("modal", "is_open"),
    [Input("howto-open", "n_clicks"), Input("howto-close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


if __name__ == "__main__":
    app.run_server(debug=True, dev_tools_props_check=False)
