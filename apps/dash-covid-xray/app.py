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

app = dash.Dash(__name__, update_title=None)
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
    dragmode="drawclosedpath", newshape_line_color="cyan"
)
slicer1.graph.config.update(
    modeBarButtonsToAdd=["drawclosedpath", "eraseshape",]
)

slicer2 = VolumeSlicer(app, img, axis=1, spacing=spacing, thumbnail=False)
slicer2.graph.figure.update_layout(dragmode="drawrect", newshape_line_color="cyan")
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


app.layout = html.Div(
    [
        html.Div(
            id="banner",
            children=[
                html.H2(
                    "Exploration and annotation of CT images",
                    id="title",
                    className="seven columns",
                ),
                html.Img(id="logo", src=app.get_asset_url("dash-logo-new.png"),),
            ],
            className="twelve columns app-background",
        ),
        html.Div(
            [
                slicer1.graph,
                slicer1.slider,
                html.H6(
                    [
                        html.Span(
                            "1 - Slide to find the occlusion and draw a path around its contour ❓",
                            id="tooltip-target-1",
                        ),
                    ]
                ),
                dbc.Tooltip(
                    "Draw a rough path which encloses the occlusion at all heights",
                    target="tooltip-target-1",
                ),
                *slicer1.stores,
            ],
            className="app-background",
        ),
        html.Div(
            [
                slicer2.graph,
                slicer2.slider,
                html.H6(
                    [
                        html.Span(
                            "2 - Draw a rectangle to determine the min and max height of the occlusion ❓",
                            id="tooltip-target-2",
                        ),
                    ],
                ),
                dbc.Tooltip(
                    "Only the min and max height of the rectangle are used, the width is ignored",
                    target="tooltip-target-2",
                ),
                *slicer2.stores,
            ],
            className="app-background",
        ),
        html.Div(
            [
                dcc.Graph(
                    id="graph-histogram",
                    figure=px.bar(
                        x=hi[1],
                        y=hi[0],
                        title="Histogram of intensity values - please select first a volume of interest by annotating slices",
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
            ],
            className="app-background",
        ),
        html.Div(
            [dcc.Graph(id="graph-helper", figure=fig_mesh),],
            className="app-background",
        ),
        dcc.Store(id="annotations", data={}),
        dcc.Store(id="occlusion-surface", data={}),
    ],
    className="twelve columns",
)

t3 = time()
print("layout definition", t3 - t2)


@app.callback(
    Output("graph-histogram", "figure"), [Input("annotations", "data")],
)
def update_histo(annotations):
    if (
        annotations is None
        or annotations.get("x") is None
        or annotations.get("z") is None
    ):
        return dash.no_update
    # Horizontal mask for the xy plane (z-axis)
    path = path_to_coords(annotations["z"]["path"])
    rr, cc = draw.polygon(path[:, 1] / spacing[1], path[:, 0] / spacing[2])
    if len(rr) == 0 or len(cc) == 0:
        return dash.no_update
    mask = np.zeros(img.shape[1:])
    mask[rr, cc] = 1
    mask = ndimage.binary_fill_holes(mask)
    # top and bottom, the top is a lower number than the bottom because y values
    # increase moving down the figure
    top, bottom = sorted([int(annotations["x"][c] / spacing[0]) for c in ["y0", "y1"]])
    intensities = med_img[top:bottom, mask].ravel()
    if len(intensities) == 0:
        return dash.no_update
    hi = exposure.histogram(intensities)
    fig = px.bar(
        x=hi[1],
        y=hi[0],
        title="3 - Histogram of intensity values - select a range of values to segment the occlusion <br> Hover on slices to find the typical values of the occlusion",
        labels={"x": "intensity", "y": "count"},
    )
    fig.update_layout(dragmode="select", title_font=dict(size=20, color="blue"))
    return fig


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


if __name__ == "__main__":
    app.run_server(debug=True, dev_tools_props_check=False)
