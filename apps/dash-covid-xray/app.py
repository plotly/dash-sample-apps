import numpy as np
import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import plotly.graph_objects as go
import plotly.express as px
from dash_canvas.utils import array_to_data_url
from nilearn import image
from skimage import draw, filters, exposure, measure
from scipy import ndimage
from time import time

# ------------- I/O and data massaging ---------------------------------------------------

img = image.load_img("assets/radiopaedia_org_covid-19-pneumonia-7_85703_0-dcm.nii")
mat = img.affine
img = img.get_data()
img_1 = np.copy(np.moveaxis(img, -1, 0))[:, ::-1]
img_2 = np.copy(np.moveaxis(img, -1, 1))

l_h = img.shape[-1]
l_lat = img.shape[0]

size_factor = abs(mat[2, 2] / mat[0, 0])

slices_1 = [array_to_data_url(img_1[i]) for i in range(img_1.shape[0])]
slices_2 = [
    array_to_data_url(img_2[i]) for i in range(img_2.shape[0])[::-1]
]  # vertical


med_img = filters.median(img_1, selem=np.ones((1, 3, 3), dtype=np.bool))
very_filtered_img = filters.median(img_1, selem=np.ones((1, 7, 7), dtype=np.bool))

verts, faces, _, _ = measure.marching_cubes(very_filtered_img, 200, step_size=3)
x, y, z = verts.T
i, j, k = faces.T
fig_mesh = go.Figure()
fig_mesh.add_trace(go.Mesh3d(x=z, y=y, z=x, opacity=0.2, i=k, j=j, k=i))

hi = exposure.histogram(med_img)


def path_to_indices(path):
    """From SVG path to numpy array of coordinates, each row being a (row, col) point
    """
    indices_str = [
        el.replace("M", "").replace("Z", "").split(",") for el in path.split("L")
    ]
    return np.array(indices_str, dtype=float)


def make_figure(
    img, row=None, col=None, size_factor=1, dragmode="drawclosedpath",
):
    """
    Build plotly figure with image slice and visual cues at specific positions.
    """
    fig = px.imshow(img, binary_string=True)
    height, width = img.shape
    height = int(height * size_factor)
    fig.update_traces(dy=size_factor)
    fig.update_layout(
        template=None,
        margin=dict(t=10, b=0, l=0, r=0),
        dragmode=dragmode,
        newshape_line_color="cyan",
        height=400,
    )
    fig.update_xaxes(
        showgrid=False,
        # range=(0, width),
        showticklabels=False,
        zeroline=False,
    )
    fig.update_yaxes(
        showgrid=False,
        scaleanchor="x",
        # range=(height, 0),
        showticklabels=False,
        zeroline=False,
    )
    # Visual cues of position of other slice
    if row is None:
        row = height // 2
    if col is None:
        col = width // 2
    fig.add_shape(
        type="line",
        x0=0,
        x1=width // 20,
        y0=row,
        y1=row,
        line=dict(width=3, color="red"),
        fillcolor="pink",
    )
    fig.add_shape(
        type="line",
        x0=width - width // 20,
        x1=width,
        y0=row,
        y1=row,
        line=dict(width=3, color="red"),
        fillcolor="pink",
    )
    return fig


def largest_connected_component(mask):
    labels, _ = ndimage.label(mask)
    sizes = np.bincount(labels.ravel())[1:]
    return labels == (np.argmax(sizes) + 1)


def mask_to_color(img):
    mask = np.zeros(img.shape + (4,), dtype=np.uint8)
    mask[img] = [255, 0, 0, 100]
    return mask


app = dash.Dash(__name__)

server = app.server

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
                dcc.Graph(
                    id="graph",
                    figure=make_figure(img_1[len(img_1) // 2],),
                    config={
                        "modeBarButtonsToAdd": [
                            "drawline",
                            "drawclosedpath",
                            "drawrect",
                            "eraseshape",
                        ],
                        "scrollZoom": True,
                    },
                ),
                dcc.Slider(
                    id="slider",
                    min=0,
                    max=len(img_1) - 1,
                    step=1,
                    value=len(img_1) // 2,
                    updatemode="drag",
                ),
                html.H6(
                    "1 - Slide to find the occlusion and draw a path around its contour"
                ),
            ],
            className="app-background",
        ),
        html.Div(
            [
                dcc.Graph(
                    id="graph-2",
                    figure=make_figure(
                        img_2[len(img_2) // 2],
                        # width=630,
                        # height=45 * size_factor,
                        dragmode="drawrect",
                        size_factor=size_factor,
                    ),
                    config={"scrollZoom": True},
                ),
                dcc.Slider(
                    id="slider-2",
                    min=0,
                    max=len(img_2) - 1,
                    step=1,
                    value=len(img_2) // 2,
                    updatemode="drag",
                ),
                html.H6(
                    "2 - Draw a rectangle to determine the min and max height of the occlusion"
                ),
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
        dcc.Store(id="small-slices", data=slices_1),
        dcc.Store(id="small-slices-2", data=slices_2),
        dcc.Store(id="z-pos", data=l_h // 2),
        dcc.Store(id="x-pos", data=l_lat // 2),
        dcc.Store(id="annotations", data={}),
        dcc.Store(id="segmentation-slices", data={}),
        dcc.Store(id="segmentation-slices-2", data={}),
    ],
    className="twelve columns",
)


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
    # Horizontal mask
    path = path_to_indices(annotations["z"]["path"])
    rr, cc = draw.polygon(path[:, 1], path[:, 0])
    mask = np.zeros((l_lat, l_lat))
    mask[rr, cc] = 1
    mask = ndimage.binary_fill_holes(mask)
    # top and bottom
    top = int(annotations["x"]["y0"] / size_factor)
    bottom = int(annotations["x"]["y1"] / size_factor)
    intensities = med_img[top:bottom, mask].ravel()
    hi = exposure.histogram(intensities)
    fig = px.bar(
        x=hi[1],
        y=hi[0],
        title="3 - Histogram of intensity values - select a range of values to segment the occlusion",
        labels={"x": "intensity", "y": "count"},
    )
    fig.update_layout(dragmode="select")
    fig.update_layout(images=[])
    return fig


@app.callback(
    [
        Output("segmentation-slices", "data"),
        Output("segmentation-slices-2", "data"),
        Output("graph-helper", "figure"),
    ],
    [Input("graph-histogram", "selectedData")],
    [State("annotations", "data"), State("graph-helper", "figure")],
)
def update_segmentation_slices(selected, annotations, fig_mesh3d):
    if (
        annotations is None
        or annotations.get("x") is None
        or annotations.get("z") is None
    ):
        return (dash.no_update,) * 3
    if selected is not None and "range" in selected:
        v_min, v_max = selected["range"]["x"]
        t_start = time()
        img_mask = np.logical_and(med_img > v_min, med_img <= v_max)
        # Horizontal mask
        path = path_to_indices(annotations["z"]["path"])
        rr, cc = draw.polygon(path[:, 1], path[:, 0])
        mask = np.zeros((l_lat, l_lat))
        mask[rr, cc] = 1
        mask = ndimage.binary_fill_holes(mask)
        # top and bottom
        top = int(annotations["x"]["y0"] / size_factor)
        bottom = int(annotations["x"]["y1"] / size_factor)
        img_mask = np.logical_and(med_img > v_min, med_img <= v_max)
        img_mask[:top] = False
        img_mask[bottom:] = False
        img_mask[top:bottom, np.logical_not(mask)] = False
        img_mask = largest_connected_component(img_mask)
        img_mask_color = mask_to_color(img_mask)
        t_end = time()
        print("build the mask", t_end - t_start)
        t_start = time()
        # Update 3d viz
        verts, faces, _, _ = measure.marching_cubes(img_mask, 0.5, step_size=3)
        t_end = time()
        print("marching cubes", t_end - t_start)
        x, y, z = verts.T
        i, j, k = faces.T
        fig = go.Figure(fig_mesh3d)
        trace = go.Mesh3d(x=z, y=y, z=x, color="red", opacity=0.8, i=k, j=j, k=i)
        if len(fig.data) > 1:
            fig.data[1]["x"] = z
            fig.data[1]["y"] = y
            fig.data[1]["z"] = x
            fig.data[1]["i"] = k
            fig.data[1]["j"] = j
            fig.data[1]["k"] = i
        else:
            fig.add_trace(trace)
        t_start = time()
        # Build lists of binary strings for segmented slices
        slices_dict = {
            str(i): array_to_data_url(img_mask_color[i]) for i in range(top, bottom)
        }
        slices_dict_2 = {
            str(i): array_to_data_url(img_mask_color[:, i])
            for i in range(img_mask_color.shape[1])
        }
        t_end = time()
        print("binary string", t_end - t_start)
        return (slices_dict, slices_dict_2, fig)
    else:
        return (dash.no_update,) * 3


@app.callback(
    Output("annotations", "data"),
    [Input("graph", "relayoutData"), Input("graph-2", "relayoutData"),],
    [State("annotations", "data")],
)
def update_store(relayout, relayout2, annotations):
    if relayout is None:
        return dash.no_update
    if relayout is not None and "shapes" in relayout and len(relayout["shapes"]) > 1:
        shape = relayout["shapes"][-1]
        annotations["z"] = shape
    if relayout2 is not None and "shapes" in relayout2 and len(relayout2["shapes"]) > 1:
        shape = relayout2["shapes"][-1]
        annotations["x"] = shape
    return annotations


app.clientside_callback(
    """
function(n_slider, n_slider_2, seg_slices, slices, fig, annotations){
        var size_factor = 8.777;
        zpos = n_slider;
        xpos = n_slider_2;
        let fig_ = {...fig};
        fig_.data[0].source = slices[zpos];
        if (fig_.data.length  == 1 && zpos in seg_slices){
            fig_.data.push({...fig.data[0]});
            fig_.data[1].source = seg_slices[zpos];
            fig_.data[1].hoverinfo = 'skip';
            fig_.data[1].hovertemplate = '';
        }
        if (fig_.data.length  > 1 && zpos in seg_slices){
            fig_.data[1].source = seg_slices[zpos];
            fig_.data[1].hoverinfo = 'skip';
            fig_.data[1].hovertemplate = '';
        }
        fig_.layout.shapes[0].y0 = xpos;
        fig_.layout.shapes[0].y1 = xpos;
        fig_.layout.shapes[1].y0 = xpos;
        fig_.layout.shapes[1].y1 = xpos;
        if ('z' in annotations){
            fig_.layout.shapes.push(annotations['z']);
            }
        return fig_;
    }
""",
    output=Output("graph", "figure"),
    inputs=[
        Input("slider", "value"),
        Input("slider-2", "value"),
        Input("segmentation-slices", "data"),
    ],
    state=[
        State("small-slices", "data"),
        State("graph", "figure"),
        State("annotations", "data"),
    ],
)


app.clientside_callback(
    """
function(n_slider, n_slider_2, seg_slices_2, slices_2, fig_2, annotations){
        var size_factor = 8.777;
        zpos = n_slider;
        xpos = n_slider_2;
        let fig_2_ = {...fig_2};
        fig_2_.data[0].source = slices_2[xpos];
        if (fig_2_.data.length  == 1 && xpos in seg_slices_2){
            fig_2_.data.push({...fig_2.data[0]});
            fig_2_.data[1].source = seg_slices_2[xpos];
            fig_2_.data[1].hoverinfo = 'skip';
            fig_2_.data[1].hovertemplate = '';
        }
        if (fig_2_.data.length  > 1 && zpos in seg_slices_2){
            fig_2_.data[1].source = seg_slices_2[xpos];
            fig_2_.data[1].hoverinfo = 'skip';
            fig_2_.data[1].hovertemplate = '';
        }
        fig_2_.layout.shapes[0].y0 = zpos * size_factor;
        fig_2_.layout.shapes[0].y1 = zpos * size_factor;
        fig_2_.layout.shapes[1].y0 = zpos * size_factor;
        fig_2_.layout.shapes[1].y1 = zpos * size_factor;
        if ('x' in annotations){
            fig_2_.layout.shapes.push(annotations['x']);
            }
        return fig_2_;
    }
""",
    output=Output("graph-2", "figure"),
    inputs=[
        Input("slider", "value"),
        Input("slider-2", "value"),
        Input("segmentation-slices-2", "data"),
    ],
    state=[
        State("small-slices-2", "data"),
        State("graph-2", "figure"),
        State("annotations", "data"),
    ],
)


if __name__ == "__main__":
    app.run_server(debug=True)
