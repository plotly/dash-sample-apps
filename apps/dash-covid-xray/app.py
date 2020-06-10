import numpy as np
import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import plotly.graph_objects as go
from dash_canvas.utils import array_to_data_url
from nilearn import image
from skimage import draw, segmentation
from skimage import img_as_ubyte

img = image.load_img("assets/radiopaedia_org_covid-19-pneumonia-7_85703_0-dcm.nii")
mat = img.affine
img = img.get_data()
img_1 = np.copy(np.moveaxis(img, -1, 0))
img_2 = np.copy(np.moveaxis(img, -1, 1))

l_h = img.shape[-1]
l_lat = img.shape[0]

size_factor = abs(mat[2, 2] / mat[0, 0])

slices_1 = [array_to_data_url(img_1[i]) for i in range(img_1.shape[0])]
slices_2 = [array_to_data_url(img_2[i]) for i in range(img_2.shape[0])]  # vertical


def path_to_indices(path):
    """From SVG path to numpy array of coordinates, each row being a (row, col) point
    """
    indices_str = [
        el.replace("M", "").replace("Z", "").split(",") for el in path.split("L")
    ]
    return np.array(indices_str, dtype=float)


def make_figure(filename_uri, width, height, row=None, col=None, size_factor=1):
    fig = go.Figure()
    # Add trace
    fig.add_trace(go.Scatter(x=[], y=[]))
    # Add images
    fig.add_layout_image(
        dict(
            source=filename_uri,
            xref="x",
            yref="y",
            x=0,
            y=0,
            sizex=width,
            sizey=height * size_factor,
            sizing="stretch",
            layer="below",
            opacity=1,
        )
    )
    fig.update_layout(template=None, margin=dict(t=10, b=0, l=0, r=0))
    fig.update_xaxes(
        showgrid=False, range=(0, width), showticklabels=False, zeroline=False
    )
    fig.update_yaxes(
        showgrid=False,
        scaleanchor="x",
        range=(height, 0),
        showticklabels=False,
        zeroline=False,
    )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)
    if row is None:
        row = height // 2
    if col is None:
        col = width // 2
    fig.add_shape(
        type="line",
        x0=width // 2 - 50,
        x1=width // 2 + 50,
        y0=row,
        y1=row,
        line=dict(width=5, color="red"),
        fillcolor="pink",
    )
    fig.update_layout(dragmode="drawclosedpath", newshape_line_color="cyan", height=400)
    return fig


app = dash.Dash(__name__)

server = app.server

app.layout = html.Div(
    [        html.Div(
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
                    figure=make_figure(
                        slices_1[len(img_1) // 2], width=630, height=630
                    ),
                ),
                dcc.Slider(
                    id="slider",
                    min=0,
                    max=len(img_1) - 1,
                    step=1,
                    value=len(img_1) // 2,
                    updatemode="drag",
                ),
                html.H6(id="slider-display"),
                html.Button("interpolate", id="interp-button"),
                html.H6(id="volume-display"),
            ],
            className="app-background",
            #style={"width": "50%", "display": "inline-block"},
        ),
        html.Div(
            [
                dcc.Graph(
                    id="graph-2",
                    figure=make_figure(
                        slices_2[len(img_2) // 2], width=630, height=45 * size_factor
                    ),
                ),
                dcc.Slider(
                    id="slider-2",
                    min=0,
                    max=len(img_2) - 1,
                    step=1,
                    value=len(img_2) // 2,
                    updatemode="drag",
                ),
                html.H6(id="slider-2-display"),
            ],
            className="app-background",
            #style={"width": "50%", "display": "inline-block"},
        ),
        dcc.Store(id="small-slices", data=slices_1),
        dcc.Store(id="small-slices-2", data=slices_2),
        dcc.Store(id="z-pos", data=l_h // 2),
        dcc.Store(id="x-pos", data=l_lat // 2),
        dcc.Store(id="annotations", data={}),
        dcc.Store(id="interpolated-annotations", data={}),
        dcc.Store(id="last-shape", data={}),
    ],
    className="twelve columns"
)


@app.callback(
    [Output("interpolated-annotations", "data"),
        Output("volume-display", "children")],
    [Input("interp-button", "n_clicks")],
    [State("annotations", "data"),],
)
def interpolate(n_clicks, annotations):
    if n_clicks is None:
        return dash.no_update
    interps = {}
    keys = np.sort([int(key) for key in annotations])
    z_min = keys[0]
    z_max = keys[-1]
    volume = 0
    for z in np.arange(z_min + 1, z_max, dtype=np.uint8):
        i = np.searchsorted(keys, z)
        t = (z - keys[i - 1]) / (keys[i] - keys[i - 1])
        path1 = path_to_indices(annotations[str(keys[i])]["path"])
        rr, cc = draw.polygon(path1[:, 1], path1[:, 0])
        poly1 = np.zeros((l_lat, l_lat))
        poly1[rr, cc] = 1
        path2 = path_to_indices(annotations[str(keys[i - 1])]["path"])
        rr, cc = draw.polygon(path2[:, 1], path2[:, 0])
        poly2 = np.zeros((l_lat, l_lat))
        poly2[rr, cc] = 1
        poly = (t * poly1 + (1 - t) * poly2) > 0.5
        interps[str(z)] = array_to_data_url(img_as_ubyte(segmentation.mark_boundaries(img_1[z], poly, mode='thick')))
        volume += poly.sum()
    volume *= np.abs(np.linalg.det(mat)) / 1000
    result = f"The volume of the occlusion is {volume:.0f} cm3"
    return interps, result


@app.callback(
    [Output("annotations", "data"), Output("last-shape", "data"),],
    [Input("graph", "relayoutData"), Input("slider", "value")],
    [State("annotations", "data"), State("last-shape", "data"),],
)
def update_store(relayout, z_pos, annotations, last_shape):
    if relayout is None:
        return (dash.no_update,) * 2
    if (
        relayout is not None
        and "shapes" in relayout
        and len(relayout["shapes"]) > 1
        and not (relayout["shapes"][-1]["path"] == last_shape)
    ):
        shape = relayout["shapes"][-1]
        last_shape = relayout["shapes"][-1]["path"]
        annotations[str(z_pos)] = shape
    print("relayout", relayout.keys(), z_pos, relayout)
    return annotations, last_shape


app.clientside_callback(
    """
function(n_slider, n_slider_2, slices, fig, annotations, interps){
        var size_factor = 8.777;
        zpos = n_slider;
        xpos = n_slider_2;
        let fig_ = {...fig};
        console.log("ok 1");
        fig_.layout.images[0].source = slices[zpos];
        fig_.layout.shapes = [fig.layout.shapes[0]];
        fig_.layout.shapes[0].y0 = xpos;
        fig_.layout.shapes[0].y1 = xpos;
        if (n_slider.toString() in interps){
            console.log("in interp range");
            fig_.layout.images[0].source = interps[n_slider.toString()];
        }
        if (n_slider.toString() in annotations){
            console.log("bla");
            fig_.layout.shapes.push(annotations[n_slider.toString()]);
        }
        console.log(fig_);
        return fig_;
    }
""",
    output=Output("graph", "figure"),
    inputs=[Input("slider", "value"), Input("slider-2", "value"),],
    state=[
        State("small-slices", "data"),
        State("graph", "figure"),
        State("annotations", "data"),
        State("interpolated-annotations", "data"),
    ],
)

app.clientside_callback(
    """
function(n_slider, n_slider_2, slices_2, fig_2){
        var size_factor = 8.777;
        zpos = n_slider;
        xpos = n_slider_2;
        let fig_2_ = {...fig_2};
        fig_2_.layout.images[0].source = slices_2[xpos];
        fig_2_.layout.shapes[0].y0 = zpos * size_factor;
        fig_2_.layout.shapes[0].y1 = zpos * size_factor;
        return fig_2_;
    }
""",
    output=Output("graph-2", "figure"),
    inputs=[Input("slider", "value"), Input("slider-2", "value"),],
    state=[State("small-slices-2", "data"), State("graph-2", "figure"),],
)


if __name__ == "__main__":
    app.run_server(debug=True)
