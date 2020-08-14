import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from skimage import data
from PIL import Image
import png
from io import BytesIO
import base64

img = data.retina()
img_flip = np.ascontiguousarray(img[::-1])


def _array_to_b64str(img, backend='pil', compression=4):
    if img.ndim == 2:
        mode = 'L'
    elif img.ndim == 3 and img.shape[-1] == 3:
        mode = 'RGB'
    elif img.ndim == 3 and img.shape[-1] == 4:
        mode = 'RGBA'
    else:
        raise ValueError(
                "Invalid image shape"
                )
    if backend == 'png':
        ndim = img.ndim
        sh = img.shape
        if ndim == 3:
            img = img.reshape((sh[0], sh[1] * sh[2]))
        w = png.Writer(sh[1], sh[0], greyscale=(ndim == 2), compression=compression)
        img_png = png.from_array(img, mode=mode)
        prefix = "data:image/png;base64,"
        with BytesIO() as stream:
            w.write(stream, img_png.rows)
            base64_string = prefix + base64.b64encode(stream.getvalue()).decode("utf-8")
    else:
        pil_img = Image.fromarray(img)
        prefix = "data:image/png;base64,"
        with BytesIO() as stream:
            pil_img.save(stream, format='png', compress_level=compression)
            base64_string = prefix + base64.b64encode(stream.getvalue()).decode("utf-8")
    return base64_string


def make_image_fig(img_array, backend='pil', compression=-1):
    """
    This is the code we would like to get rid of :-)
    """
    height, width, _ = img_array.shape
    img_str = _array_to_b64str(img_array, backend=backend, compression=compression)
    fig = go.Figure(go.Scatter(x=[], y=[]))
    fig.update_layout(template=None)
    fig.add_layout_image(
        dict(
            source=img_str,
            xref="x",
            yref="y",
            x=0,
            y=0,
            sizex=width,
            sizey=height,
            sizing="stretch",
            layer="below",
            opacity=1,
        )
    )
    fig.update_xaxes(
        showgrid=False, showticklabels=False, zeroline=False, range=(0, width)
    )
    fig.update_yaxes(
        showgrid=False,
        scaleanchor="x",
        showticklabels=False,
        zeroline=False,
        range=(height, 0),
    )
    return fig


app = dash.Dash(__name__)

server = app.server

app.layout = html.Div(
    [
        html.Div(
            [
                html.H4("Image trace"),
                dcc.Graph(id="graph", figure=px.imshow(img)),
                html.Button("Refresh", id="button"),
            ],
            style={"width": "45%", "display": "inline-block"},
        ),
        html.Div(
            [
                html.H4("Layout image"),
                dcc.Graph(id="graph2", figure=make_image_fig(img)),
                dcc.Dropdown(
                    id='mode',
                    options=[{'label':'pil', 'value':'pil'}, {'label':'png', 'value':'png'}],
                    value='png'
                    ),
                html.H6('Compression level, between 0 and 9'),
                dcc.Input(
                    id='level',
                    type='number',
                    value=5
                    ),
                html.Button("Refresh", id="button2"),
            ],
            style={"width": "45%", "display": "inline-block"},
        ),
        dcc.Store(id='changes', data=0)
    ]
)


@app.callback(Output("graph", "figure"), [Input("button", "n_clicks")])
def update_graph(n):
    if n is None:
        return dash.no_update
    print(n)
    if n % 2 == 0:
        im = img
    else:
        im = img_flip
    return px.imshow(im)


@app.callback([Output("graph2", "figure"),
               Output("changes", "data")],
            [Input("button2", "n_clicks"),
             Input("level", "value"),
             Input("mode", "value")],
             [State("changes", 'data')])
def update_other_graph(n, compression, backend, changes):
    if n is None and compression is None and backend is None:
        return dash.no_update, dash.no_update
    if changes % 2 == 1:
        im = img
    else:
        im = img_flip
    return make_image_fig(im, compression=compression, backend=backend), changes + 1


if __name__ == "__main__":
    app.run_server(debug=True)
