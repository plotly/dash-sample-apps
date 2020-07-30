import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from skimage import data
from dash_canvas.utils import array_to_data_url

img = data.astronaut()
img = np.tile(img, (2, 2, 1))
img_flip = np.ascontiguousarray(img[::-1])


def make_image_fig(img_array):
    """
    This is the code we would like to get rid of :-)
    """
    height, width, _ = img_array.shape
    img_str = array_to_data_url(img_array)
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
                html.Button("Refresh", id="button2"),
            ],
            style={"width": "45%", "display": "inline-block"},
        ),
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


@app.callback(Output("graph2", "figure"), [Input("button2", "n_clicks")])
def update_other_graph(n):
    if n is None:
        return dash.no_update
    print(n)
    if n % 2 == 0:
        im = img
    else:
        im = img_flip
    return make_image_fig(im)


if __name__ == "__main__":
    app.run_server(debug=True)
