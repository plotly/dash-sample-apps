import plotly.graph_objects as go
from skimage import io
from PIL import Image
import numpy as np
import os


def make_figure(filename_uri, mode="layout", dragmode="drawrect", show_axes=True):
    if mode == "layout":
        fig = go.Figure()

        # Add trace
        fig.add_trace(go.Scatter(x=[], y=[]))
        filename = os.path.join("assets", os.path.basename(filename_uri))
        im = Image.open(filename)
        width, height = im.size
        # Add images
        fig.add_layout_image(
            dict(
                source=filename_uri,
                xref="x",
                yref="y",
                x=0,
                y=0,
                sizex=width,
                sizey=height,
                sizing="contain",
                layer="below",
            )
        )
        fig.update_layout(template=None)
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
    else:
        im = io.imread(filename_uri[1:])
        fig = go.Figure(go.Image(z=im))
    fig.update_layout(
        margin=dict(t=0, b=0), dragmode=dragmode,
    )
    if not show_axes:
        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=False)
    return fig
