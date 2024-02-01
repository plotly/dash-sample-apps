import base64
import PIL.Image
import io
import numpy as np
import plotly.graph_objects as go
import skimage.util
from plotly.utils import ImageUriValidator


def path_to_img_ndarray(path):
    with open(path, "rb") as fp:
        img = PIL.Image.open(fp)
        img_array = np.array(img)
        return skimage.util.img_as_float(img_array)


def str_to_pil_img(s):
    return PIL.Image.open(io.BytesIO(s))


def str_to_img_ndarrary(s):
    """
    Returns image in ndarray.
    ndarray will have dimensions (width,height,number_of_color_channels)
    This format is compatible with skimage
    """
    img = str_to_pil_img(s)
    img_array = np.array(img)
    return skimage.util.img_as_float(img_array)


def base64_to_img_array(data):
    decoded_img = base64.b64decode(data)
    return str_to_img_ndarrary(decoded_img)


def str_to_img_trace(s):
    img_array = str_to_img_ndarrary(s)
    img_trace = go.Image(z=img_array)
    return img_trace


def base64_to_img_trace(data):
    decoded_img = base64.b64decode(data)
    return str_to_img_trace(decoded_img)


def dummy_fig():
    fig = go.Figure(go.Scatter(x=[], y=[]))
    fig.update_layout(template=None)
    fig.update_xaxes(showgrid=False, showticklabels=False, zeroline=False)
    fig.update_yaxes(
        showgrid=False, scaleanchor="x", showticklabels=False, zeroline=False
    )
    return fig


def _pilim_if_path(im):
    if type(im) == type(str()):
        return PIL.Image.open(im)
    return im


def add_layout_images_to_fig(fig, images, update_ranges=True):
    """ images is a sequence of PIL Image objects """
    if len(images) <= 0:
        return fig
    for im in images:
        # if image is a path to an image, load the image to get its size
        width, height = _pilim_if_path(im).size
        # Add images
        fig.add_layout_image(
            dict(
                source=im,
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
    if update_ranges:
        width, height = [
            max([_pilim_if_path(im).size[i] for im in images]) for i in range(2)
        ]
        # TODO showgrid,showticklabels,zeroline should be passable to this
        # function
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
    return fig


def img_array_to_pil_image(ia):
    ia = skimage.util.img_as_ubyte(ia)
    img = PIL.Image.fromarray(ia)
    return img


def pil_image_to_layout_image_fig(img):
    width, height = img.size

    fig = dummy_fig()
    fig.add_layout_image(
        dict(
            source=img,
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
    fig.update_xaxes(showgrid=False, range=(0, width))
    fig.update_yaxes(showgrid=False, range=(height, 0))
    return fig


def img_array_to_layout_image_fig(ia):
    """ Returns a figure containing a layout image for faster rendering in the browser. """
    img = img_array_to_pil_image(ia)
    return pil_image_to_layout_image_fig(img)


def pil_image_to_uri(img):
    return ImageUriValidator.pil_image_to_uri(img)


def img_array_to_uri(img_array):
    imgf = img_array_to_pil_image(img_array)
    uri = pil_image_to_uri(imgf)
    return uri


def img_array_to_mime_bytes(img_array):
    uri = img_array_to_uri(img_array)
    mime, contents = uri.split(";")
    typ, cont = contents.split(",")
    byt = base64.b64decode(cont)
    return (mime, byt)
