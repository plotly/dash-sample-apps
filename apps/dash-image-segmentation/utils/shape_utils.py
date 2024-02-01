from cairosvg import svg2png
import skimage
import PIL.Image
import io
import numpy as np


def shape_to_svg_code(shape, fig=None, width=None, height=None):
    """
    fig is the plotly.py figure which shape resides in (to get width and height)
    and shape is one of the shapes the figure contains.
    """
    if fig is not None:
        # get width and height
        wrange = next(fig.select_xaxes())["range"]
        hrange = next(fig.select_yaxes())["range"]
        width, height = [max(r) - min(r) for r in [wrange, hrange]]
    else:
        if width is None or height is None:
            raise ValueError("If fig is None, you must specify width and height")
    fmt_dict = dict(
        width=width,
        height=height,
        stroke_color=shape["line"]["color"],
        stroke_width=shape["line"]["width"],
        path=shape["path"],
    )
    return """
<svg
    width="{width}"
    height="{height}"
    viewBox="0 0 {width} {height}"
>
<path
    stroke="{stroke_color}"
    stroke-width="{stroke_width}"
    d="{path}"
    fill-opacity="0"
/>
</svg>
""".format(
        **fmt_dict
    )


def shape_to_png(fig=None, shape=None, width=None, height=None, write_to=None):
    """
    Like svg2png, if write_to is None, returns a bytestring. If it is a path
    to a file it writes to this file and returns None.
    """
    svg_code = shape_to_svg_code(fig=fig, shape=shape, width=width, height=height)
    r = svg2png(bytestring=svg_code, write_to=write_to)
    return r


def shapes_to_mask(shape_args, shape_layers):
    """
    Returns numpy array (type uint8) with number of rows equal to maximum height
    of all shapes's bounding boxes and number of columns equal to their number
    of rows.
    shape_args is a list of dictionaries whose keys are the parameters to the
    shape_to_png function.
    The mask is taken to be all the pixels that are non-zero in the resulting
    image from rendering the shape.
    shape_layers is either a number or an array
    if a number, all the layers have the same number in the mask
    if an array, must be the same length as shape_args and each entry is an
    integer in [0...255] specifying the layer number. Note that the convention
    is that 0 means no mask, so generally the layer numbers will be non-zero.
    """
    images = []
    for sa in shape_args:
        pngbytes = shape_to_png(**sa)
        images.append(PIL.Image.open(io.BytesIO(pngbytes)))

    mwidth, mheight = [max([im.size[i] for im in images]) for i in range(2)]
    mask = np.zeros((mheight, mwidth), dtype=np.uint8)
    if type(shape_layers) != type(list()):
        layer_numbers = [shape_layers for _ in shape_args]
    else:
        layer_numbers = shape_layers
    imarys = []
    for layer_num, im in zip(layer_numbers, images):
        # layer 0 is reserved for no mask
        imary = skimage.util.img_as_ubyte(np.array(im))
        imary = np.sum(imary, axis=2)
        imary.resize((mheight, mwidth))
        imarys.append(imary)
        mask[imary != 0] = layer_num
    return mask
