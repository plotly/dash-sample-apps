import plot_common
import image_utils
import numpy as np


def test_image_figure(shape=(300, 500), color="#002EA7"):
    """ Make a figure containing an image that is just a constant color for
    testing. """
    fig = plot_common.dummy_fig()
    im = np.ones(shape, dtype="uint8")
    imc = image_utils.label_to_colors(im, ["#000000", color], alpha=255)
    imcu = plot_common.img_array_to_uri(imc)
    fig = plot_common.add_layout_images_to_fig(fig, [imcu])
    # and we make it so you can draw for fun
    fig.update_layout(
        {
            "dragmode": "drawopenpath",
            "shapes": [],
            "newshape.line.color": "purple",
            "newshape.line.width": 5,
            "margin": dict(l=0, r=0, b=0, t=0, pad=4),
        }
    )
    return fig
