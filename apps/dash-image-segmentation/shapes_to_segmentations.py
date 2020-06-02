import PIL.Image
import numpy as np
import skimage
import skimage.util
import skimage.io
import skimage.color
import io
import shape_utils
from image_segmentation import trainable_segmentation
import plotly.express as px


def img_to_ubyte_array(img):
    """
    PIL.Image.open is used so that a io.BytesIO object containing the image data
    can be passed as img and parsed into an image. Passing a path to an image
    for img will also work.
    """
    ret_ = skimage.util.img_as_ubyte(np.array(PIL.Image.open(img)))
    return ret_


def fromhex(n):
    return int(n, base=16)


def label_to_colors(
    img, colormap=px.colors.qualitative.Light24, alpha=128, color_class_offset=0
):
    """
    Take MxN matrix containing integers representing labels and return an MxNx4
    matrix where each label has been replaced by a color looked up in colormap.
    colormap entries must be strings like plotly.express style colormaps.
    alpha is the value of the 4th channel
    color_class_offset allows adding a value to the color class index to force
    use of a particular range of colors in the colormap. This is useful for
    example if 0 means 'no class' but we want the color of class 1 to be
    colormap[0].
    """
    colormap = [
        tuple([fromhex(h[s : s + 2]) for s in range(0, len(h), 2)])
        for h in [c.replace("#", "") for c in colormap]
    ]
    cimg = np.zeros(img.shape[:2] + (3,), dtype="uint8")
    minc = np.min(img)
    maxc = np.max(img)
    for c in range(minc, maxc + 1):
        cimg[img == c] = colormap[(c + color_class_offset) % len(colormap)]
    return np.concatenate(
        (cimg, alpha * np.ones(img.shape[:2] + (1,), dtype="uint8")), axis=2
    )


def grey_labels(img):
    minc = np.min(img)
    maxc = np.max(img)
    img -= minc
    img += 1
    img *= 255 // (maxc - minc + 1)
    return img


def compute_segmentations(
    shapes,
    img_path="assets/segmentation_img.jpg",
    segmenter_args={},
    shape_layers=None,
    label_to_colors_args={},
):

    # load original image
    img = img_to_ubyte_array(img_path)

    # load labels
    label_imgs = []
    for shape in shapes:
        lab = img_to_ubyte_array(
            io.BytesIO(
                shape_utils.shape_to_png(
                    shape=shape, width=img.shape[1], height=img.shape[0], write_to=None
                )
            )
        )
        label_imgs.append(lab)

    shape_args = [
        {"width": img.shape[1], "height": img.shape[0], "shape": shape}
        for shape in shapes
    ]
    if (shape_layers is None) or (len(shape_layers) != len(shapes)):
        shape_layers = [(n + 1) for n, _ in enumerate(shapes)]
    mask = shape_utils.shapes_to_mask(shape_args, shape_layers)

    # do segmentation and return this
    seg, clf = trainable_segmentation(img, mask, **segmenter_args)
    color_seg = label_to_colors(seg, **label_to_colors_args)
    # color_seg is a 3d tensor representing a colored image whereas seg is a
    # matrix whose entries represent the classes
    return (color_seg, seg, clf)


def blend_image_and_classified_regions(img, classr):
    """
    If img has an alpha channel, it is ignored.
    If classr has an alpha channel, the images are combined as
        out_img = img * (1 - alpha) + classr * alpha
    If classr doesn't have an alpha channel, just classr is returned.
    Both images are converted to ubyte before blending and the alpha channel is
    divided by 255 to get the scalar.
    The returned image has no alpha channel.
    """
    img = skimage.img_as_ubyte(img)
    classr = skimage.img_as_ubyte(classr)
    img = img[:, :, :3]
    if classr.shape[2] < 4:
        return classr
    alpha = (classr[:, :, 3] / 255)[:, :, None]
    classr = classr[:, :, :3]
    out_img = img * (1 - alpha) + classr * alpha
    out_img = np.round(out_img)
    out_img[out_img > 255] = 255
    out_img[out_img < 0] = 0
    return out_img.astype("uint8")


def blend_image_and_classified_regions_pil(img, classr):
    img = np.array(img)
    classr = np.array(classr)
    out_img = blend_image_and_classified_regions(img, classr)
    return PIL.Image.fromarray(out_img)
