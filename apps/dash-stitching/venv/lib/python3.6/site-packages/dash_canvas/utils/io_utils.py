from PIL import Image
from io import BytesIO
import base64


def array_to_data_url(img, dtype=None):
    """
    Converts numpy array to data string, using Pillow.

    The returned image string has the right format for the ``image_content``
    property of DashCanvas.

    Parameters
    ==========

    img : numpy array

    Returns
    =======

    image_string: str
    """
    if dtype is not None:
        img = img.astype(dtype)
    pil_img = Image.fromarray(img)
    buff = BytesIO()
    pil_img.save(buff, format="png")
    prefix = b'data:image/png;base64,'
    image_string = (prefix + base64.b64encode(buff.getvalue())).decode("utf-8")
    return image_string


def image_string_to_PILImage(image_string):
    """
    Converts image string to PIL image object.
    """
    return Image.open(BytesIO(base64.b64decode(image_string[22:])))
