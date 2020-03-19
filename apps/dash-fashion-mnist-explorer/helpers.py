import io
import base64
import os
import datetime
import gzip
from io import BytesIO

import numpy as np
from PIL import Image, ImageOps, ImageEnhance

label_mapping = {
    0: "T-Shirt / Top",
    1: "Trouser",
    2: "Pullover",
    3: "Dress",
    4: "Coat",
    5: "Sandals",
    6: "Shirt",
    7: "Sneaker",
    8: "Bag",
    9: "Ankle Boots",
}


def load_mnist(path, subset="train"):
    if subset == "train":
        s = "train"
    elif subset == "test":
        s = "t10k"

    labels_path = os.path.join(path, f"{s}-labels-idx1-ubyte.gz")
    images_path = os.path.join(path, f"{s}-images-idx3-ubyte.gz")

    with gzip.open(labels_path, "rb") as lbpath:
        labels = np.frombuffer(lbpath.read(), dtype=np.uint8, offset=8)

    with gzip.open(images_path, "rb") as imgpath:
        images = np.frombuffer(imgpath.read(), dtype=np.uint8, offset=16).reshape(
            len(labels), 784
        )

    return images, labels


def parse_image(contents, filename, date):
    # Take uploaded image, from dcc upload, convert to np array, and reshape
    # for nn
    content_type, content_string = contents.split(",")
    im = Image.open(io.BytesIO(base64.b64decode(content_string)))
    im2 = im.copy()
    resized = np.array(
        ImageOps.invert(ImageOps.fit(im2, (28, 28), Image.ANTIALIAS).convert("L"))
    ).reshape(28, 28, 1)
    return im, resized


def numpy_to_b64(array, scalar=True):
    # Convert from 0-1 to 0-255
    if scalar:
        array = np.uint8(255 * array)

    array[np.where(array == 0)] = 255

    im_pil = Image.fromarray(array)

    enhancer = ImageEnhance.Sharpness(im_pil)
    enhanced_im = enhancer.enhance(10.0)

    buff = BytesIO()
    im_pil.save(buff, format="png")
    im_b64 = base64.b64encode(buff.getvalue()).decode("utf-8")

    return "data:image/png;base64," + im_b64


def create_img(arr, shape=(28, 28)):
    arr = arr.reshape(shape).astype(np.float64)
    image_b64 = numpy_to_b64(arr)
    return image_b64
