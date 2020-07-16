import skimage.io
import io
import base64


def img_from_b64(b64str):
    imdata = base64.b64decode(b64str)
    im = skimage.io.imread(io.BytesIO(imdata))
    return im


def img_from_mime(mimestr):
    b64str = mimestr.split(",")[-1]
    return img_from_b64(b64str)


def b64_from_img(img, format_str="png"):
    bio = io.BytesIO()
    # use 'pil' plugin because it supports writing image files to memory
    skimage.io.imsave(bio, img, plugin="pil", format_str=format_str)
    b64str = base64.b64encode(bio.getvalue()).decode()
    return b64str


def mime_from_img(img, format_str="png"):
    mimefmt = "data:image/%s;base64" % (format_str,)
    b64str = b64_from_img(img, format_str=format_str)
    mimestr = ",".join([mimefmt, b64str])
    return mimestr


def mime_from_img_path(img_path, format_str="png"):
    img = skimage.io.imread(img_path)
    return mime_from_img(img, format_str=format_str)
