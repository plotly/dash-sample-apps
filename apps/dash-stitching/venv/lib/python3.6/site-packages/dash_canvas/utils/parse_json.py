import numpy as np
import json
from skimage import draw, morphology
from scipy import ndimage


def _indices_of_path(path, scale=1):
    """
    Retrieve pixel indices (integer values). 

    Parameters
    ----------

    path: SVG-like path formatted as JSON string
        The path is formatted like
        ['M', x0, y0],
        ['Q', xc1, yc1, xe1, ye1],
        ['Q', xc2, yc2, xe2, ye2],
        ...
        ['L', xn, yn]
        where (xc, yc) are for control points and (xe, ye) for end points.

    Notes
    -----

    I took a weight of 1 and it seems fine from visual inspection.
    """
    rr, cc = [], []
    for (Q1, Q2) in zip(path[:-2], path[1:-1]):
        # int(round()) is for Python 2 compatibility
        inds = draw.bezier_curve(int(round(Q1[-1] / scale)), 
                                 int(round(Q1[-2] / scale)), 
                                 int(round(Q2[2] / scale)), 
                                 int(round(Q2[1] / scale)), 
                                 int(round(Q2[4] / scale)), 
                                 int(round(Q2[3] / scale)), 1)
        rr += list(inds[0])
        cc += list(inds[1])
    return rr, cc


def parse_jsonstring(string, shape=None, scale=1):
    """
    Parse JSON string to draw the path saved by react-sketch.

    Up to now only path objects are processed (created with Pencil tool).

    Parameters
    ----------

    data : str
        JSON string of data
    shape: tuple, optional
        shape of returned image.

    Returns
    -------

    mask: ndarray of bools
        binary array where the painted regions are one.
    """
    if shape is None:
        shape = (500, 500)
    mask = np.zeros(shape, dtype=np.bool)
    try:
        data = json.loads(string)
    except:
        return mask
    scale = 1
    for obj in data['objects']:
        if obj['type'] == 'image':
            scale = obj['scaleX']
        elif obj['type'] == 'path':
            scale_obj = obj['scaleX']
            inds = _indices_of_path(obj['path'], scale=scale / scale_obj)
            radius = round(obj['strokeWidth'] / 2. / scale)
            mask_tmp = np.zeros(shape, dtype=np.bool)
            mask_tmp[inds[0], inds[1]] = 1
            mask_tmp = ndimage.binary_dilation(mask_tmp,
                                                  morphology.disk(radius))
            mask += mask_tmp
    return mask


def parse_jsonstring_line(string):
    """
    Return geometry of line objects.

    Parameters
    ----------

    data : str
        JSON string of data

    """
    try:
        data = json.loads(string)
    except:
        return None
    scale = 1
    props = []
    for obj in data['objects']:
        if obj['type'] == 'image':
            scale = obj['scaleX']
        elif obj['type'] == 'line':
            length = np.sqrt(obj['width']**2 + obj['height']**2)
            scale_factor = obj['scaleX'] / scale
            props.append([scale_factor * length,
                          scale_factor * obj['width'],
                          scale_factor * obj['height'],
                          scale_factor * obj['left'],
                          scale_factor * obj['top']])
    return (np.array(props)).astype(np.int)


def parse_jsonstring_rectangle(string):
    """
    Return geometry of rectangle objects.

    Parameters
    ----------

    data : str
        JSON string of data

    """
    try:
        data = json.loads(string)
    except:
        return None
    scale = 1
    props = []
    for obj in data['objects']:
        if obj['type'] == 'image':
            scale = obj['scaleX']
        elif obj['type'] == 'rect':
            scale_factor = obj['scaleX'] / scale
            props.append([scale_factor * obj['width'],
                          scale_factor * obj['height'],
                          scale_factor * obj['left'],
                          scale_factor * obj['top']])
    return (np.array(props)).astype(np.int)


def parse_jsonfile(filename, shape=None):
    with open(filename, 'r') as fp:
        string = json.load(fp)
    return parse_jsonstring(string, shape=shape)


