from dash_canvas.utils import watershed_segmentation, modify_segmentation
from skimage import data, segmentation, morphology, measure
import numpy as np
from scipy import ndimage


def test_watershed_segmentation():
    img = np.zeros((20, 20))
    img[2:6, 2:6] = 1
    img[10:15, 10:15] = 1
    mask = np.zeros_like(img, dtype=np.uint8)
    mask[4, 4] = 1
    mask[12, 12] = 2
    res = watershed_segmentation(img, mask, sigma=0.1)
    assert np.all(res[2:6, 2:6] == 1)
    assert np.all(res[10:15, 10:15] == 2)


def test_split_segmentation():
    img = np.zeros((100, 100), dtype=np.uint8)
    img[:40, 55:] = 1
    img[40:, :30] = 2
    img[40:, 30:65] = 3
    img[40:, 65:] = 4 
    img = ndimage.rotate(img, 20)
    img = measure.label(img)
    img = morphology.remove_small_objects(img, 20)
    img = segmentation.relabel_sequential(img)[0]

    mask = np.zeros_like(img)
    mask[2:53, 75] = 1
    mask[100, 17:60] = 1

    # Labels start at 1
    seg = modify_segmentation(img, measure.label(mask), mode='split')
    assert len(np.unique(seg)) == len(np.unique(img)) + 2

    # Labels start at 0
    seg = modify_segmentation(img + 1, measure.label(mask))
    assert len(np.unique(seg)) == len(np.unique(img)) + 2


def test_merge_segmentation():
    img = np.zeros((20, 20), dtype=np.uint8)
    img[:10, :10] = 1
    img[10:, :10] = 2
    mask = np.zeros_like(img)
    mask[:, 5] = 1
    seg = modify_segmentation(img, mask, mode='merge')
    assert np.all(np.unique(seg) == np.array([0, 1]))
