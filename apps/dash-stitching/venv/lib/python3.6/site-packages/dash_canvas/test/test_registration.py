from dash_canvas.utils import register_tiles
import numpy as np
from skimage import data, color
import matplotlib.pyplot as plt

def test_stitching_one_row():
    im = data.moon()
    l = 256

    n_rows = 1
    n_cols = im.shape[1] // l

    init_i, init_j = 0, 0

    imgs = np.empty((n_rows, n_cols, l, l))

    overlap_h = [5, 50]

    i = 0
    for j in range(n_cols):
        sub_im = im[init_i:init_i + l, init_j:init_j + l]
        imgs[i, j] = sub_im
        init_j += l - overlap_h[1]

    stitch = register_tiles(imgs, n_rows, n_cols, overlap_global=0.2, pad=l//2)
    delta = im[:l, :stitch.shape[1]].astype(np.float) - stitch.astype(np.float)
    assert np.all(delta == 0)
    stitch = register_tiles(imgs, n_rows, n_cols, overlap_global=0.2, pad=l//2,
                            blending=False)
    delta = im[:l, :stitch.shape[1]].astype(np.float) - stitch.astype(np.float)
    assert np.all(delta == 0)

    # local_overlap
    stitch = register_tiles(imgs, n_rows, n_cols, overlap_global=0.5, pad=l//2,
                            overlap_local={(0, 1):[0, 45]})
    delta = im[:l, :stitch.shape[1]].astype(np.float) - stitch.astype(np.float)
    assert np.all(delta == 0)


def test_stitching_two_rows():
    im = data.moon()
    l = 256
    # two rows
    n_rows = 2
    n_cols = im.shape[1] // l
    init_i, init_j = 0, 0
    overlap_h = [5, 50]
    overlap_v = 40

    imgs = np.empty((n_rows, n_cols, l, l))
    for i in range(n_rows):
        for j in range(n_cols):
            sub_im = im[init_i:init_i + l, init_j:init_j + l]
            imgs[i, j] = sub_im
            init_j += l - overlap_h[1]
        init_j = 0
        init_i += l - overlap_v

    stitch = register_tiles(imgs, n_rows, n_cols, overlap_global=0.2, pad=l//2)
    delta = im[:stitch.shape[0], :stitch.shape[1]].astype(np.float) - stitch.astype(np.float)
    print(delta.mean())
    assert np.all(delta == 0)


def test_stitching_color():
    im = color.gray2rgb(data.moon())
    l = 256

    n_rows = 1
    n_cols = im.shape[1] // l

    init_i, init_j = 0, 0

    imgs = np.empty((n_rows, n_cols, l, l, 3))

    overlap_h = [5, 50]

    i = 0
    for j in range(n_cols):
        sub_im = im[init_i:init_i + l, init_j:init_j + l]
        imgs[i, j] = sub_im
        init_j += l - overlap_h[1]

    stitch = register_tiles(imgs, n_rows, n_cols, overlap_global=0.2, pad=l//2)
    delta = im[:l, :stitch.shape[1]].astype(np.float) - stitch.astype(np.float)
    assert np.all(delta == 0)

