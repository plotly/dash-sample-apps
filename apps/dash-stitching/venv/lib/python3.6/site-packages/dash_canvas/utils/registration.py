import numpy as np
from skimage import io, measure, feature
from scipy import ndimage


def autocrop(img):
    """
    Remove zero-valued rectangles at the border of the image.

    Parameters
    ----------

    img: ndarray
        Image to be cropped
    """
    slices = ndimage.find_objects(img > 0)[0]
    return img[slices]


def _blending_mask(shape):
    mask = np.zeros(shape, dtype=np.int)
    mask[1:-1, 1:-1] = 1
    return ndimage.distance_transform_cdt(mask) + 1


def register_tiles(imgs, n_rows, n_cols, overlap_global=None,
                   overlap_local=None, pad=None, blending=True):
    """
    Stitch together overlapping tiles of a mosaic, using Fourier-based
    registration to estimate the shifts between neighboring tiles.

    Parameters
    ----------

    imgs: array of tiles, of shape (n_rows, n_cols, l_r, l_r) with (l_c, l_r)
        the shape of individual tiles.
    n_rows: int
        number of rows of the mosaic.
    n_cols : int
        number of columns of the mosaic.
    overlap_global : float
        Fraction of overlap between tiles.
    overlap_local : dictionary
        Local overlaps between pairs of tiles. overlap_local[(i, j)] is a pair
        of (x, y) shifts giving the 2D shift vector between tiles i and j.
        Indices (i, j) are the raveled indices of the tile numbers.
    pad : int
        Value of the padding used at the border of the stitched image. An
        autocrop is performed at the end to remove the unnecessary padding.

    Notes
    -----

    Fourier-based registration is used in this function
    (skimage.feature.register_translation).
    """
    if pad is None:
        pad = 200
    l_r, l_c = imgs.shape[2:4]
    if overlap_global is None:
        overlap_global = 0.15
    overlap_value = int(float(overlap_global) * l_r)
    imgs = imgs.astype(np.float)
    if blending:
        blending_mask = _blending_mask((l_r, l_c))
    else:
        blending_mask = np.ones((l_r, l_c))

    if imgs.ndim == 4:
        canvas = np.zeros((2 * pad + n_rows * l_r, 2 * pad + n_cols * l_c),
                      dtype=imgs.dtype)
    else:
        canvas = np.zeros((2 * pad + n_rows * l_r, 2 * pad + n_cols * l_c, 3),
                      dtype=imgs.dtype)
        blending_mask = np.dstack((blending_mask, )*3)
    weights = np.zeros_like(canvas)
    init_r, init_c = pad, pad
    weighted_img = imgs[0, 0] * blending_mask
    canvas[init_r:init_r + l_r, init_c:init_c + l_c] = weighted_img
    weights[init_r:init_r + l_r, init_c:init_c + l_c] = blending_mask
    shifts = np.empty((n_rows, n_cols, 2), dtype=np.int)
    shifts[0, 0] = init_r, init_c

    for i_rows in range(n_rows):
        # Shifts between rows
        if i_rows >= 1:
            index_target = np.ravel_multi_index((i_rows, 0), (n_rows, n_cols))
            index_orig = index_target - n_cols
            try:
                overlap = overlap_local[(index_orig, index_target)]
            except (KeyError, TypeError):
                overlap = np.array([overlap_value, 0])
            init_r, init_c = shifts[i_rows - 1, 0]
            init_r += l_r
            shift_vert = feature.register_translation(
                    imgs[i_rows - 1, 0, -overlap[0]:, :(l_c - overlap[1])],
                    imgs[i_rows, 0, :overlap[0], -(l_c - overlap[1]):])[0]
            init_r += int(shift_vert[0])  - overlap[0]
            init_c += int(shift_vert[1]) - overlap[1]
            shifts[i_rows, 0] = init_r, init_c
            # Fill canvas and weights
            weighted_img = imgs[i_rows, 0] * blending_mask
            canvas[init_r:init_r + l_r, init_c:init_c + l_c] += weighted_img
            weights[init_r:init_r + l_r, init_c:init_c + l_c] += blending_mask
        # Shifts between columns
        for j_cols in range(n_cols - 1):
            index_orig = np.ravel_multi_index((i_rows, j_cols),
                                              (n_rows, n_cols))
            index_target = index_orig + 1
            try:
                overlap = overlap_local[(index_orig, index_target)]
            except (KeyError, TypeError):
                overlap = np.array([0, overlap_value])
            init_c += l_c
            shift_horiz = feature.register_translation(
                imgs[i_rows, j_cols, - (l_r - overlap[0]):, -overlap[1]:],
                imgs[i_rows, j_cols + 1, : l_r - overlap[0], :overlap[1]])[0]
            init_r += int(shift_horiz[0]) + overlap[0]
            init_c += int(shift_horiz[1]) - overlap[1]
            shifts[i_rows, j_cols + 1] = init_r, init_c
            # Fill canvas and weights
            weighted_img = imgs[i_rows, j_cols + 1] * blending_mask
            canvas[init_r:init_r + l_r, init_c:init_c + l_c] += weighted_img
            weights[init_r:init_r + l_r, init_c:init_c + l_c] += blending_mask

    canvas /= (weights + 1.e-5)
    return autocrop(np.rint(canvas).astype(np.uint8))
