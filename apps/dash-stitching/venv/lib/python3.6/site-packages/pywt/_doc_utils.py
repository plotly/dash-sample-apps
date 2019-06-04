"""Utilities used to generate various figures in the documentation."""
from itertools import product

import numpy as np
from matplotlib import pyplot as plt

__all__ = ['wavedec_keys', 'wavedec2_keys', 'draw_2d_wp_basis',
           'draw_2d_fswavedecn_basis', 'pad', 'boundary_mode_subplot']


def wavedec_keys(level):
    """Subband keys corresponding to a wavedec decomposition."""
    approx = ''
    coeffs = {}
    for lev in range(level):
        for k in ['a', 'd']:
            coeffs[approx + k] = None
        approx = 'a' * (lev + 1)
        if lev < level - 1:
            coeffs.pop(approx)
    return list(coeffs.keys())


def wavedec2_keys(level):
    """Subband keys corresponding to a wavedec2 decomposition."""
    approx = ''
    coeffs = {}
    for lev in range(level):
        for k in ['a', 'h', 'v', 'd']:
            coeffs[approx + k] = None
        approx = 'a' * (lev + 1)
        if lev < level - 1:
            coeffs.pop(approx)
    return list(coeffs.keys())


def _box(bl, ur):
    """(x, y) coordinates for the 4 lines making up a rectangular box.

    Parameters
    ==========
    bl : float
        The bottom left corner of the box
    ur : float
        The upper right corner of the box

    Returns
    =======
    coords : 2-tuple
        The first and second elements of the tuple are the x and y coordinates
        of the box.
    """
    xl, xr = bl[0], ur[0]
    yb, yt = bl[1], ur[1]
    box_x = [xl, xr,
             xr, xr,
             xr, xl,
             xl, xl]
    box_y = [yb, yb,
             yb, yt,
             yt, yt,
             yt, yb]
    return (box_x, box_y)


def _2d_wp_basis_coords(shape, keys):
    # Coordinates of the lines to be drawn by draw_2d_wp_basis
    coords = []
    centers = {}  # retain center of boxes for use in labeling
    for key in keys:
        offset_x = offset_y = 0
        for n, char in enumerate(key):
            if char in ['h', 'd']:
                offset_x += shape[0] // 2**(n + 1)
            if char in ['v', 'd']:
                offset_y += shape[1] // 2**(n + 1)
        sx = shape[0] // 2**(n + 1)
        sy = shape[1] // 2**(n + 1)
        xc, yc = _box((offset_x, -offset_y),
                      (offset_x + sx, -offset_y - sy))
        coords.append((xc, yc))
        centers[key] = (offset_x + sx // 2, -offset_y - sy // 2)
    return coords, centers


def draw_2d_wp_basis(shape, keys, fmt='k', plot_kwargs={}, ax=None,
                     label_levels=0):
    """Plot a 2D representation of a WaveletPacket2D basis."""
    coords, centers = _2d_wp_basis_coords(shape, keys)
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig = ax.get_figure()
    for coord in coords:
        ax.plot(coord[0], coord[1], fmt)
    ax.set_axis_off()
    ax.axis('square')
    if label_levels > 0:
        for key, c in centers.items():
            if len(key) <= label_levels:
                ax.text(c[0], c[1], key,
                        horizontalalignment='center',
                        verticalalignment='center')
    return fig, ax


def _2d_fswavedecn_coords(shape, levels):
    coords = []
    centers = {}  # retain center of boxes for use in labeling
    for key in product(wavedec_keys(levels), repeat=2):
        (key0, key1) = key
        offsets = [0, 0]
        widths = list(shape)
        for n0, char in enumerate(key0):
            if char in ['d']:
                offsets[0] += shape[0] // 2**(n0 + 1)
        for n1, char in enumerate(key1):
            if char in ['d']:
                offsets[1] += shape[1] // 2**(n1 + 1)
        widths[0] = shape[0] // 2**(n0 + 1)
        widths[1] = shape[1] // 2**(n1 + 1)
        xc, yc = _box((offsets[0], -offsets[1]),
                      (offsets[0] + widths[0], -offsets[1] - widths[1]))
        coords.append((xc, yc))
        centers[(key0, key1)] = (offsets[0] + widths[0] / 2,
                                 -offsets[1] - widths[1] / 2)
    return coords, centers


def draw_2d_fswavedecn_basis(shape, levels, fmt='k', plot_kwargs={}, ax=None,
                             label_levels=0):
    """Plot a 2D representation of a WaveletPacket2D basis."""
    coords, centers = _2d_fswavedecn_coords(shape, levels)
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig = ax.get_figure()
    for coord in coords:
        ax.plot(coord[0], coord[1], fmt)
    ax.set_axis_off()
    ax.axis('square')
    if label_levels > 0:
        for key, c in centers.items():
            lev = np.max([len(k) for k in key])
            if lev <= label_levels:
                ax.text(c[0], c[1], key,
                        horizontalalignment='center',
                        verticalalignment='center')
    return fig, ax


def pad(x, pad_widths, mode):
    """Extend a 1D signal using a given boundary mode.

    Like numpy.pad but supports all PyWavelets boundary modes.
    """
    if np.isscalar(pad_widths):
        pad_widths = (pad_widths, pad_widths)

    if x.ndim > 1:
        raise ValueError("This padding function is only for 1D signals.")

    if mode in ['symmetric', 'reflect']:
        xp = np.pad(x, pad_widths, mode=mode)
    elif mode in ['periodic', 'periodization']:
        if mode == 'periodization' and x.size % 2 == 1:
            raise ValueError("periodization expects an even length signal.")
        xp = np.pad(x, pad_widths, mode='wrap')
    elif mode == 'zeros':
        xp = np.pad(x, pad_widths, mode='constant', constant_values=0)
    elif mode == 'constant':
        xp = np.pad(x, pad_widths, mode='edge')
    elif mode == 'smooth':
        xp = np.pad(x, pad_widths, mode='linear_ramp',
                    end_values=(x[0] + pad_widths[0] * (x[0] - x[1]),
                                x[-1] + pad_widths[1] * (x[-1] - x[-2])))
    elif mode == 'antisymmetric':
        # implement by flipping portions symmetric padding
        npad_l, npad_r = pad_widths
        xp = np.pad(x, pad_widths, mode='symmetric')
        r_edge = npad_l + x.size - 1
        l_edge = npad_l
        # width of each reflected segment
        seg_width = x.size
        # flip reflected segments on the right of the original signal
        n = 1
        while r_edge <= xp.size:
            segment_slice = slice(r_edge + 1,
                                  min(r_edge + 1 + seg_width, xp.size))
            if n % 2:
                xp[segment_slice] *= -1
            r_edge += seg_width
            n += 1

        # flip reflected segments on the left of the original signal
        n = 1
        while l_edge >= 0:
            segment_slice = slice(max(0, l_edge - seg_width), l_edge)
            if n % 2:
                xp[segment_slice] *= -1
            l_edge -= seg_width
            n += 1
    elif mode == 'antireflect':
        npad_l, npad_r = pad_widths
        # pad with zeros to get final size
        xp = np.pad(x, pad_widths, mode='constant', constant_values=0)

        # right and left edge of the original signal within the padded one
        r_edge = npad_l + x.size - 1
        l_edge = npad_l
        # values of the right and left edge of the original signal
        rv1 = x[-1]
        lv1 = x[0]
        # width of each reflected segment
        seg_width = x.size - 1

        # Generate all reflected segments on the right of the signal.
        # odd reflections of the signal involve these coefficients
        xr_odd = x[-2::-1]
        # even reflections of the signal involve these coefficients
        xr_even = x[1:]
        n = 1
        while r_edge <= xp.size:
            segment_slice = slice(r_edge + 1,
                                  min(r_edge + 1 + seg_width, xp.size))
            orig_sl = slice(segment_slice.stop - segment_slice.start)
            rv = xp[r_edge]
            if n % 2:
                xp[segment_slice] = rv - (xr_odd[orig_sl] - rv1)
            else:
                xp[segment_slice] = rv + (xr_even[orig_sl] - lv1)
            r_edge += seg_width
            n += 1

        # Generate all reflected segments on the left of the signal.
        # odd reflections of the signal involve these coefficients
        xl_odd = x[-1:0:-1]
        # even reflections of the signal involve these coefficients
        xl_even = x[:-1]
        n = 1
        while l_edge >= 0:
            segment_slice = slice(max(0, l_edge - seg_width), l_edge)
            orig_sl = slice(segment_slice.start - segment_slice.stop, None)
            lv = xp[l_edge]
            if n % 2:
                xp[segment_slice] = lv - (xl_odd[orig_sl] - lv1)
            else:
                xp[segment_slice] = lv + (xl_even[orig_sl] - rv1)
            l_edge -= seg_width
            n += 1
    return xp


def boundary_mode_subplot(x, mode, ax, symw=True):
    """Plot an illustration of the boundary mode in a subplot axis."""

    # if odd-length, periodization replicates the last sample to make it even
    if mode == 'periodization' and len(x) % 2 == 1:
        x = np.concatenate((x, (x[-1], )))

    npad = 2 * len(x)
    t = np.arange(len(x) + 2 * npad)
    xp = pad(x, (npad, npad), mode=mode)

    ax.plot(t, xp, 'k.')
    ax.set_title(mode)

    # plot the original signal in red
    if mode == 'periodization':
        ax.plot(t[npad:npad + len(x) - 1], x[:-1], 'r.')
    else:
        ax.plot(t[npad:npad + len(x)], x, 'r.')

    # add vertical bars indicating points of symmetry or boundary extension
    o2 = np.ones(2)
    left = npad
    if symw:
        step = len(x) - 1
        rng = range(-2, 4)
    else:
        left -= 0.5
        step = len(x)
        rng = range(-2, 4)
    if mode in ['smooth', 'constant', 'zeros']:
        rng = range(0, 2)
    for rep in rng:
        ax.plot((left + rep * step) * o2, [xp.min() - .5, xp.max() + .5], 'k-')
