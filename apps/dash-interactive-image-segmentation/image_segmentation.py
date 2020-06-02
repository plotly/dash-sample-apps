"""
Trainable segmentation using local features and random forests
==============================================================

A pixel-based segmentation is computed here using local features based on
local intensity, edges and textures at different scales. A user-provided
mask is used to identify different regions. The pixels of the mask are used
to train a random-forest classifier [1]_ from scikit-learn. Unlabeled pixels
are then labeled from the prediction of the classifier.

This segmentation algorithm is called trainable segmentation in other software
such as ilastik [2]_ or ImageJ [3]_ (where it is also called "weka
segmentation").

.. [1] https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
.. [2] https://www.ilastik.org/documentation/pixelclassification/pixelclassification
.. [3] https://imagej.net/Trainable_Weka_Segmentation#Training_features_.282D.29

(Based on plot_trainable_segmentation.py from scikit-image)
"""


from itertools import combinations_with_replacement
import itertools
import numpy as np
from skimage import filters, feature, img_as_float32
from sklearn.ensemble import RandomForestClassifier
from time import time


def _features_sigma(img, sigma, intensity=True, edges=True, texture=True):
    """Features for a single value of the Gaussian blurring parameter ``sigma``
    """
    features = []
    img_blur = filters.gaussian(img, sigma)
    if intensity:
        features.append(img_blur)
    if edges:
        features.append(filters.sobel(img_blur))
    if texture:
        H_elems = [
            np.gradient(np.gradient(img_blur)[ax0], axis=ax1)
            for ax0, ax1 in combinations_with_replacement(range(img.ndim), 2)
        ]
        eigvals = feature.hessian_matrix_eigvals(H_elems)
        for eigval_mat in eigvals:
            features.append(eigval_mat)
    return features


def _compute_features_gray(
    img, intensity=True, edges=True, texture=True, sigma_min=0.5, sigma_max=16
):
    """Features for a single channel image. ``img`` can be 2d or 3d.
    """
    # computations are faster as float32
    img = img_as_float32(img)
    sigmas = np.logspace(
        np.log2(sigma_min),
        np.log2(sigma_max),
        num=int(np.log2(sigma_max) - np.log2(sigma_min) + 1),
        base=2,
        endpoint=True,
    )
    n_sigmas = len(sigmas)
    all_results = [
        _features_sigma(img, sigma, intensity=intensity, edges=edges, texture=texture)
        for sigma in sigmas
    ]
    return list(itertools.chain.from_iterable(all_results))


def compute_features(
    img,
    multichannel=True,
    intensity=True,
    edges=True,
    texture=True,
    sigma_min=0.5,
    sigma_max=16,
):
    """Features for a single- or multi-channel image.
    """
    if img.ndim == 3 and multichannel:
        all_results = (
            _compute_features_gray(
                img[..., dim],
                intensity=intensity,
                edges=edges,
                texture=texture,
                sigma_min=sigma_min,
                sigma_max=sigma_max,
            )
            for dim in range(img.shape[-1])
        )
        features = list(itertools.chain.from_iterable(all_results))
    else:
        features = _compute_features_gray(
            img,
            intensity=intensity,
            edges=edges,
            texture=texture,
            sigma_min=sigma_min,
            sigma_max=sigma_max,
        )
    return np.array(features)


def trainable_segmentation(
    img,
    mask=None,
    multichannel=True,
    intensity=True,
    edges=True,
    texture=True,
    sigma_min=0.5,
    sigma_max=16,
    downsample=10,
    clf=None,
    verbose=False,
):
    """
    Segmentation using labeled parts of the image and a random forest classifier.
    """
    t1 = time()
    features = compute_features(
        img,
        multichannel=multichannel,
        intensity=intensity,
        edges=edges,
        texture=texture,
        sigma_min=sigma_min,
        sigma_max=sigma_max,
    )
    t2 = time()
    if clf is None:
        if mask is None:
            raise ValueError("If no classifier clf is passed, you must specify a mask.")
        training_data = features[:, mask > 0].T
        training_labels = mask[mask > 0].ravel()
        data = features[:, mask == 0].T
        t3 = time()
        clf = RandomForestClassifier(n_estimators=100, n_jobs=-1)
        clf.fit(training_data[::downsample], training_labels[::downsample])
        result = np.copy(mask)
    else:
        # we have to flatten all but the first dimension of features
        data = features.reshape((features.shape[0], np.product(features.shape[1:]))).T
        t3 = time()
    t4 = time()
    labels = clf.predict(data)
    if mask is None:
        result = labels.reshape(img.shape[:2])
    else:
        result[mask == 0] = labels
    t5 = time()
    if verbose:
        print("trainable_segmentation timings:")
        print("\tcompute features", t2 - t1)
        print("\tfit", t4 - t3)
        print("\tpredict", t5 - t4)
    return result, clf
