import numpy as np
from skimage import segmentation, morphology, measure, color, feature, filters
from skimage import img_as_float
from scipy import ndimage
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

# ------------------- Modification of segmentation ----------------------

def _split_labels(labels, mask, img, erosion_width=10, compactness=0):
    """
    Divide already labeled array ``labels`` according to array of annotations
    ``mask``.

    Parameters
    ----------

    labels : array of ints
        Array of labels.
    mask : array of ints
        Array with annotations.
    img : array, default None
        Image used for the segmentation.
    erosion_width : int, optional
        Number of pixels to erode on both sides of annotations. Decrease if
        boundaries should correspond to the exact annotations, increase if
        boundaries should rather follow strong image gradients.
    compactness : float, optional (default 0)
        Parameter of the compact waterhsed algorithm
        ``skimage.segmentation.watershed``, should be zero (for normal
        watershed) or greater than zero (for compact watershed).
    """
    out = np.copy(labels)
    bounding_boxes = ndimage.find_objects(labels)
    max_label = labels.max()
    annot_indices = np.unique(mask)[1:]
    count = max_label + 1
    shift = 0 if labels.min() > 0 else 1
    for annot_index in annot_indices:
        obj_label = (np.argmax(np.bincount(labels[mask == annot_index])[shift:]
                                )  + shift)
        # Get subarrays
        box = bounding_boxes[obj_label - 1]
        img_box = img[box]
        labels_box = labels[box]
        # Prepare watershed
        gradient_img = - ndimage.gaussian_gradient_magnitude(img_box, 2)
        mask_box = np.ones(img_box.shape, dtype=np.uint8)
        mask_box[mask[box] == annot_index] = 0
        mask_box = morphology.binary_erosion(mask_box,
                                             morphology.disk(erosion_width))
        masked_region = labels_box == obj_label
        mask_box[np.logical_not(masked_region)] = 0
        mask_box = measure.label(mask_box)
        res = segmentation.watershed(gradient_img, mask_box,
                                     compactness=compactness,
                                     mask=masked_region)
        out[box][res == 1] = count  # only modify one of the regions
        count += 1
    return out


def _merge_labels(labels, mask, skip_zero=True):
    """
    Merge objects covered by the same annotation, given an array of 
    labels defining objects and a labeled mask of annotations.

    Parameters
    ----------

    labels : array of ints
        Array of labels.
    mask : array of ints
        Array with annotations.

    """
    annot_indices = np.unique(mask)[1:]
    label_indices = np.arange(labels.max() + 1, dtype=np.int)
    for index in annot_indices:
        object_indices = np.unique(labels[mask == index])
        if skip_zero:
            object_indices = np.setdiff1d(object_indices, [0])
        label_indices[object_indices] = object_indices[0]
    relabeled, _, _ = segmentation.relabel_sequential(label_indices)
    return relabeled[labels]


def modify_segmentation(labels, mask, img=None, mode='split'):
    """
    Modify already labeled image according to annotations. In 'split' mode,
    each annotation is used to divide a label in two subregions. In 'merge'
    mode, objects covered by the same annotation are merged together.

    Parameters
    ----------

    labels : array of ints
        Array of labels.
    mask : binary of int array
        Array with annotations.
    img : array, default None
        Image used for the segmentation.
    mode: string, 'split' or 'merge'

    Returns
    -------

    out : array of ints
        New labels.
    """
    labels = np.asarray(labels)
    mask = measure.label(mask)
    if img is None:
        img = np.zeros_like(mask)
    
    if mode == 'split':
        return _split_labels(labels, mask, img)
    elif mode == 'merge':
        return _merge_labels(labels, mask)
    else:
        raise ValueError('mode should be either split or merge')


# ------------- Segmentation from markers ----------------------------

def watershed_segmentation(img, mask, sigma=4):
    """
    Watershed segmentation of image using annotations as label markers.

    The image used for the watershed is minus the gradient of the original
    image, convoluted with a Gaussian for more robustness.

    Parameters
    ----------

    img : ndarray
        image to be segmented
    mask : ndarray of ints
        binary array, each connected component corresponds to a different
        object to be segmented
    sigma : float
        standard deviation of Gaussian convoluting the gradient. Increase
        for smoother boundaries.
    """
    if img.ndim > 2:
        img = color.rgb2gray(img)
    labels = measure.label(mask)
    gradient_img = - ndimage.gaussian_gradient_magnitude(img, sigma)
    output = segmentation.watershed(gradient_img, labels)
    return output


def random_walker_segmentation(img, mask, beta=5000):
    """
    Random-walker segmentation of image using annotations as label markers.

    Parameters
    ----------

    img : ndarray
        image to be segmented
    mask : ndarray of ints
        binary array, each connected component corresponds to a different
        object to be segmented
    beta : float
        beta parameter of random walker, smaller beta results in more compact
        objects, larger beta in boundaries following gradients more
        accurately.

    Returns
    -------

    output : ndarray of ints
        image of labels

    Notes
    -----

    see skimage.segmentation.random_walker
    """
    if img.ndim > 2:
        img = color.rgb2gray(img)
    labels = measure.label(mask)
    output = segmentation.random_walker(img, labels, beta=beta, mode='cg_mg')
    return output


def _compute_features_gabor(im):
    gabor_frequencies = np.logspace(-5, 0.5, num=6, base=2)
    thetas = [0, np.pi/4., np.pi/2]
    nb_fq = len(gabor_frequencies) * len(thetas)
    im = np.atleast_3d(im)
    im_gabor = np.empty((im.shape[-1], nb_fq) + im.shape[:2])
    for ch in range(im.shape[-1]):
        img = img_as_float(im[..., ch])
        for i_fq, fq in enumerate(gabor_frequencies):
            for i_th, theta in enumerate(thetas):
                tmp = filters.gabor(img, fq, theta=theta)
                im_gabor[ch, len(thetas) * i_fq + i_th] = \
                                    np.abs(tmp[0] + 1j * tmp[1])
    return im_gabor


def random_forest_segmentation(img, mask, mode='gabor'):
    """
    Segmentation of image using machine learning. Annotations are used 
    as training set for a random forest classifier. Features used for the
    classification are Gabor or Daisy features.

    Parameters
    ----------

    img : ndarray
        Image to be segmented
    mask : ndarray of ints
        Binary array, each connected component corresponds to a different
        object to be segmented
    mode : string
        Type of features used in classification, can be 'daisy' or 'gabor'.
    """
    labels = measure.label(mask)
    if mode == 'daisy':
        if img.ndim > 2:
            img = color.rgb2gray(img)
        radius = 15
        features = feature.daisy(img, step=1, radius=radius, rings=2,
                                 histograms=4)
        crop_labels = labels[radius:-radius, radius:-radius]
    elif mode == 'gabor':
        features = _compute_features_gabor(img)
        nb_ch, nb_fq, sh_1, sh_2 = features.shape
        features = features.reshape((nb_ch * nb_fq, sh_1, sh_2))
        features = np.moveaxis(features, 0, -1)
        crop_labels = labels
    X_train = features[crop_labels > 0, :]
    Y_train = crop_labels[crop_labels > 0]
    rf = RandomForestClassifier(n_estimators=100)
    rf.fit(X_train, Y_train)
    output = rf.predict(features.reshape(-1, features.shape[2]))
    output = output.reshape(crop_labels.shape)
    output[crop_labels > 0] = crop_labels[crop_labels > 0]
    return output


def segmentation_generic(img, mask, mode='watershed'):
    """
    Segmentation of image using annotations as markers (seeds), using
    different algorithms.

    Parameters
    ----------

    img : ndarray
        Image to be segmented
    mask : ndarray of ints
        Binary array, each connected component corresponds to a different
        object to be segmented
    mode : 'string', optional
        Algorithm to be used, can be 'watershed', 'random_walker' or
        'random_forest'.
    """
    if mode=='watershed':
        return watershed_segmentation(img, mask)
    elif mode=='random_walker':
        return random_walker_segmentation(img, mask)
    elif mode=='random_forest':
        return random_forest_segmentation(img, mask)
    else:
        raise NotImplementedError


# -------------------- Background removal ----------------------------


def superpixel_color_segmentation(im, mask, mode='bbox', remove_holes='all',
        object='single', keep_label=True):
    if not mask.shape == im.shape:
        big_mask = np.zeros(im.shape[:2], dtype=np.uint8)
        big_mask[:mask.shape[0], :mask.shape[1]] = mask
        mask = big_mask
    if mask.max() == 1:
        mask = annotation_to_background_mask(mask, mode=mode)
    px = segmentation.felzenszwalb(im, scale=4)
    colors = np.stack([ndimage.mean(im[..., i], px,
                       index=np.arange(0, px.max() + 1))
                       for i in range(3)]).T
    nb_phases = len(np.unique(mask)) - 1
    indices = [np.unique(px[mask == i]) for i in range(1, nb_phases + 1)]
    target = np.concatenate([(i + 1) * np.ones(len(indices[i]))
                             for i in range(len(indices))])
    training_set = np.vstack([colors[indices_i] for indices_i in indices])
    clf = KNeighborsClassifier()
    clf.fit(training_set, target)
    res = clf.predict(colors)
    if keep_label:
        for i in range(1, nb_phases + 1):
            res[indices[i - 1]] = i
    res_img = res[px]
    mask_res = res_img == 1
    if object == 'single':
        labels = measure.label(mask_res)
        index = np.argmax(np.bincount(labels.ravel())[1:]) + 1
        mask_res = labels == index

    if remove_holes == 'all':
        return ndimage.binary_fill_holes(mask_res)
    elif np.isscalar(remove_holes):
        return morphology.remove_small_holes(mask_res, remove_holes)


def annotation_to_background_mask(mask_annotation, mode='bbox'):
    mask = np.copy(mask_annotation).astype(np.uint8)
    if mode == 'bbox':
        sls = ndimage.find_objects(mask_annotation)[0]
        mask[:int(0.9 * sls[0].start)] = 2
        mask[int(1.1 * sls[0].stop):] = 2
        mask[:, :int(0.9 * sls[1].start)] = 2
        mask[:, int(1.1 * sls[1].stop):] = 2
    if mode == 'cvxhull':
        cvx = morphology.convex_hull_image(mask_annotation)
        cvx = morphology.binary_dilation(cvx, np.ones((11, 11)))
        mask[np.logical_not(cvx)] = 2
    return mask
