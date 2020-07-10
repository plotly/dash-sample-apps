The code of this app is [available here](https://github.com/plotly/dash-sample-apps/blob/master/apps/dash-image-segmentation/app.py). 

This app uses machine learning in order to compute the segmentation of an
image, given user-provided annotations. In addition to plotly and Dash, the app
 uses off-the-shelf algorithms and estimators from PyData packages, namely
[scikit-image](https://scikit-image.org/) and [scikit-learn](https://scikit-learn.org).

### User annotations

To build the training set, we use the [new shape drawing capabilities
of plotly.py](https://eoss-image-processing.github.io/2020/05/06/shape-drawing.html)
and in particular the `drawopenpath` dragmode which can used to draw
"squiggles" on parts of the image which you want to label. You can adjust the width of the squiggle with a Dash
`dcc.Slider` to make it possible to annotate features of different sizes. Each time a new annotation is drawn, it is captured by the plotly figure's
[relayoutData event, which triggers a callback](https://dash.plotly.com/interactive-graphing). You should select a different color for every type of object you wish to segment.

### Computation of features

In machine learning, a sample is represented as a vector of *features*. Deep
learning models learn features directly from the data and are very popular for
image processing. Nevertheless, they require a large training set (a large
number of images) and their training is very resource-intensive. We instead use local features which,
for each pixel, represent
- the average intensity in a small region around the pixel
- the average magnitude of gradients in the same region
- measures of local texture in this region

Such features are computed by first convolving the image of interest with a Gaussian
kernel, and then measuring the local color intensity, gradient intensity, or the
eigenvalues of the Hessian matrix. Conveniently, these operations are provided
by the [`filters` module of `scikit-image`](https://scikit-image.org/docs/stable/api/skimage.filters.html)
and are relatively fast, since they operate on local neighbourhoods.

## Model training and prediction

Features are extracted for the annotated pixels, and passed to a scikit-learn [Random Forest Classifier](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html). This estimator belongs to the class of [ensemble methods](https://scikit-learn.org/stable/modules/ensemble.html), where the predictions by several base estimators are combined to improve the generalizability or robustness of the prediction. After the model is trained, its prediction is computed on unlabeled pixels, resulting in a segmentation of the image. It is possible to add more annotations to improve the segmentation if some pixels are wrongly classified. Furthermore, one can download the estimator in order to classify new images of the same type (for example, a time series).


