## What this app can do
This dash app allows you to annotate automatically segmented brain regions. 
Here you can see a brain that contains a brain tumour. The dataset is taken from the [MICCAI Brain Tumor Segmentation Challenge](http://braintumorsegmentation.org/). 

You may want to annotate the area of the brain where the tumour is located and then download this annotation for further use.
This app allows you to:
- interactively move through sliced representations of the brain
- see an overlay of the automatically computed over-segmentation of the brain into superpixels (regions of the brain)
- highlight one or several superpixels
- see the highlighted superpixels rendered in a 3D glass brain representation
- download a file with that contains the highlighted superpixels

## How to use this app
When you start the app, you will see two graphs that show axial (Top View) and saggital (Side View) slices through the same brain. 
Overlayed on the brain in orange is an automatically computed segmentation of the brain into areas of similar intensity. 
You can toggle this overlay on and off by clicking on the "Hide/Show Segmentation" button.
These segmented areas are called superpixels, because they extend in all 3 dimensions of the brain volume. 

You can move through the slices of the brain by dragging the sliders under each graph.

To annotate an area of the brain (superpixel) you can mark the area by drawing a short line on it with your mouse cursor. 
All areas that you draw on in this way will then be highlighted in the Top View and Side View representations.

![Image](assets/howto-screenshot.png)

You can undo or redo your last drawing with the corresponding buttons. To see the annotated superpixels in an interactive 3D rendering,
click on the "3D View" button. 
The "Download Selected Partitions" button will download the highlighted areas as a Nifti file and you can download the underlying brain volume with the "Download Brain Volume" button.

You can open, plot, and explore the downloaded .nii files with one of the [plotting functions from the nilearn library](https://nilearn.github.io/plotting/index.html) or with a viewer program like [MRIcroGL](https://www.mccauslandcenter.sc.edu/mricrogl/home).

## Where to find out more about this app
You can find the source code of this app [on Github](https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-3d-image-partitioning).

To find out how to build an app like this yourself, learn more about [Dash](https://plot.ly/dash).

To learn more about how the segmentation is computed, visit [scikit-image](https://scikit-image.org/).
 
