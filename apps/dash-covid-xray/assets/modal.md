## What this app can do

This app displays 3-D chest tomography data of a patient with Covid-19. The imaging data show 
[Ground Glass Opacities (GGO) in the lung](https://en.wikipedia.org/wiki/Ground-glass_opacity#COVID-19), a common imaging 
finding among Covid-19 patients. The purpose of this app is to extract the spatial extent of the GGOs in the imaging
data for the 
[eventual purpose of training machine learning models](https://eoss-image-processing.github.io/2020/12/16/ct-app.html) 
to do this automatically.
The imaging data in this app are displayed with [`dash-slicer`](https://dash.plotly.com/slicer), a Dash component that 
provides interactive slicing views of volumetric data. The data used in this app come from the open dataset of
the [COVID-19 image data collection](https://github.com/ieee8023/covid-chestxray-dataset).

## How to use this app
With this app you can:

- Interactively explore the 3-D chest tomography data with two linked `dash-slicer` viewers. Note that the position of 
  the left (axial) slice of the lung is displayed with a blue line in the right (sagittal) slice of the lung. Similarly, 
  the position of the right slice is displayed with an orange line in the left viewer.
- Draw an outline of the GGOs on the axial viewer. Make sure that this outline encompasses all GGO across all of the 
  axial slices.
- Indicate the vertical span of GGOs in the image by drawing a rectangle on the saggital viewer.
- Selecting a range of image intensity values that reflect GGOs in the outlined region of interest by drawing a 
  selection in the histogram. Note that you first need to define the region of interest by drawing the outline and 
  height of GGOs in the axial and sagittal viewers.
  
## Where to learn more
If you want to learn how to build apps like this, check out:
- [Dash by plotly](https://plotly.com/dash/) for a python based framework to build powerful dashboard apps
- [Dash Slicer](https://dash.plotly.com/slicer) to learn more on how to make interactive 3D image slicers in Dash
- [scikit-image](https://scikit-image.org/docs/stable/user_guide.html) for a scientific image processing python library
- [this app on github](https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-covid-xray) to check out the code used to make this web app.
- [our blog post on this app](https://eoss-image-processing.github.io/2020/12/16/ct-app.html)



