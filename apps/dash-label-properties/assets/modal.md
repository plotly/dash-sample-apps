## What this app does
This app shows how rich meta-data can be extracted from an image segmentation and linked with 
a datatable to explore the data. 

You see an image of a type of 
[white blood cells](https://en.wikipedia.org/wiki/Agranulocyte). This image is segmented into regions with
[scikit-image](https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_regionprops.html) and we compute
 a number of properties of these regions such as their area, their perimeter length and so on. Each region is then 
 annotated on the image with a color coded contour and is also entered together with its computed properties in 
 the datatable. The annotated image and the datatable are then dynamically linked together, so that interacting with 
 one changes the displayed information in the other. 

## How to use this app
In this app, you can:
- hover over the annotated regions in the image to display information about them and highlight the corresponding
line in the datatable
- control what information is displayed for each region by selecting the corresponding column in the datatable
- choose the region property that the colorscale is computed on
- click on a row in the datatable to draw an outline around the corresponding region
- filter the datatable to constrain the annotated regions
- remove rows from the datatable to delete the annotation of the corresponding regions from the image

See also:
![short video on app function]("assets/label-prop-demo.gif")
## Where to learn more
If you want to learn how to build apps like this, check out:
- [Dash by plotly](https://plotly.com/dash/) for a python based framework to build powerful dashboard apps
- [DataTable filtering syntax](https://dash.plot.ly/datatable/filtering) to learn more on how to filter datatable entries with plotly
- [scikit-image](https://scikit-image.org/docs/stable/user_guide.html) for a scientific image processing python library
- [this app on github](https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-label-properties) to check out the code used to make this web app.
