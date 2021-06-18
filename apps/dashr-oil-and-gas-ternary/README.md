## dashr-oil-gas-ternary

`dashr-oil-gas-ternary` creates a dashboard for mineral composition evaluations from natural gas wells.  
By pairing a geographic map of well locations with ternary diagram, this app helps geologists to find out results of gas wells with relationship to formations or rock types.
This is a demo of Dash interactive Python framework developed by [Plotly](https//plot.ly/).

![Animated](../dash-oil-gas-ternary/assets/Screencast.gif)

## Screenshots
![initial](assets/Screenshot.png)

## Requirements
We suggest you to create a separate virtual environment running Python 3 for this app, and install all of the required dependencies there. Run in Terminal/Command Prompt:

```
git clone https://github.com/plotly/dash-sample-apps.git
Rscript app.R
```

## How to use the app
Run this app locally by:
```
RScript app.R
```
Open http://127.0.0.1:8050/ in your browser.

Select a subset of data points from the well map, ternary map or bar graph to visualize cross-filtering to other plots.
Selection can be done by clicking on individual data points, or by using the lasso tool to capture multiple data points or bars. Hold the SHIFT
key while clicking and dragging for multi-region selection.

## Data source
Historical shale gas production data for selected Devonian Ohio Shale wells http://www.uky.edu/KGS/emsweb/kyogfaq/kyogfaq10.html

