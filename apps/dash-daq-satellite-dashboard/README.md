# Dash-DAQ-Satellite-Dashboard

## Introduction
A Dash application that tracks satellites and displays live data captured by them.

### Satellite
Artificial satellites are objects placed into orbit for various tasks, such as surveillance and transferring radio data 
across the world. It's important to monitor satellites to ensure that they can accomplish their jobs, so information such as
their position and elevation are, for example, useful to determine whether or not a satellite is deviating from its original path
due to unforeseen circumstances.

### dash-daq
[Dash-DAQ](http://dash-daq.netlify.com/#about) is a data acquisition and control package built on top of Plotly's 
[Dash](https://plot.ly/products/dash/). It comprises a robust set of controls that make it simpler to integrate data 
acquisition and controls into your Dash applications.


## Requirements
We suggest you to create a virtual environment for python3 to run this app. To do so, run:
```bash
python3 python3 -m virtualenv [your environment name]
```
```bash
source activate [your environment name]
```
To install all of this app-specific required packages to this environment, simply run:
```bash
pip install -r requirements.txt
```


## How to use the app
To run this application locally, simply type the following into the command line:
```bash
python app.py
```
--

![Satellite Dashboard](assets/screenshot.png)

### Controls
* Satellite dropdown: Select which satellite to track.
* Histogram: Data is updated every 2 seconds, and to view the histogram for a desired data type, simply click on the
corresponding Dash component.
* Path toggle: Show and hide the expected satellite path.
* Time toggle: Display data from the past hour or the past minute. 


### Resources
#### Dash
Dash abstracts away all of the technologies and protocols required to build an interactive web-based application, and 
is a simple and effective way to bind a user interface around your Python code. To learn more about Dash, check out our 
[documentation](https://dash.plot.ly/). 
