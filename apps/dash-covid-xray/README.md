# Dash Exploration of COVID-19 chest X-ray CT images

## About this app

This app shows how to explore 3-D chest tomography data using Dash. 

The data used in this app come from the open dataset of
https://github.com/ieee8023/covid-chestxray-dataset

## How to run this app

(The following instructions apply to Windows command line.)

To run this app first clone repository and then open a terminal to the app folder.

```
git clone https://github.com/plotly/dash-sample-apps.git
cd dash-sample-apps/apps/dash-covid-xray
```

Create and activate a new virtual environment (recommended) by running
the following:

On Windows

```
virtualenv venv 
\venv\scripts\activate
```

Or if using linux

```bash
python3 -m venv myvenv
source myvenv/bin/activate
```

Install the requirements:

```
pip install -r requirements.txt
```
Run the app:

```
python app.py
```
You can run the app on your browser at http://127.0.0.1:8050


## Resources

To learn more about Dash, please visit [documentation](https://plot.ly/dash).

