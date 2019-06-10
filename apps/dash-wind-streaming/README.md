# Dash Wind Streaming

## About this app

This app queries a SQL database every second and uses the data to update the wind speed diagram and the wind direction diagram. 
The wind speed values are then binned in real time to generate the wind histogram plot.

Original repo: [plotly/dash-wind-streaming](https://github.com/plotly/dash-wind-streaming)


## How to run this app

(The following instructions apply to Posix/bash. Windows users should check
[here](https://docs.python.org/3/library/venv.html).)

First, clone this repository and open a terminal inside the root folder.

Create and activate a new virtual environment (recommended) by running
the following:

```bash
python3 -m venv myvenv
source myvenv/bin/activate
```

Install the requirements:

```bash
pip install -r requirements.txt
```
Run the app:

```bash
python app.py
```
Open a browser at http://127.0.0.1:8050

## Screenshots

![demo.gif](demo.gif)

## Resources

- To learn more about Dash, check out our [documentation](https://plot.ly/dash).