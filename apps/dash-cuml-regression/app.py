import os
import time

import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import numpy as np
import cudf
from dash.dependencies import Input, Output, State


# Load data
data_dir = os.environ.get("DATA_DIR", "")
path = os.path.join(data_dir, 'loan_clean.csv')
gdf = cudf.io.csv.read_csv(path)
gdf = gdf.dropna().drop(columns='policy_code')


# Reusable components
def Dropdown(label, **kwargs):
    return dbc.FormGroup(
        [dbc.Label(label), dcc.Dropdown(**kwargs)]
    )


def Slider(label, **kwargs):
    return dbc.FormGroup(
        [dbc.Label(label), dcc.Slider(**kwargs)]
    )

def RadioItems(label, **kwargs):
    return dbc.FormGroup(
        [dbc.Label(label), dcc.RadioItems(**kwargs)]
    )


# Define app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

model_controls = [
    html.B('Model Controls'),
    html.Hr(),
    Dropdown(
        id='model-selected',
        label='Model',
        value='Linear Regression',
        options=[
            {'label': 'Linear Regression', 'value': 'Linear Regression'}
        ]
    )
]

data_controls = [
    html.B('Dataset Controls'),
    html.Hr(),
    Dropdown(
        id='predicted-feature',
        label='Predicted Feature',
        value=gdf.columns[0],
        options=[{'label': x, 'value': x} for x in gdf.columns]
    ),
    RadioItems(
        id='data-size',
        label='Data Size',
        value=500000,
        labelStyle={'padding-right': '10px'},
        options=[
            {'label': '50k', 'value': 50000},
            {'label': '100k', 'value': 100000},
            {'label': '200k', 'value': 200000},
            {'label': '500k', 'value': 500000},
            {'label': '2M', 'value': 2000000},
        ]
    ),
    Slider(
        id='train-split',
        label='Train Split',
        value=0.8, min=0, max=1, step=0.01,
        marks={i/10: f"{i/10:.1f}" for i in range(0, 11, 2)},
    )
]


# Define Layout
app.layout = dbc.Container(fluid=True, children=[
    html.H1('Dash cuML Regression Demo'),
    html.Hr(),
    dbc.Col(dbc.Card(model_controls, body=True), md=4),
    dbc.Col(dbc.Card(data_controls, body=True), md=4),
])


if __name__ == "__main__":
    app.run_server(debug=True)
