# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.


import base64
import io

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State

import plotly.express as px
import pandas as pd
import skrf as rf

# from skrf.instances import air50 as med
import numpy as np

external_stylesheets = [
    "https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css",
    "https://codepen.io/chriddyp/pen/bWLwgP.css",
]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(
    [
        dcc.Upload(id="upload-ntwk", children=[html.Button("Upload File")]),
        html.Div(id="ntwk-name", children=["Div"]),
    ]
)


@app.callback(
    Output("ntwk-name", "children"),
    Input("upload-ntwk", "contents"),
    State("upload-ntwk", "filename"),
    State("upload-ntwk", "last_modified"),
)
def update_file(contents, filename, dates):
    if contents is not None:
        print(filename)
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string)
        dut = rf.Network()
        str_io = io.StringIO(decoded.decode("utf-8"))
        str_io.name = filename  # due to bug in skrf
        # dut.read_touchstone(str_io)
        dut = rf.Network(str_io)
        # dut = rf.util.get_fid(str_io).name
        return str(dut)
    else:
        return "upload file"


if __name__ == "__main__":
    app.run_server(debug=True)
