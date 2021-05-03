# Because I have to learn Dash

# run with python chess_app.py and visit
# http://127.0.0.1:8050/ in your web browser.


# Imports
import dash
import ast

import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc


import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
from chessboard import *
from styles import *

# Read the .csv file with the preprocessed data.
df = pd.read_csv(
    "chess_app.csv",
    dtype={"pawns": int, "knights": int, "bishops": int, "rooks": int, "queens": int},
    converters={"wKing_sqr": ast.literal_eval, "bKing_sqr": ast.literal_eval},
)


# Optionally, produces a .csv of such a dataframe.
# board_output(wKing_sqr).to_csv("wKing_Heatmap.csv")
df = board_output(df["wKing_sqr"])

# FILLER STUFF ~ LET'S KEEP THIS FILE CLEAN, have other .py files with everything!
x_coords = ["A", "B", "C", "D", "E", "F", "G", "H"]
replacer = {i + 1: x for i, x in enumerate(x_coords)}
df = (
    df.stack()
    .reset_index()
    .rename(columns={"level_0": "rows", "level_1": "cols", 0: "freq"})
)
df.iloc[:, 0:2] = df.iloc[:, 0:2].apply(lambda x: x + 1)
df["letters"] = df.cols.replace(replacer)

# Set stylesheets and app.
external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "Chess Analytics"

chessboard = getChessboard()

chessboard.add_trace(getHeatmap(dataframe=df))

# Defining app layout
app.layout = html.Div(
    style={"backgroundColor": colors["background"]},
    children=[
        html.Div(
            [html.H2("Chess App"), html.Img(src="/assets/chess-app.jpg"),],
            className="banner",
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Label("Range Slider"),
                        dcc.RangeSlider(
                            id="rangeslider",
                            marks={i: "label {}".format(i) for i in range(0, 5)},
                            min=0,
                            max=5,
                            value=[0, 1],
                            step=1,
                        ),
                    ],
                    className="five columns",
                ),
                html.Div(
                    [
                        dcc.Graph(
                            id="id_chessboard",
                            figure=chessboard,
                            config={
                                "displayModeBar": False,
                                "scrollZoom": False,
                                "showAxisDragHandles": False,
                            },
                        )
                    ],
                    className=" six columns",
                ),
            ],
            className="row",
        ),
    ],
)
# Statring the dash app
if __name__ == "__main__":
    app.run_server(debug=True)
"""
app.layout = html.Div([
    html.Label('Range Slider'),
    dcc.RangeSlider(
        id="rangeslider",
        marks={i: 'label {}'.format(i) for i in range(0, 5)},
        min=0,
        max=5,
        value=[0, 1],
        step=1
    ),
])
"""
"""
fig = px.imshow(df_test, color_continuous_scale=px.colors.sequential.Burgyl, labels=dict(
    x="", y="", color="No of checkmates"),
    x=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
    y=['1', '2', '3', '4', '5', '6', '7', '8'])
"""
