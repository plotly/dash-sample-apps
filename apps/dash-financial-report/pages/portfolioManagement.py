import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

from utils import Header, make_dash_table
import pandas as pd
import pathlib

# get relative data folder
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../data").resolve()

df_equity_char = pd.read_csv(DATA_PATH.joinpath("df_equity_char.csv"))
df_equity_diver = pd.read_csv(DATA_PATH.joinpath("df_equity_diver.csv"))


def create_layout(app):
    return html.Div(
        [
            Header(app),
            # page 3
            html.Div(
                [
                    # Row 1
                    html.Div(
                        [
                            html.Div(
                                [html.H6(["Portfolio"], className="subtitle padded")],
                                className="twelve columns",
                            )
                        ],
                        className="rows",
                    ),
                    # Row 2
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.P(["Stock style"], style={"color": "#7a7a7a"}),
                                    dcc.Graph(
                                        id="graph-5",
                                        figure={
                                            "data": [
                                                go.Scatter(
                                                    x=["1"],
                                                    y=["1"],
                                                    hoverinfo="none",
                                                    marker={"opacity": 0},
                                                    mode="markers",
                                                    name="B",
                                                )
                                            ],
                                            "layout": go.Layout(
                                                title="",
                                                annotations=[
                                                    {
                                                        "x": 0.990130093458,
                                                        "y": 1.00181709504,
                                                        "align": "left",
                                                        "font": {
                                                            "family": "Raleway, sans-serif",
                                                            "size": 7,
                                                            "color": "#7a7a7a",
                                                        },
                                                        "showarrow": False,
                                                        "text": "<b>Market<br>Cap</b>",
                                                        "xref": "x",
                                                        "yref": "y",
                                                    },
                                                    {
                                                        "x": 1.00001816013,
                                                        "y": 1.35907755794e-16,
                                                        "font": {
                                                            "family": "Raleway, sans-serif",
                                                            "size": 7,
                                                            "color": "#7a7a7a",
                                                        },
                                                        "showarrow": False,
                                                        "text": "<b>Style</b>",
                                                        "xref": "x",
                                                        "yanchor": "top",
                                                        "yref": "y",
                                                    },
                                                ],
                                                autosize=False,
                                                width=200,
                                                height=150,
                                                hovermode="closest",
                                                margin={
                                                    "r": 30,
                                                    "t": 20,
                                                    "b": 20,
                                                    "l": 30,
                                                },
                                                shapes=[
                                                    {
                                                        "fillcolor": "#f9f9f9",
                                                        "line": {
                                                            "color": "#ffffff",
                                                            "width": 0,
                                                        },
                                                        "type": "rect",
                                                        "x0": 0,
                                                        "x1": 0.33,
                                                        "xref": "paper",
                                                        "y0": 0,
                                                        "y1": 0.33,
                                                        "yref": "paper",
                                                    },
                                                    {
                                                        "fillcolor": "#f2f2f2",
                                                        "line": {
                                                            "color": "#ffffff",
                                                            "dash": "solid",
                                                            "width": 0,
                                                        },
                                                        "type": "rect",
                                                        "x0": 0.33,
                                                        "x1": 0.66,
                                                        "xref": "paper",
                                                        "y0": 0,
                                                        "y1": 0.33,
                                                        "yref": "paper",
                                                    },
                                                    {
                                                        "fillcolor": "#f9f9f9",
                                                        "line": {
                                                            "color": "#ffffff",
                                                            "width": 0,
                                                        },
                                                        "type": "rect",
                                                        "x0": 0.66,
                                                        "x1": 0.99,
                                                        "xref": "paper",
                                                        "y0": 0,
                                                        "y1": 0.33,
                                                        "yref": "paper",
                                                    },
                                                    {
                                                        "fillcolor": "#f2f2f2",
                                                        "line": {
                                                            "color": "#ffffff",
                                                            "width": 0,
                                                        },
                                                        "type": "rect",
                                                        "x0": 0,
                                                        "x1": 0.33,
                                                        "xref": "paper",
                                                        "y0": 0.33,
                                                        "y1": 0.66,
                                                        "yref": "paper",
                                                    },
                                                    {
                                                        "fillcolor": "#f9f9f9",
                                                        "line": {
                                                            "color": "#ffffff",
                                                            "width": 0,
                                                        },
                                                        "type": "rect",
                                                        "x0": 0.33,
                                                        "x1": 0.66,
                                                        "xref": "paper",
                                                        "y0": 0.33,
                                                        "y1": 0.66,
                                                        "yref": "paper",
                                                    },
                                                    {
                                                        "fillcolor": "#f2f2f2",
                                                        "line": {
                                                            "color": "#ffffff",
                                                            "width": 0,
                                                        },
                                                        "type": "rect",
                                                        "x0": 0.66,
                                                        "x1": 0.99,
                                                        "xref": "paper",
                                                        "y0": 0.33,
                                                        "y1": 0.66,
                                                        "yref": "paper",
                                                    },
                                                    {
                                                        "fillcolor": "#f9f9f9",
                                                        "line": {
                                                            "color": "#ffffff",
                                                            "width": 0,
                                                        },
                                                        "type": "rect",
                                                        "x0": 0,
                                                        "x1": 0.33,
                                                        "xref": "paper",
                                                        "y0": 0.66,
                                                        "y1": 0.99,
                                                        "yref": "paper",
                                                    },
                                                    {
                                                        "fillcolor": " #97151c",
                                                        "line": {
                                                            "color": "#ffffff",
                                                            "width": 0,
                                                        },
                                                        "type": "rect",
                                                        "x0": 0.33,
                                                        "x1": 0.66,
                                                        "xref": "paper",
                                                        "y0": 0.66,
                                                        "y1": 0.99,
                                                        "yref": "paper",
                                                    },
                                                    {
                                                        "fillcolor": "#f9f9f9",
                                                        "line": {
                                                            "color": "#ffffff",
                                                            "width": 0,
                                                        },
                                                        "type": "rect",
                                                        "x0": 0.66,
                                                        "x1": 0.99,
                                                        "xref": "paper",
                                                        "y0": 0.66,
                                                        "y1": 0.99,
                                                        "yref": "paper",
                                                    },
                                                ],
                                                xaxis={
                                                    "autorange": True,
                                                    "range": [
                                                        0.989694747864,
                                                        1.00064057995,
                                                    ],
                                                    "showgrid": False,
                                                    "showline": False,
                                                    "showticklabels": False,
                                                    "title": "<br>",
                                                    "type": "linear",
                                                    "zeroline": False,
                                                },
                                                yaxis={
                                                    "autorange": True,
                                                    "range": [
                                                        -0.0358637178721,
                                                        1.06395696354,
                                                    ],
                                                    "showgrid": False,
                                                    "showline": False,
                                                    "showticklabels": False,
                                                    "title": "<br>",
                                                    "type": "linear",
                                                    "zeroline": False,
                                                },
                                            ),
                                        },
                                        config={"displayModeBar": False},
                                    ),
                                ],
                                className="four columns",
                            ),
                            html.Div(
                                [
                                    html.P(
                                        "Calibre Index Fund seeks to track the performance of\
                        a benchmark index that measures the investment return of large-capitalization stocks."
                                    ),
                                    html.P(
                                        "Learn more about this portfolio's investment strategy and policy."
                                    ),
                                ],
                                className="eight columns middle-aligned",
                                style={"color": "#696969"},
                            ),
                        ],
                        className="row ",
                    ),
                    # Row 3
                    html.Br([]),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H6(
                                        ["Equity characteristics as of 01/31/2018"],
                                        className="subtitle padded",
                                    ),
                                    html.Table(
                                        make_dash_table(df_equity_char),
                                        className="tiny-header",
                                    ),
                                ],
                                className=" twelve columns",
                            )
                        ],
                        className="row ",
                    ),
                    # Row 4
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H6(
                                        ["Equity sector diversification"],
                                        className="subtitle padded",
                                    ),
                                    html.Table(
                                        make_dash_table(df_equity_diver),
                                        className="tiny-header",
                                    ),
                                ],
                                className=" twelve columns",
                            )
                        ],
                        className="row ",
                    ),
                ],
                className="sub_page",
            ),
        ],
        className="page",
    )
