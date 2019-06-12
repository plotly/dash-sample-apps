import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go

from utils import Header, make_dash_table
import pandas as pd

df_equity_char = pd.read_csv('data/df_equity_char.csv')
df_equity_diver = pd.read_csv('data/df_equity_diver.csv')

layout = html.Div(
    [  # page 3
        html.Div(
            [
                Header(),
                # Row 1
                html.Div(
                    [
                        html.Div(
                            [
                                html.H6(
                                    ["Portfolio"],
                                    className="gs-header gs-table-header padded",
                                )
                            ],
                            className="twelve columns",
                        )
                    ],
                    className="row ",
                ),
                # Row 2
                html.Div(
                    [
                        html.Div(
                            [
                                html.Strong(["Stock style"]),
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
                                                        "family": "Raleway",
                                                        "size": 9,
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
                                                        "family": "Raleway",
                                                        "size": 9,
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
                                            margin={"r": 30, "t": 20,
                                                    "b": 20, "l": 30},
                                            shapes=[
                                                {
                                                    "fillcolor": "#f9f9f9",
                                                    "line": {
                                                        "color": "#ffffff",
                                                        "width": 0,
                                                    },
                                                    # "opacity": 1,
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
                                                    # "opacity": 0.3,
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
                                                    # "opacity": 0.3,
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
                                                    # "opacity": 0.3,
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
                                                    # "opacity": 0.3,
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
                                                    # "opacity": 0.3,
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
                                                    # "opacity": 0.3,
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
                                                    # "opacity": 0.9,
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
                                                    # "opacity": 0.3,
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
                                    className="gs-header gs-table-header tiny-header",
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
                                    className="gs-header gs-table-header tiny-header",
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
        )
    ],
    className="page",
)
