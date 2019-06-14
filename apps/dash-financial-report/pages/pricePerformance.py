import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go
from utils import Header, make_dash_table
import pandas as pd

df_current_prices = pd.read_csv("data/df_current_prices.csv")
df_hist_prices = pd.read_csv("data/df_hist_prices.csv")
df_avg_returns = pd.read_csv("data/df_avg_returns.csv")
df_after_tax = pd.read_csv("data/df_after_tax.csv")
df_recent_returns = pd.read_csv("data/df_recent_returns.csv")

df_graph = pd.read_csv("data/df_graph.csv")

layout = html.Div(
    [  # page 2
        html.Div(
            [
                Header(),
                # Row
                html.Div(
                    [
                        html.Div(
                            [
                                html.H6(
                                    ["Current Prices"], className="subtitle padded"
                                ),
                                html.Table(make_dash_table(df_current_prices)),
                            ],
                            className="six columns",
                        ),
                        html.Div(
                            [
                                html.H6(
                                    ["Historical Prices"], className="subtitle padded"
                                ),
                                html.Table(make_dash_table(df_hist_prices)),
                            ],
                            className="six columns",
                        ),
                    ],
                    className="row ",
                ),
                # Row 2
                html.Div(
                    [
                        html.Div(
                            [
                                html.H6("Performance", className="subtitle padded"),
                                dcc.Graph(
                                    id="graph-4",
                                    figure={
                                        "data": [
                                            go.Scatter(
                                                x=df_graph["Date"],
                                                y=df_graph["Calibre Index Fund"],
                                                line={"color": "#97151c"},
                                                mode="lines",
                                                name="Calibre Index Fund",
                                            ),
                                            go.Scatter(
                                                x=df_graph["Date"],
                                                y=df_graph[
                                                    "MSCI EAFE Index Fund (ETF)"
                                                ],
                                                line={"color": "#b5b5b5"},
                                                mode="lines",
                                                name="MSCI EAFE Index Fund (ETF)",
                                            ),
                                        ],
                                        "layout": go.Layout(
                                            autosize=False,
                                            width=700,
                                            height=200,
                                            font={"family": "Raleway", "size": 10},
                                            margin={"r": 40, "t": 40, "b": 30, "l": 40},
                                            showlegend=True,
                                            titlefont={"family": "Raleway", "size": 10},
                                            xaxis={
                                                "autorange": True,
                                                "range": ["2007-12-31", "2018-03-06"],
                                                "rangeselector": {
                                                    "buttons": [
                                                        {
                                                            "count": 1,
                                                            "label": "1Y",
                                                            "step": "year",
                                                            "stepmode": "backward",
                                                        },
                                                        {
                                                            "count": 3,
                                                            "label": "3Y",
                                                            "step": "year",
                                                            "stepmode": "backward",
                                                        },
                                                        {
                                                            "count": 5,
                                                            "label": "5Y",
                                                            "step": "year",
                                                        },
                                                        {
                                                            "count": 10,
                                                            "label": "10Y",
                                                            "step": "year",
                                                            "stepmode": "backward",
                                                        },
                                                        {"label": "All", "step": "all"},
                                                    ]
                                                },
                                                "showline": True,
                                                "type": "date",
                                                "zeroline": False,
                                            },
                                            yaxis={
                                                "autorange": True,
                                                "range": [18.6880162434, 278.431996757],
                                                "showline": True,
                                                "type": "linear",
                                                "zeroline": False,
                                            },
                                        ),
                                    },
                                    config={"displayModeBar": False},
                                ),
                            ],
                            className="twelve columns",
                        )
                    ],
                    className="row ",
                ),
                # Row 3
                html.Div(
                    [
                        html.Div(
                            [
                                html.H6(
                                    [
                                        "Average annual returns--updated monthly as of 02/28/2018"
                                    ],
                                    className="subtitle padded",
                                ),
                                html.Table(
                                    make_dash_table(df_avg_returns),
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
                                    [
                                        "After-tax returns--updated quarterly as of 12/31/2017"
                                    ],
                                    className="subtitle padded",
                                ),
                                html.Table(
                                    make_dash_table(df_after_tax),
                                    className="tiny-header",
                                ),
                            ],
                            className=" twelve columns",
                        )
                    ],
                    className="row ",
                ),
                # Row 5
                html.Div(
                    [
                        html.Div(
                            [
                                html.H6(
                                    ["Recent investment returns"],
                                    className="subtitle padded",
                                ),
                                html.Table(
                                    make_dash_table(df_recent_returns),
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
