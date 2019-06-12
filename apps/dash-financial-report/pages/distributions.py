import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go

from utils import Header, make_dash_table
import pandas as pd

df_dividend = pd.read_csv('data/df_dividend.csv')
df_realized = pd.read_csv('data/df_realized.csv')
df_unrealized = pd.read_csv('data/df_unrealized.csv')

layout = html.Div(
    [  # page 5
        html.Div(
            [
                Header(),
                # Row 1
                html.Div(
                    [
                        html.Div(
                            [
                                html.H6(
                                    ["Distributions"],
                                    className="gs-header gs-table-header padded",
                                ),
                                html.Strong(
                                    [
                                        "Distributions for this fund are scheduled quaterly"
                                    ]
                                ),
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
                                html.Br([]),
                                html.H6(
                                    ["Dividend and capital gains distributions"],
                                    className="gs-header gs-table-header tiny-header",
                                ),
                                html.Table(
                                    make_dash_table(df_dividend),
                                    className="tiny-header",
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
                                    ["Realized/unrealized gains as of 01/31/2018"],
                                    className="gs-header gs-table-header tiny-header",
                                )
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
                            [html.Table(make_dash_table(df_realized))],
                            className="six columns",
                        ),
                        html.Div(
                            [html.Table(make_dash_table(df_unrealized))],
                            className="six columns",
                        ),
                    ],
                    className="row ",
                ),
            ],
            className="sub_page",
        )
    ],
    className="page",
)
