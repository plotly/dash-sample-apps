import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go

from utils import Header, make_dash_table
import pandas as pd


df_expenses = pd.read_csv('data/df_expenses.csv')
df_minimums = pd.read_csv('data/df_minimums.csv')


layout = html.Div(
    [  # page 4
        html.Div(
            [
                Header(),
                # Row 1
                html.Div(
                    [
                        html.Div(
                            [
                                html.H6(
                                    ["Expenses"],
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
                                html.Strong(),
                                html.Table(make_dash_table(df_expenses)),
                                html.H6(
                                    ["Minimums"],
                                    className="gs-header gs-table-header padded",
                                ),
                                html.Table(make_dash_table(df_minimums)),
                            ],
                            className="six columns",
                        ),
                        html.Div(
                            [
                                html.Br([]),
                                html.Strong(
                                    "Fees on $10,000 invested over 10 years", style={"color": "#3a3a3a"}),
                                dcc.Graph(
                                    id="graph-6",
                                    figure={
                                        "data": [
                                            go.Bar(
                                                x=["Category Average",
                                                    "This fund"],
                                                y=["2242", "329"],
                                                marker={
                                                    "color": "#97151c"},
                                                name="A",
                                            ),
                                            go.Bar(
                                                x=["This fund"],
                                                y=["1913"],
                                                marker={"color": " #dddddd"},
                                                name="B",
                                            ),
                                        ],
                                        "layout": go.Layout(
                                            annotations=[
                                                {
                                                    "x": -0.0111111111111,
                                                    "y": 2381.92771084,
                                                    "font": {
                                                        "color": "#7a7a7a",
                                                        "family": "Raleway",
                                                        "size": 10,
                                                    },
                                                    "showarrow": False,
                                                    "text": "$2,242",
                                                    "xref": "x",
                                                    "yref": "y",
                                                },
                                                {
                                                    "x": 0.995555555556,
                                                    "y": 509.638554217,
                                                    "font": {
                                                        "color": "#7a7a7a",
                                                        "family": "Raleway",
                                                        "size": 10,
                                                    },
                                                    "showarrow": False,
                                                    "text": "$329",
                                                    "xref": "x",
                                                    "yref": "y",
                                                },
                                                {
                                                    "x": 0.995551020408,
                                                    "y": 1730.32432432,
                                                    "font": {
                                                        "color": "#7a7a7a",
                                                        "family": "Raleway",
                                                        "size": 10,
                                                    },
                                                    "showarrow": False,
                                                    "text": "You save<br><b>$1,913</b>",
                                                    "xref": "x",
                                                    "yref": "y",
                                                },
                                            ],
                                            autosize=False,
                                            height=150,
                                            width=340,
                                            bargap=0.4,
                                            barmode="stack",
                                            hovermode="closest",
                                            margin={"r": 40, "t": 20,
                                                    "b": 20, "l": 40},
                                            showlegend=False,
                                            title="",
                                            xaxis={
                                                "autorange": True,
                                                "range": [-0.5, 1.5],
                                                "showline": True,
                                                "tickfont": {
                                                    "family": "Raleway",
                                                    "size": 10,
                                                },
                                                "title": "",
                                                "type": "category",
                                                "zeroline": False,
                                            },
                                            yaxis={
                                                "autorange": False,
                                                "mirror": False,
                                                "nticks": 3,
                                                "range": [0, 3000],
                                                "showgrid": True,
                                                "showline": True,
                                                "tickfont": {
                                                    "family": "Raleway",
                                                    "size": 10,
                                                },
                                                "tickprefix": "$",
                                                "title": "",
                                                "type": "linear",
                                                "zeroline": False,
                                            },
                                        ),
                                    },
                                    config={"displayModeBar": False},
                                ),
                            ],
                            className="six columns",
                        ),
                    ],
                    className="row ",
                ),
                # Row 3
                html.Div(
                    [
                        html.Div(
                            [
                                html.H6(
                                    ["Fees"],
                                    className="gs-header gs-table-header padded",
                                ),
                                html.Br([]),
                                html.Div(
                                    [
                                        html.Div(
                                            [html.Strong(["Purchase fee"])],
                                            className="three columns right-aligned",
                                        ),
                                        html.Div(
                                            [html.P(["None"])], className="nine columns"
                                        ),
                                    ],
                                    className="row ",
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            [html.Strong(["Redemption fee"])],
                                            className="three columns right-aligned",
                                        ),
                                        html.Div(
                                            [html.P(["None"])], className="nine columns"
                                        ),
                                    ],
                                    className="row ",
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            [html.Strong(["12b-1 fee"])],
                                            className="three columns right-aligned",
                                        ),
                                        html.Div(
                                            [html.P(["None"])], className="nine columns"
                                        ),
                                    ],
                                    className="row ",
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            [html.Strong(
                                                ["Account service fee"])],
                                            className="three columns right-aligned",
                                        ),
                                        html.Div(
                                            [
                                                html.Strong(
                                                    [
                                                        "Nonretirement accounts, traditional IRAs, Roth IRAs, UGMAs/UTMAs, SEP-IRAs, and education savings accounts (ESAs)"
                                                    ]
                                                ),
                                                html.P(
                                                    [
                                                        "We charge a $20 annual account service fee for each Brokerage Account, as well as each individual mutual fund holding with a balance of less than $10,000 in an account. This fee does not apply if you sign up for account and choose electronic delivery of statements, confirmations, and fund reports and prospectuses. This fee also does not apply to members of Flagship Select™, Flagship®, Voyager Select®, and Voyager® Services."
                                                    ]
                                                ),
                                                html.Br([]),
                                                html.Strong(["SIMPLE IRAs"]),
                                                html.P(
                                                    [
                                                        "We charge participants a $25 annual account service fee for each fund they hold in their SIMPLE IRA. This fee does not apply to members of Flagship Select, Flagship, Voyager Select, and Voyager Services."
                                                    ]
                                                ),
                                                html.Br([]),
                                                html.Strong(
                                                    ["403(b)(7) plans"]),
                                                html.P(
                                                    [
                                                        "We charge participants a $15 annual account service fee for each fund they hold in their 403(b)(7) account. This fee does not apply to members of Flagship Select, Flagship, Voyager Select, and Voyager Services."
                                                    ]
                                                ),
                                                html.Br([]),
                                                html.Strong(
                                                    ["Individual 401(k) plans"]
                                                ),
                                                html.P(
                                                    [
                                                        "We charge participants a $20 annual account service fee for each fund they hold in their Individual 401(k) account. This fee will be waived for all participants in the plan if at least 1 participant qualifies for Flagship Select, Flagship, Voyager Select, and Voyager Services"
                                                    ]
                                                ),
                                                html.Br([]),
                                            ],
                                            className="nine columns",
                                        ),
                                    ],
                                    className="row ",
                                ),
                            ],
                            className="twelve columns",
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
