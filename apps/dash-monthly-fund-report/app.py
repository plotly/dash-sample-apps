import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pathlib

from plotly import graph_objs as go
from datetime import datetime as dt
import json
import pandas as pd
import os


app = dash.Dash(__name__)
server = app.server


df_fund_data = pd.read_csv("https://plot.ly/~jackp/17534.csv")
df_fund_data.head()
df_perf_summary = pd.read_csv("https://plot.ly/~jackp/17530.csv")
df_perf_summary.head()
df_cal_year = pd.read_csv("https://plot.ly/~jackp/17528.csv")
df_cal_year.head()
df_perf_pc = pd.read_csv("https://plot.ly/~jackp/17532.csv")


def make_dash_table(df):
    """ Return a dash definition of an HTML table for a Pandas dataframe """
    table = []
    for index, row in df.iterrows():
        html_row = []
        for i in range(len(row)):
            html_row.append(html.Td([row[i]]))
        table.append(html.Tr(html_row))
    return table


modifed_perf_table = make_dash_table(df_perf_summary)


modifed_perf_table.insert(
    0,
    html.Tr(
        [
            html.Td([]),
            html.Td(["Cumulative"], colSpan=4, style={"text-align": "center"}),
            html.Td(["Annualised"], colSpan=4, style={"text-align": "center"}),
        ],
        style={"background": "white", "font-weight": "600"},
    ),
)

df_fund_info = pd.read_csv("https://plot.ly/~jackp/17544.csv")
df_fund_characteristics = pd.read_csv("https://plot.ly/~jackp/17542.csv")
df_fund_facts = pd.read_csv("https://plot.ly/~jackp/17540.csv")
df_bond_allocation = pd.read_csv("https://plot.ly/~jackp/17538.csv")


PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()

df_sector_allocation = pd.read_csv(DATA_PATH.joinpath("sector-allocation.csv"))
df_performance = pd.read_csv(DATA_PATH.joinpath("performance.csv"))
df_currency_weights = pd.read_csv(DATA_PATH.joinpath("currency-weights.csv"))
df_credit_allocation = pd.read_csv(DATA_PATH.joinpath("credit-allocation.csv"))


def header():
    return html.Div(
        [
            html.Div(
                className="row",
                children=html.Img(className="logo", src=app.get_asset_url("logo.png")),
            ),
            html.Div(
                className="row",
                children=[
                    html.Div(
                        className="header-left",
                        children=[
                            html.H1("Dash Monthly Fund Report"),
                            html.H2("A sub-fund of Dash Monthly Fund, SICAV"),
                        ],
                    ),
                    html.Div(
                        className="header-right",
                        children=[
                            html.H1(
                                children=[
                                    html.Span("03", className="light-blue"),
                                    html.Span("17", className="blue"),
                                ]
                            ),
                            html.H6("Monthly Fund Update"),
                        ],
                    ),
                ],
            ),
        ]
    )


# Describe the layout, or the UI, of the app
app.layout = html.Div(
    [
        # Page 1
        html.Div(
            [
                # Subpage 1
                html.Div(
                    [
                        # Row 1
                        header(),
                        html.Br([]),
                        # Row 2
                        html.Div(
                            className="spec-row",
                            children=[
                                html.Div(
                                    className="six columns div-investor-profile",
                                    children=[
                                        html.H5("Investor Profile"),
                                        html.H6("Investor objective"),
                                        html.P("Capital appreciation and income."),
                                        html.H6(
                                            "Position in your overall investment portfolio"
                                        ),
                                        html.P(
                                            "The fund can complement your portfolio."
                                        ),
                                    ],
                                ),
                                html.Div(
                                    className="six columns div-fund-designed-for",
                                    children=[
                                        html.H4("The fund is designed for:"),
                                        html.P(
                                            "The fund is designed for investors who are looking for a flexible \
                                global investment and sub-investment grade fixed income portfolio \
                                that has the ability to alter its exposure with an emphasis on interest \
                                rates, currencies and credit markets and that seeks to generate returns \
                                through different market conditions with a riskier investment strategy \
                                than GS Strategic Absolute Return Bond I Portfolio."
                                        ),
                                    ],
                                ),
                            ],
                        ),
                        # Row 3
                        html.Div(
                            className="row",
                            children=[
                                # Left Side
                                html.Div(
                                    className="six columns",
                                    children=[
                                        # Performace %
                                        html.Div(
                                            children=[
                                                html.H6("Performance (%)"),
                                                html.Table(make_dash_table(df_perf_pc)),
                                            ]
                                        ),
                                        # Fund Data
                                        html.Div(
                                            children=[
                                                html.H6("Fund Data"),
                                                html.Table(
                                                    make_dash_table(df_fund_data)
                                                ),
                                            ]
                                        ),
                                    ],
                                ),
                                # Right Side
                                html.Div(
                                    className="six columns div-graphs",
                                    children=[
                                        html.Div(
                                            [
                                                html.H6("Performance (Indexed)"),
                                                dcc.Graph(
                                                    figure={
                                                        "data": [
                                                            {
                                                                "uid": "4cd1a4",
                                                                "line": {
                                                                    "color": "#119df",
                                                                    "width": 3,
                                                                },
                                                                "mode": "lines",
                                                                "name": "Absolute Return Bond II Portfolio Base Shares",
                                                                "type": "scatter",
                                                                "x": df_performance[
                                                                    "x"
                                                                ],
                                                                "y": df_performance[
                                                                    "y1"
                                                                ],
                                                            },
                                                            {
                                                                "uid": "f7fed3",
                                                                "line": {
                                                                    "dash": "dash",
                                                                    "color": "#2a3f5f",
                                                                    "width": 3,
                                                                },
                                                                "mode": "lines",
                                                                "name": "3 Month Libor (USD)",
                                                                "type": "scatter",
                                                                "x": df_performance[
                                                                    "x"
                                                                ],
                                                                "y": df_performance[
                                                                    "y2"
                                                                ],
                                                                "connectgaps": True,
                                                            },
                                                        ],
                                                        "layout": {
                                                            "xaxis": {
                                                                "type": "date",
                                                                "ticks": "outside",
                                                                "title": "",
                                                                "showline": True,
                                                                "showgrid": True,
                                                                "tickfont": {
                                                                    "color": "rgb(68, 68, 68)"
                                                                },
                                                                "gridcolor": "#BDC1C470",
                                                                "tickformat": "%b %Y",
                                                            },
                                                            "yaxis": {
                                                                "type": "linear",
                                                                "ticks": "outside",
                                                                "range": [80, 135],
                                                                "nticks": 11,
                                                                "showline": True,
                                                                "showgrid": True,
                                                                "gridcolor": "#BDC1C470",
                                                                "fixedrange": True,
                                                            },
                                                            "legend": {
                                                                "x": 0.6,
                                                                "y": 0,
                                                                "bgcolor": "#ecf7fd70",
                                                                "font": dict(size=7.5),
                                                            },
                                                            "margin": dict(
                                                                b=40, l=35, r=0, t=10
                                                            ),
                                                            "height": 200,
                                                            "hovermode": "closest",
                                                            "showlegend": True,
                                                        },
                                                    }
                                                ),
                                                html.P(
                                                    "This is an actively managed fund that is not designed to track its reference benchmark. \
                                                    Therefore the performance of the fund and the performance of its reference benchmark \
                                                    may diverge. In addition stated reference benchmark returns do not reflect any management \
                                                    or other charges to the fund, whereas stated returns of the fund do."
                                                ),
                                                html.P(
                                                    "Past performance does not guarantee future results, which may vary. \
                                                    The value of investments and the income derived from investments will fluctuate and \
                                                    can go down as well as up. A loss of capital may occur."
                                                ),
                                            ]
                                        )
                                    ],
                                ),
                            ],
                        ),
                        # Row 3
                        html.Div(
                            className="row",
                            children=[
                                html.Div(
                                    children=[
                                        html.H6("Performance Summary (%)"),
                                        html.Table(
                                            modifed_perf_table, className="reversed"
                                        ),
                                        html.H6("Calendar Year Performance (%)"),
                                        html.Table(make_dash_table(df_cal_year)),
                                    ]
                                )
                            ],
                        ),
                    ],
                    className="subpage",
                )
            ],
            className="page",
        ),
        # Page 2
        html.Div(
            [
                # Subpage 2
                html.Div(
                    children=[
                        # Row 1
                        header(),
                        html.Br([]),
                        # Row 2
                        html.Div(
                            className="row",
                            children=[
                                # Left Column
                                html.Div(
                                    className="five columns",
                                    children=[
                                        html.H6("Financial Information"),
                                        html.Table(make_dash_table(df_fund_info)),
                                        html.H6("Fund Characteristics"),
                                        html.Table(
                                            make_dash_table(df_fund_characteristics)
                                        ),
                                        html.H6("Fund Facts"),
                                        html.Table(make_dash_table(df_fund_facts)),
                                        html.H6("Country Bond Allocation (%)"),
                                        html.Table(make_dash_table(df_bond_allocation)),
                                    ],
                                ),
                                # Right Column
                                html.Div(
                                    className="seven columns div-graphs",
                                    children=[
                                        html.H6("Sector Allocation (%)"),
                                        dcc.Graph(
                                            figure={
                                                "data": [
                                                    {
                                                        "uid": "353874",
                                                        "mode": "markers",
                                                        "name": "B",
                                                        "type": "bar",
                                                        "x": df_sector_allocation["x"],
                                                        "y": df_sector_allocation["y"],
                                                        "marker": {"color": "#119dff"},
                                                    }
                                                ],
                                                "layout": {
                                                    "width": 380,
                                                    "height": 220,
                                                    "font": dict(size=9),
                                                    "xaxis": {
                                                        "type": "category",
                                                        "range": [-1, 12],
                                                    },
                                                    "yaxis": {"type": "linear"},
                                                    "margin": dict(t=0, r=0, b=90, l=0),
                                                    "hovermode": "closest",
                                                    "bargroupgap": 0.2,
                                                },
                                            }
                                        ),
                                        html.H6("Top 10 Currency Weights (%)"),
                                        dcc.Graph(
                                            figure={
                                                "data": [
                                                    {
                                                        "uid": "80eb70",
                                                        "name": "Col1",
                                                        "type": "bar",
                                                        "x": df_currency_weights["x"],
                                                        "y": df_currency_weights["y"],
                                                        "marker": {"color": "#119dff"},
                                                        "orientation": "h",
                                                    }
                                                ],
                                                "layout": {
                                                    "title": "",
                                                    "width": 380,
                                                    "font": dict(size=9),
                                                    "xaxis": {
                                                        "type": "linear",
                                                        "range": [
                                                            -62.13380540229584,
                                                            164.54230264362104,
                                                        ],
                                                        "ticks": "outside",
                                                        "nticks": 6,
                                                        "showline": True,
                                                        "zeroline": False,
                                                        "showgrid": False,
                                                        "autorange": True,
                                                        "ticksuffix": "%",
                                                    },
                                                    "yaxis": {
                                                        "type": "category",
                                                        "range": [-0.5, 11.5],
                                                        "showgrid": True,
                                                        "autorange": True,
                                                    },
                                                    "height": 250,
                                                    "margin": dict(
                                                        b=40, l=100, r=5, t=40
                                                    ),
                                                },
                                            }
                                        ),
                                        html.H6("Credit Allocation (%)"),
                                        dcc.Graph(
                                            figure={
                                                "data": [
                                                    {
                                                        "uid": "a8c61b",
                                                        "name": "GS Strategic<br>Absolute<br>Return<br>Bond II<br>Portfolio",
                                                        "type": "bar",
                                                        "xsrc": "jackp:17802:62c223",
                                                        "x": df_credit_allocation["x1"],
                                                        "ysrc": "jackp:17802:d66b98",
                                                        "y": df_credit_allocation["y"],
                                                        "marker": {"color": "#c2ebff"},
                                                        "visible": True,
                                                        "orientation": "h",
                                                    },
                                                    {
                                                        "uid": "f84602",
                                                        "name": "3 Month<br>USD Libor",
                                                        "type": "bar",
                                                        "xsrc": "jackp:17802:4a6d89",
                                                        "x": df_credit_allocation["x2"],
                                                        "ysrc": "jackp:17802:d66b98",
                                                        "y": df_credit_allocation["y"],
                                                        "marker": {"color": "#119dff"},
                                                        "orientation": "h",
                                                    },
                                                ],
                                                "layout": {
                                                    "title": "",
                                                    "width": 380,
                                                    "font": dict(size=9),
                                                    "xaxis": {
                                                        "type": "linear",
                                                        "range": [0, 100],
                                                        "ticks": "outside",
                                                        "nticks": 11,
                                                        "showline": True,
                                                        "showgrid": False,
                                                        "ticksuffix": "%",
                                                    },
                                                    "yaxis": {
                                                        "type": "category",
                                                        "showgrid": True,
                                                        "autorange": True,
                                                    },
                                                    "height": 195,
                                                    "legend": {
                                                        "x": 0.6,
                                                        "y": 0.75,
                                                        "bgcolor": "#ecf7fd70",
                                                    },
                                                    "margin": dict(
                                                        b=40, l=60, r=0, t=10, pad=0
                                                    ),
                                                },
                                            }
                                        ),
                                    ],
                                ),
                            ],
                        ),
                    ],
                    className="subpage",
                )
            ],
            className="page",
        ),
    ]
)

if __name__ == "__main__":
    app.server.run()
