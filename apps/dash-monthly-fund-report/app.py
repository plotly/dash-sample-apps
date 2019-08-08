import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html

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

# Describe the layout, or the UI, of the app
app.layout = html.Div(
    [
        # Page 1
        html.Div(
            [
                html.A(["Print PDF"], className="button no-print"),
                # Subpage 1
                html.Div(
                    [
                        # Header
                        html.Div(
                            className="row",
                            children=html.Img(
                                className="logo", src=app.get_asset_url("logo.png")
                            ),
                        ),
                        html.Div(
                            className="row",
                            children=[
                                html.Div(
                                    className="nine columns header-left",
                                    children=[
                                        html.H1(
                                            "Goldman Sachs Strategic Absolute Return Bond II Portfolio"
                                        ),
                                        html.H2(
                                            "A sub-fund of Goldman Sachs Funds, SICAV"
                                        ),
                                    ],
                                ),
                                html.Div(
                                    className="three columns header-right",
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
                                html.Div(
                                    className="six columns",
                                    children=[
                                        html.H6("Fund Data"),
                                        html.Table(make_dash_table(df_fund_data)),
                                    ],
                                )
                            ],
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H6(
                                            ["Performance (Indexed)"],
                                            className="gs-header gs-table-header padded",
                                        ),
                                        html.Iframe(
                                            src="https://plot.ly/~jackp/17553.embed?modebar=false&link=false&autosize=true",
                                            style={
                                                "border": "0",
                                                "width": "100%",
                                                "height": "250",
                                            },
                                        ),
                                    ],
                                    className="eight columns",
                                )
                            ],
                            className="row ",
                        ),
                        # Row 2.5
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H6(
                                            "Performance (%)",
                                            className="gs-header gs-text-header padded",
                                        ),
                                        html.Table(
                                            make_dash_table(df_perf_pc),
                                            className="tiny-header",
                                        ),
                                    ],
                                    className="four columns",
                                ),
                                html.Div(
                                    [
                                        html.P(
                                            "This is an actively managed fund that is not designed to track its reference benchmark. \
                        Therefore the performance of the fund and the performance of its reference benchmark \
                        may diverge. In addition stated reference benchmark returns do not reflect any management \
                        or other charges to the fund, whereas stated returns of the fund do."
                                        ),
                                        html.Strong(
                                            "Past performance does not guarantee future results, which may vary. \
                        The value of investments and the income derived from investments will fluctuate and \
                        can go down as well as up. A loss of capital may occur."
                                        ),
                                    ],
                                    className="eight columns",
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
                                            "Fund Data",
                                            className="gs-header gs-text-header padded",
                                        ),
                                        html.Table(make_dash_table(df_fund_data)),
                                    ],
                                    className="four columns",
                                ),
                                html.Div(
                                    [
                                        html.H6(
                                            "Performance Summary (%)",
                                            className="gs-header gs-table-header padded",
                                        ),
                                        html.Table(
                                            modifed_perf_table, className="reversed"
                                        ),
                                        html.H6(
                                            "Calendar Year Performance (%)",
                                            className="gs-header gs-table-header padded",
                                        ),
                                        html.Table(make_dash_table(df_cal_year)),
                                    ],
                                    className="eight columns",
                                ),
                            ],
                            className="row ",
                        ),
                    ],
                    className="subpage",
                ),
            ],
            className="page",
        ),
        html.Div(
            [  # page 2
                html.A(
                    ["Print PDF"],
                    className="button no-print",
                    style={"position": "absolute", "top": "-40", "right": "0"},
                ),
                html.Div(
                    [  # subpage 2
                        # Row 1 (Header)
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H5(
                                            "Goldman Sachs Strategic  Absolute Return Bond II Portfolio"
                                        ),
                                        html.H6(
                                            "A sub-fund of Goldman Sachs Funds, SICAV",
                                            style={"color": "#7F90AC"},
                                        ),
                                    ],
                                    className="nine columns padded",
                                ),
                                html.Div(
                                    [
                                        html.H1(
                                            [
                                                html.Span(
                                                    "03", style={"opacity": "0.5"}
                                                ),
                                                html.Span("17"),
                                            ]
                                        ),
                                        html.H6("Monthly Fund Update"),
                                    ],
                                    className="three columns gs-header gs-accent-header padded",
                                    style={"float": "right"},
                                ),
                            ],
                            className="row gs-header gs-text-header",
                        ),
                        # Row 2
                        html.Div(
                            [
                                # Data tables on this page:
                                # ---
                                # df_fund_info = pd.read_csv('https://plot.ly/~jackp/17544/.csv')
                                # df_fund_characteristics = pd.read_csv('https://plot.ly/~jackp/17542/.csv')
                                # df_fund_facts = pd.read_csv('https://plot.ly/~jackp/17540/.csv')
                                # df_bond_allocation = pd.read_csv('https://plot.ly/~jackp/17538/')
                                # Column 1
                                html.Div(
                                    [
                                        html.H6(
                                            "Financial Information",
                                            className="gs-header gs-text-header padded",
                                        ),
                                        html.Table(make_dash_table(df_fund_info)),
                                        html.H6(
                                            "Fund Characteristics",
                                            className="gs-header gs-text-header padded",
                                        ),
                                        html.Table(
                                            make_dash_table(df_fund_characteristics)
                                        ),
                                        html.H6(
                                            "Fund Facts",
                                            className="gs-header gs-text-header padded",
                                        ),
                                        html.Table(make_dash_table(df_fund_facts)),
                                    ],
                                    className="four columns",
                                ),
                                # Column 2
                                html.Div(
                                    [
                                        html.H6(
                                            "Sector Allocation (%)",
                                            className="gs-header gs-table-header padded",
                                        ),
                                        html.Iframe(
                                            src="https://plot.ly/~jackp/17560.embed?modebar=false&link=false&autosize=true",
                                            style={"border": "0"},
                                            width="100%",
                                            height="300",
                                        ),
                                        html.H6(
                                            "Country Bond Allocation (%)",
                                            className="gs-header gs-table-header padded",
                                        ),
                                        html.Table(make_dash_table(df_bond_allocation)),
                                    ],
                                    className="four columns",
                                ),
                                # Column 3
                                html.Div(
                                    [
                                        html.H6(
                                            "Top 10 Currency Weights (%)",
                                            className="gs-header gs-table-header padded",
                                        ),
                                        html.Iframe(
                                            src="https://plot.ly/~jackp/17555.embed?modebar=false&link=false&autosize=true",
                                            style={"border": "0"},
                                            width="100%",
                                            height="300",
                                        ),
                                        html.H6(
                                            "Credit Allocation (%)",
                                            className="gs-header gs-table-header padded",
                                        ),
                                        html.Iframe(
                                            src="https://plot.ly/~jackp/17557.embed?modebar=false&link=false&autosize=true",
                                            style={"border": "0"},
                                            width="100%",
                                            height="300",
                                        ),
                                    ],
                                    className="four columns",
                                ),
                            ],
                            className="row",
                        ),
                    ],
                    className="subpage",
                ),
            ],
            className="page",
        ),
    ]
)

if "DYNO" in os.environ:
    app.scripts.append_script(
        {
            "external_url": "https://cdn.rawgit.com/chriddyp/ca0d8f02a1659981a0ea7f013a378bbd/raw/e79f3f789517deec58f41251f7dbb6bee72c44ab/plotly_ga.js"
        }
    )


external_js = [
    "https://code.jquery.com/jquery-3.2.1.min.js",
    "https://cdn.rawgit.com/plotly/dash-app-stylesheets/a3401de132a6d0b652ba11548736b1d1e80aa10d/dash-goldman-sachs-report-js.js",
]

for js in external_js:
    app.scripts.append_script({"external_url": js})


if __name__ == "__main__":
    app.server.run()
