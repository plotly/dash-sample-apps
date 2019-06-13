# -*- coding: utf-8 -*-
import json
import base64
import datetime
import requests

import pandas as pd
import flask
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.plotly as py
import plotly.graph_objs as go

from dash.dependencies import Input, Output, State
from plotly import tools


server = flask.Flask(__name__)
app = dash.Dash(__name__, server=server, meta_tags=[{"name": "viewport", "content": "width=device-width"}])


# Loading historical tick data
EURUSD = pd.read_csv("pairs/EURUSD.csv", index_col=1, parse_dates=["Date"])
USDJPY = pd.read_csv("pairs/USDJPY.csv", index_col=1, parse_dates=["Date"])
GBPUSD = pd.read_csv("pairs/GBPUSD.csv", index_col=1, parse_dates=["Date"])
USDCHF = pd.read_csv("pairs/USDCHF.csv", index_col=1, parse_dates=["Date"])

# Currency pairs
currencies = ["EURUSD", "USDCHF", "USDJPY", "GBPUSD"]

# Generate HTML table to display the news
# Display a maximum of ten rows
def generate_news_table(dataframe, max_rows=10):
    return html.Div(
        children=[
            html.P("Headlines"),
            html.P("Last update : " + datetime.datetime.now().strftime("%H:%M:%S")),
            html.Table(
                className="table-news",
                children=[
                    html.Tr(
                        [
                            html.Td(
                                children=[
                                    html.A(
                                        className="td-link",
                                        children=dataframe.iloc[i]["title"],
                                        href=dataframe.iloc[i]["url"],
                                        target="_blank",
                                    )
                                ]
                            )
                        ]
                    )
                    for i in range(min(len(dataframe), max_rows))
                ],
            ),
        ]
    )


# API Call to update news
def update_news():
    r = requests.get(
        "https://newsapi.org/v2/top-headlines?sources=bbc-news&apiKey=da8e2e705b914f9f86ed2e9692e66012"
    )
    json_data = r.json()["articles"]
    df = pd.DataFrame(json_data)
    df = pd.DataFrame(df[["title", "url"]])
    return generate_news_table(df)


# Get Current Time and create header
def get_header(t=datetime.datetime.now()):
    return html.Div(
        className="header-bid-ask",
        children=[
            html.P(t.strftime("%H:%M:%S"), id="live_clock", className="four columns"),
            html.P("Bid", className="four columns"),
            html.P("Ask", className="four columns"),
        ],
    )


# color of Bid & Ask rates
def get_color(a, b):
    if a == b:
        return "white"
    elif a > b:
        return "#45df7e"
    else:
        return "#da5657"


# returns row of given index for given currency pair dataset
# for modal returns ten previous rows
def get_ask_bid(currency_pair, index, modal=False):
    if modal == False:
        return globals()[currency_pair].iloc[index]
    return globals()[currency_pair].iloc[index - 10 : index]


# returns dataset row with nearest datetime to current time
def first_ask_bid(currency_pair, t):
    t = t.replace(year=2016, month=1, day=5)
    items = globals()[currency_pair]
    dates = items.index.to_pydatetime()
    index = min(dates, key=lambda x: abs(x - t))
    df_row = items.loc[index]
    int_index = items.index.get_loc(index)
    return [df_row, int_index]  # returns dataset row and index of row


# Creates HTML Bid and Ask (Buy/Sell buttons)
def get_row(data):
    index = data[1]
    current_row = data[0]
    return html.Details(
        children=[
            html.Summary(
                className="row",
                children=[
                    html.Div(
                        children=[
                            html.P(
                                current_row[0],  # currency pair name
                                id=current_row[0],
                                className="four columns",
                                style={"textAlign": "center"},
                            ),
                            html.P(
                                current_row[1].round(5),  # Bid value
                                id=current_row[0] + "bid",
                                className="four columns",
                                style={"textAlign": "center"},
                            ),
                            html.P(
                                current_row[2].round(5),  # Ask value
                                className="four columns",
                                id=current_row[0] + "ask",
                                style={"textAlign": "center"},
                            ),
                            html.Div(
                                index,
                                id=current_row[0]
                                + "index",  # we save index of row in hidden div
                                style={"display": "none"},
                            ),
                        ],
                        id=current_row[0] + "row",
                        className="row eleven columns",
                        style={"height": "25", "float": "right"},
                    )
                ],
            ),
            html.Div(
                className="row",
                children=[
                    html.Div(
                        className="button-buy-sell-chart",
                        children=[
                            html.Button(
                                id=current_row[0] + "Buy",
                                children="Buy/Sell",
                                n_clicks=0,
                            )
                        ]
                    ),
                    html.Div(
                        className="button-buy-sell-chart-right",
                        children=[
                            html.Button(
                                id=current_row[0] + "Button_chart",
                                children="Chart",
                                n_clicks=1 if current_row[0] in ["EURUSD", "USDCHF"] else 0,
                            )
                        ] 
                    )
                ],
            ),
        ],
        id=current_row[0] + "row_div",
        n_clicks=0,
    )


# replace ask_bid row without buttons
def replace_row(currency_pair, index, bid, ask):
    index = index + 1  # index of new data row
    new_row = (
        get_ask_bid(currency_pair, index)
        if index != len(globals()[currency_pair])
        else first_ask_bid(currency_pair, datetime.datetime.now())
    )  # if not the end of the dataset we retrieve next dataset row

    return [
        html.P(
            currency_pair,  # currency pair name
            id=currency_pair,
            className="four columns",
            style={"textAlign": "center"},
        ),
        html.P(
            new_row[1].round(5),  # Bid value
            id=new_row[0] + "bid",
            className="four columns",
            style={"textAlign": "center", "color": get_color(new_row[1], bid)},
        ),
        html.P(
            new_row[2].round(5),  # Ask value
            className="four columns",
            id=new_row[0] + "ask",
            style={"textAlign": "center", "color": get_color(new_row[2], ask)},
        ),
        html.Div(
            index, id=currency_pair + "index", style={"display": "none"}
        ),  # save index in hidden div
    ]


# Returns Top cell bar for header area
def get_top_bar_cell(cellTitle, cellValue, color="white"):
    return html.Div(
        className="two columns",
        children=[
            html.P(className="p-top-bar", children=cellTitle),
            html.P(cellValue, id=cellTitle),
        ],
    )


# Returns HTML Top Bar for app layout
def get_top_bar(
    balance=50000.00, equity=50000, margin=0, fm=50000, m_level="%", open_pl=0
):
    color_open_pl = get_color(float(open_pl), 0)
    return [
        get_top_bar_cell("Balance", balance),
        get_top_bar_cell("Equity", equity),
        get_top_bar_cell("Margin", margin),
        get_top_bar_cell("Free Margin", fm),
        get_top_bar_cell("Margin Level", m_level),
        get_top_bar_cell("Open P/L", open_pl, color=color_open_pl),
    ]


# returns currency pair OHLC data for given period
def get_OHLC_data(currency_pair, period="5Min"):
    data_frame = globals()[currency_pair]
    t = datetime.datetime.now()
    data = data_frame.loc[
        : t.strftime(
            "2016-01-05 %H:%M:%S"
        )  # all the data from the beginning until current time
    ]
    data_bid = data["Bid"]
    data_bid = data_bid.resample(period).ohlc()
    return data_bid


#######STUDIES TRACES######

# Moving average
def moving_average_trace(df, fig):
    df2 = df.rolling(window=5).mean()
    trace = go.Scatter(
        x=df2.index, y=df2["close"], mode="lines", showlegend=False, name="MA"
    )
    fig.append_trace(trace, 1, 1)  # plot in first row
    return fig


# Exponential moving average
def e_moving_average_trace(df, fig):
    df2 = df.rolling(window=20).mean()
    trace = go.Scatter(
        x=df2.index, y=df2["close"], mode="lines", showlegend=False, name="EMA"
    )
    fig.append_trace(trace, 1, 1)  # plot in first row
    return fig


# Bollinger Bands
def bollinger_trace(df, fig, window_size=10, num_of_std=5):
    price = df["close"]
    rolling_mean = price.rolling(window=window_size).mean()
    rolling_std = price.rolling(window=window_size).std()
    upper_band = rolling_mean + (rolling_std * num_of_std)
    lower_band = rolling_mean - (rolling_std * num_of_std)

    trace = go.Scatter(
        x=df.index, y=upper_band, mode="lines", showlegend=False, name="BB_upper"
    )

    trace2 = go.Scatter(
        x=df.index, y=rolling_mean, mode="lines", showlegend=False, name="BB_mean"
    )

    trace3 = go.Scatter(
        x=df.index, y=lower_band, mode="lines", showlegend=False, name="BB_lower"
    )

    fig.append_trace(trace, 1, 1)  # plot in first row
    fig.append_trace(trace2, 1, 1)  # plot in first row
    fig.append_trace(trace3, 1, 1)  # plot in first row
    return fig


# Accumulation Distribution
def accumulation_trace(df):
    df["volume"] = ((df["close"] - df["low"]) - (df["high"] - df["close"])) / (
        df["high"] - df["low"]
    )
    trace = go.Scatter(
        x=df.index, y=df["volume"], mode="lines", showlegend=False, name="Accumulation"
    )
    return trace


# Commodity Channel Index
def cci_trace(df, ndays=5):
    TP = (df["high"] + df["low"] + df["close"]) / 3
    CCI = pd.Series(
        (TP - TP.rolling(window=10, center=False).mean())
        / (0.015 * TP.rolling(window=10, center=False).std()),
        name="cci",
    )
    trace = go.Scatter(x=df.index, y=CCI, mode="lines", showlegend=False, name="CCI")
    return trace


# Price Rate of Change
def roc_trace(df, ndays=5):
    N = df["close"].diff(ndays)
    D = df["close"].shift(ndays)
    ROC = pd.Series(N / D, name="roc")
    trace = go.Scatter(x=df.index, y=ROC, mode="lines", showlegend=False, name="ROC")
    return trace


# Stochastic oscillator %K
def stoc_trace(df):
    SOk = pd.Series((df["close"] - df["low"]) / (df["high"] - df["low"]), name="SO%k")
    trace = go.Scatter(x=df.index, y=SOk, mode="lines", showlegend=False, name="SO%k")
    return trace


# Momentum
def mom_trace(df, n=5):
    M = pd.Series(df["close"].diff(n), name="Momentum_" + str(n))
    trace = go.Scatter(x=df.index, y=M, mode="lines", showlegend=False, name="MOM")
    return trace


# Pivot points
def pp_trace(df, fig):
    PP = pd.Series((df["high"] + df["low"] + df["close"]) / 3)
    R1 = pd.Series(2 * PP - df["low"])
    S1 = pd.Series(2 * PP - df["high"])
    R2 = pd.Series(PP + df["high"] - df["low"])
    S2 = pd.Series(PP - df["high"] + df["low"])
    R3 = pd.Series(df["high"] + 2 * (PP - df["low"]))
    S3 = pd.Series(df["low"] - 2 * (df["high"] - PP))
    trace = go.Scatter(x=df.index, y=PP, mode="lines", showlegend=False, name="PP")
    trace1 = go.Scatter(x=df.index, y=R1, mode="lines", showlegend=False, name="R1")
    trace2 = go.Scatter(x=df.index, y=S1, mode="lines", showlegend=False, name="S1")
    trace3 = go.Scatter(x=df.index, y=R2, mode="lines", showlegend=False, name="R2")
    trace4 = go.Scatter(x=df.index, y=S2, mode="lines", showlegend=False, name="S2")
    trace5 = go.Scatter(x=df.index, y=R3, mode="lines", showlegend=False, name="R3")
    trace6 = go.Scatter(x=df.index, y=S3, mode="lines", showlegend=False, name="S3")
    fig.append_trace(trace, 1, 1)
    fig.append_trace(trace1, 1, 1)
    fig.append_trace(trace2, 1, 1)
    fig.append_trace(trace3, 1, 1)
    fig.append_trace(trace4, 1, 1)
    fig.append_trace(trace5, 1, 1)
    fig.append_trace(trace6, 1, 1)
    return fig


## MAIN CHART TRACES (STYLE tab)
def line_trace(df):
    trace = go.Scatter(
        x=df.index, y=df["close"], mode="lines", showlegend=False, name="line"
    )
    return trace


def area_trace(df):
    trace = go.Scatter(
        x=df.index, y=df["close"], showlegend=False, fill="toself", name="area"
    )
    return trace


def bar_trace(df):
    return go.Ohlc(
        x=df.index,
        open=df["open"],
        high=df["high"],
        low=df["low"],
        close=df["close"],
        increasing=dict(line=dict(color="#888888")),
        decreasing=dict(line=dict(color="#888888")),
        showlegend=False,
        name="bar",
    )


def colored_bar_trace(df):
    return go.Ohlc(
        x=df.index,
        open=df["open"],
        high=df["high"],
        low=df["low"],
        close=df["close"],
        showlegend=False,
        name="colored bar",
    )


def candlestick_trace(df):
    return go.Candlestick(
        x=df.index,
        open=df["open"],
        high=df["high"],
        low=df["low"],
        close=df["close"],
        increasing=dict(line=dict(color="#00ff00")),
        decreasing=dict(line=dict(color="white")),
        showlegend=False,
        name="candlestick",
    )


# For buy/sell modal
def ask_modal_trace(currency_pair, index):
    df = get_ask_bid(currency_pair, index, True)  # returns ten rows
    return go.Scatter(x=df.index, y=df["Ask"], mode="lines", showlegend=False)


# For buy/sell modal
def bid_modal_trace(currency_pair, index):
    df = get_ask_bid(currency_pair, index, True)  # returns ten rows
    return go.Scatter(x=df.index, y=df["Bid"], mode="lines", showlegend=False)


# returns modal figure for a currency pair
def get_modal_fig(currency_pair, index):
    fig = tools.make_subplots(
        rows=2, shared_xaxes=True, shared_yaxes=False, cols=1, print_grid=False
    )

    fig.append_trace(ask_modal_trace(currency_pair, index), 1, 1)
    fig.append_trace(bid_modal_trace(currency_pair, index), 2, 1)

    fig["layout"]["height"] = 400
    fig["layout"]["autosize"] = True
    fig["layout"]["margin"] = {"b": 0, "r": 5, "l": 50, "t": 5}
    fig["layout"].update(paper_bgcolor="#22252b", plot_bgcolor="#22252b")
    return fig


# returns graph figure
def get_fig(currency_pair, ask, bid, type_trace, studies, period):
    df = get_OHLC_data(currency_pair, period)
    subplot_traces = [  # first row traces
        "accumulation_trace",
        "cci_trace",
        "roc_trace",
        "stoc_trace",
        "mom_trace",
    ]
    selected_subplots_studies = []
    selected_first_row_studies = []
    row = 1  # number of subplots

    if studies:
        for study in studies:
            if study in subplot_traces:
                row += 1  # increment number of rows only if the study needs a subplot
                selected_subplots_studies.append(study)
            else:
                selected_first_row_studies.append(study)

    fig = tools.make_subplots(
        rows=row,
        shared_xaxes=True,
        shared_yaxes=True,
        cols=1,
        print_grid=False,
        vertical_spacing=0.12,
    )

    # Add main trace (style) to figure
    fig.append_trace(globals()[type_trace](df), 1, 1)

    # Add trace(s) on fig's first row
    for study in selected_first_row_studies:
        fig = globals()[study](df, fig)

    row = 1
    # Plot trace on new row
    for study in selected_subplots_studies:
        row += 1
        fig.append_trace(globals()[study](df), row, 1)

    fig["layout"]["margin"] = {"b": 50, "r": 5, "l": 50, "t": 5}
    fig["layout"]["autosize"] = True
    fig["layout"]["height"] = 400
    fig["layout"]["xaxis"]["rangeslider"]["visible"] = False
    fig["layout"]["xaxis"]["tickformat"] = "%H:%M"
    fig["layout"].update(paper_bgcolor="#22252b", plot_bgcolor="#22252b")
    return fig


# updates figure
def replace_fig(currency_pair, ask, bid, type_trace, period, old_fig, studies):
    fig = get_fig(currency_pair, ask, bid, type_trace, studies, period)
    fig["layout"]["xaxis"]["range"] = old_fig["layout"]["xaxis"]["range"]
    return fig


# returns chart div
def chart_div(pair):
    return html.Div(
        id=pair + "graph_div",
        children=[
            # Menu for Currency Graph
            html.Div(
                id=pair + "menu",
                className="not_visible",
                children=[
                    # stores current menu tab
                    html.Div(
                        id=pair + "menu_tab",
                        children=["Studies"],
                        style={"display": "none"},
                    ),
                    html.Span(
                        "Style",
                        id=pair + "style_header",
                        className="btn-style",
                        n_clicks_timestamp=2,
                    ),
                    html.Span(
                        "Studies",
                        id=pair + "studies_header",
                        n_clicks_timestamp=1,
                        style={"textDecoration": "none", "cursor": "pointer"},
                    ),
                    html.Div(
                        html.Div(
                            dcc.Checklist(
                                id=pair + "studies",
                                options=[
                                    {
                                        "label": "Accumulation/D",
                                        "value": "accumulation_trace",
                                    },
                                    {
                                        "label": "Bollinger bands",
                                        "value": "bollinger_trace",
                                    },
                                    {"label": "MA", "value": "moving_average_trace"},
                                    {"label": "EMA", "value": "e_moving_average_trace"},
                                    {"label": "CCI", "value": "cci_trace"},
                                    {"label": "ROC", "value": "roc_trace"},
                                    {"label": "Pivot points", "value": "pp_trace"},
                                    {
                                        "label": "Stochastic oscillator",
                                        "value": "stoc_trace",
                                    },
                                    {
                                        "label": "Momentum indicator",
                                        "value": "mom_trace",
                                    },
                                ],
                                values=[],
                            ),
                            style={"marginTop": "30", "textAlign": "left"},
                        ),
                        id=pair + "studies_tab",
                        style={"display": "none"},
                    ),
                    html.Div(
                        id=pair + "style_tab",
                        children=[
                            dcc.RadioItems(
                                id=pair + "chart_type",
                                options=[
                                    {
                                        "label": "candlestick",
                                        "value": "candlestick_trace",
                                    },
                                    {"label": "line", "value": "line_trace"},
                                    {"label": "mountain", "value": "area_trace"},
                                    {"label": "bar", "value": "bar_trace"},
                                    {
                                        "label": "colored bar",
                                        "value": "colored_bar_trace",
                                    },
                                ],
                                value="colored_bar_trace",
                            )
                        ],
                    ),
                ],
            ),

            # Chart Top Bar
            html.Div(
                className="row chart-top-bar",
                children=[
                    html.Span(
                        id=pair + "menu_button",
                        children=f"{pair} ☰",
                        n_clicks=0,
                        style={'display':'inline-block'}
                    ),
                    # Dropdown and close button float right
                    html.Div(
                        className="graph-top-right",
                        children=[
                            dcc.Dropdown(
                                className="dropdown-period",
                                id=pair + "dropdown_period",
                                options=[
                                    {"label": "5 min", "value": "5Min"},
                                    {"label": "15 min", "value": "15Min"},
                                    {"label": "30 min", "value": "30Min"},
                                ],
                                value="15Min",
                                clearable=False,
                                style={'display':'inline-block'}
                            ),
                            html.Span(
                                id=pair + "close",
                                className="graph-close row",
                                children="×",
                                n_clicks=0,
                                style={'display':'inline-block'}
                            ),
                        ]
                    ),
                ],
            ),
            # Graph div
            html.Div(
                dcc.Graph(
                    id=pair + "chart",
                    config={"displayModeBar": False, "scrollZoom": True},
                    style={"width": "100%", "height": "100%"},
                ),
                id=pair + "graph",
                style={"height": "100%"},
            ),
        ],
        className="",
        style={"display": "none"},
    )


# bottom panel for orders
def bottom_panel():
    return html.Div(
        id="bottom_panel",
        className="row div-bottom-panel",
        children=[
            dcc.Dropdown(
                id="dropdown_positions",
                className="bottom-dropdown",
                options=[
                    {"label": "Open Positions", "value": "open"},
                    {"label": "Closed Positions", "value": "closed"},
                ],
                value="open",
                clearable=False,
            ),
            dcc.Dropdown(
                id="closable_orders",
                className="bottom-dropdown",
                placeholder="Close order",
            ),
            html.Div(
                id="bottom_content",
                className="row",
                children=[html.Table(id="orders_table")],
            ),
        ],
    )


# returns modal Buy/Sell
def modal(pair):
    return html.Div(
        id=pair + "modal",
        className="modal",
        style={"display": "none"},
        children=[
            html.Div(
                className="modal-content",
                children=[
                    html.Span(
                        "×",
                        id=pair + "closeModal",
                        style={
                            "float": "right",
                            "cursor": "pointer",
                            "marginTop": "0",
                            "marginBottom": "10",
                        },
                    ),
                    html.Span(
                        pair,
                        id="modal" + pair,
                        style={"marginBottom": "10", "color": "#45df7e"},
                    ),
                    # row div with two div
                    html.Div(
                        className="row",
                        children=[
                            # graph div
                            html.Div(
                                className="six columns",
                                children=[
                                    dcc.Graph(
                                        id=pair + "modal_graph",
                                        config={"displayModeBar": False},
                                    )
                                ]
                            ),
                            # order values div
                            html.Div(
                                className="six columns",
                                children=[
                                    html.Div(
                                        children=[
                                            html.P("Volume"),
                                            dcc.Input(
                                                id=pair + "volume",
                                                className="modal-input",
                                                type="number",
                                                value=0.1,
                                                min=0,
                                                step=0.1,
                                            ),
                                        ],
                                        #style={"marginBottom": "5"},
                                    ),
                                    html.Div(
                                        children=[
                                            html.P("Type"),
                                            dcc.RadioItems(
                                                id=pair + "trade_type",
                                                options=[
                                                    {"label": "Buy", "value": "buy"},
                                                    {"label": "Sell", "value": "sell"},
                                                ],
                                                value="buy",
                                                labelStyle={"display": "inline-block"},
                                            ),
                                        ],
                                        #style={"marginBottom": "5"},
                                    ),
                                    html.Div(
                                        children=[
                                            html.P("SL TPS"),
                                            dcc.Input(
                                                id=pair + "SL",
                                                type="number",
                                                min=0,
                                                step=1,
                                            ),
                                        ],
                                        #style={"marginBottom": "5"},
                                    ),
                                    html.Div(
                                        children=[
                                            html.P("TP TPS"),
                                            dcc.Input(
                                                id=pair + "TP",
                                                type="number",
                                                min=0,
                                                step=1,
                                            ),
                                        ],
                                        #style={"marginBottom": "5"},
                                    ),
                                ],
                            ),
                        ]
                    ),
                    html.Div(
                        html.Button("Order", id=pair + "button_order", n_clicks=0),
                        style={"textAlign": "center", "marginTop": "12"},
                    ),
                ]
            )
        ],
    )


# Dash App Layout
app.layout = html.Div(
    className="row",
    children=[
        # Interval component for live clock
        dcc.Interval(id="interval", interval=1 * 1000, n_intervals=0),
        # Interval component for ask bid updates
        dcc.Interval(id="i_bis", interval=1 * 2000, n_intervals=0),
        # Interval component for graph updates
        dcc.Interval(id="i_tris", interval=1 * 5000, n_intervals=0),
        # Interval component for graph updates
        dcc.Interval(id="i_news", interval=1 * 60000, n_intervals=0),
        # Left Panel Div
        html.Div(
            className="three columns div-left-panel",
            children=[
                html.Div(
                    className="div-info",
                    children=[
                        html.Img(className="logo", src="assets/dash-logo.png"),
                        html.P(
                            """
                            This app continually queries csv files and updates Ask and Bid prices 
                            for major currency pairs as well as Stock Charts. You can also virtually 
                            buy and sell stocks and see the profit updates.
                            """
                        ),
                    ],
                ),
                html.Div(get_header()),
                html.Div(
                    id="pairs",
                    className="div-bid-ask",
                    children=[
                        get_row(first_ask_bid(pair, datetime.datetime.now()))
                        for pair in currencies
                    ],
                ),
                html.Div(
                    className="div-news",
                    children=[html.Div(id="news", children=update_news())],
                ),
            ],
        ),
        # Right Panel Div
        html.Div(
            className="nine columns div-right-panel",
            children=[
                # Top Bar Div - Displays Balance, Equity, ... , Open P/L
                html.Div(
                    id="top_bar", className="row div-top-bar", children=get_top_bar()
                ),
                # Charts Div
                html.Div(
                    id="charts",
                    className="row",
                    children=[chart_div(pair) for pair in currencies],
                ),
                # Panel for orders
                bottom_panel(),
            ],
        ),
        # Hidden div that stores all clicked charts (EURUSD, USDCHF, etc.)
        html.Div(id="charts_clicked", style={"display": "none"}),
        # Hidden div for each pair that stores orders
        html.Div(
            children=[
                html.Div(id=pair + "orders", style={"display": "none"})
                for pair in currencies
            ]
        ),
        html.Div([modal(pair) for pair in currencies]),
        # Hidden Div that stores all orders
        html.Div(id="orders", style={"display": "none"}),
    ],
)

# Dynamic Callbacks

# Replace currency pair row
def generate_ask_bid_row_callback(pair):
    def output_callback(n, i, bid, ask):
        return replace_row(pair, int(i), float(bid), float(ask))

    return output_callback


# returns string containing clicked charts
def generate_chart_button_callback():
    def chart_button_callback(*args):
        pairs = ""
        for i in range(len(currencies)):
            if args[i] > 0:
                pair = currencies[i]
                if pairs:
                    pairs = pairs + "," + pair
                else:
                    pairs = pair
        return pairs

    return chart_button_callback


# Function to update Graph Figure
def generate_figure_callback(pair):
    def chart_fig_callback(n_i, p, t, s, pairs, a, b, old_fig):
        """
        n_i : Number of intervals for Ask/Bid updates
        p: Period chosen for chart (5, 15 or 30 min)
        t: Type of chart (candlestick, line, mountain, etc.)
        s: Type of study ()
        pairs: Currency pairs that were clicked to be displayed
        a: ask
        b: bid
        old_fig: 
        """

        if pairs is None:
            return {"layout": {}, "data": {}}

        pairs = pairs.split(",")  # list of displayed divs
        if pair not in pairs:
            return {
                "layout": {},
                "data": [],
            }  # we only update figure when the div is displayed

        if old_fig is None or old_fig == {"layout": {}, "data": []}:
            return get_fig(pair, a, b, t, s, p)

        # return replace_fig(pair, a, b, t, p, old_fig, s)

        return get_fig(pair, a, b, t, s, p)

    return chart_fig_callback


def generate_close_graph_callback():
    def close_callback(n, n2):
        if n == 0:
            if n2 == 1:
                return 1
            return 0
        return 0

    return close_callback


def generate_open_close_menu_callback():
    def open_close_menu(n, className):
        if n == 0:
            return "not_visible"
        if className == "visible":
            return "not_visible"
        else:
            return "visible"

    return open_close_menu


# Updates hidden div that stores the last clicked menu tab
def generate_active_menu_tab_callback():
    def update_current_tab_name(n_style, n_studies):
        if n_style >= n_studies:
            return "Style"
        return "Studies"

    return update_current_tab_name


# Show/hide STUDIES menu for chart
def generate_studies_content_tab_callback():
    def studies_tab(current_tab):
        if current_tab == "Studies":
            return {"display": "block", "textAlign": "left", "marginTop": "30"}
        return {"display": "none"}

    return studies_tab


# Show/hide STYLE menu for chart
def generate_style_content_tab_callback():
    def style_tab(current_tab):
        if current_tab == "Style":
            return {"display": "block", "textAlign": "left", "marginTop": "30"}
        return {"display": "none"}

    return style_tab


# Open Modal
def generate_modal_open_callback():
    def open_modal(n):
        if n > 0:
            return {"display": "block"}
        else:
            return {"display": "none"}

    return open_modal


# Close Modal
def generate_modal_close_callback():
    def close_modal(n, n2):
        return 0

    return close_modal


# set modal SL value to none
def generate_clean_sl_callback():
    def clean_sl(n):
        return 0

    return clean_sl


# set modal SL value to none
def generate_clean_tp_callback():
    def clean_tp(n):
        return 0

    return clean_tp


# Create figure for Buy/Sell Modal
def generate_modal_figure_callback(pair):
    def figure_modal(index, n, old_fig):
        if (n == 0 and old_fig is None) or n == 1:
            return get_modal_fig(pair, index)
        return old_fig  # avoid to compute new figure when the modal is hidden

    return figure_modal


# updates the pair orders div
def generate_order_button_callback(pair):
    def order_callback(n, vol, type_order, sl, tp, pair_orders, ask, bid):
        if n > 0:
            t = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            l = [] if pair_orders is None else json.loads(pair_orders)
            price = bid if type_order == "sell" else ask

            if tp != 0:
                tp = (
                    price + tp * 0.001
                    if tp != 0 and pair[3:] == "JPY"
                    else price + tp * 0.00001
                )

            if sl != 0:
                sl = price - sl * 0.001 if pair[3:] == "JPY" else price + sl * 0.00001

            order = {
                "id": pair + str(len(l)),
                "time": t,
                "type": type_order,
                "volume": vol,
                "symbol": pair,
                "tp": tp,
                "sl": sl,
                "price": price,
                "profit": 0.00,
                "status": "open",
            }
            l.append(order)

            return json.dumps(l)

        return json.dumps([])

    return order_callback


def update_orders(orders, current_bids, current_asks, id_to_close):
    for order in orders:
        if order["status"] == "open":
            type_order = order["type"]
            current_bid = current_bids[currencies.index(order["symbol"])]
            current_ask = current_asks[currencies.index(order["symbol"])]

            profit = (
                order["volume"]
                * 100000
                * ((current_bid - order["price"]) / order["price"])
                if type_order == "buy"
                else (
                    order["volume"]
                    * 100000
                    * ((order["price"] - current_ask) / order["price"])
                )
            )

            order["profit"] = "%.2f" % profit
            price = current_bid if order["type"] == "buy" else current_ask

            if order["id"] == id_to_close:
                order["status"] = "closed"
                order["close Time"] = datetime.datetime.now().strftime(
                    "%Y-%m-%d %H:%M:%S"
                )
                order["close Price"] = price

            if order["tp"] != 0 and price >= order["tp"]:
                order["status"] = "closed"
                order["close Time"] = datetime.datetime.now().strftime(
                    "%Y-%m-%d %H:%M:%S"
                )
                order["close Price"] = price

            if order["sl"] != 0 and order["sl"] >= price:
                order["status"] = "closed"
                order["close Time"] = datetime.datetime.now().strftime(
                    "%Y-%m-%d %H:%M:%S"
                )
                order["close Price"] = price

    return orders


def generate_update_orders_div_callback():
    def update_orders_callback(*args):
        orders = []
        current_orders = args[-1]
        close_id = args[-2]
        args = args[:-2]  # contains list of orders for each pair + asks + bids
        len_args = len(args)
        current_bids = args[len_args // 3 : 2 * len_args]
        current_asks = args[2 * len_args // 3 : len_args]
        args = args[: len_args // 3]
        ids = []

        if current_orders is not None:
            orders = json.loads(current_orders)
            for order in orders:
                ids.append(
                    order["id"]  # ids that allready have been added to current orders
                )

        for list_order in args:  # each currency pair has its list of orders
            if list_order != "[]":
                list_order = json.loads(list_order)
                for order in list_order:
                    if order["id"] not in ids:  # only add new orders
                        orders.append(order)
        if len(orders) == 0:
            return None

        # we update status and profit of orders
        orders = update_orders(orders, current_bids, current_asks, close_id)

        return json.dumps(orders)

    return update_orders_callback


def generate_show_hide_graph_div_callback(pair):
    def show_hide_graph_callback(charts_clicked):
        if charts_clicked is None:
            return {"display": "none"}

        charts_clicked = charts_clicked.split(",")[
            :4
        ]  # only allow to display four charts
        len_list = len(charts_clicked)

        for i in range(len_list):
            if charts_clicked[i] == pair:
                style = {"position": "relative", "float": "left", "overflow": "hidden"}

                if i == 0 or (i == 2 and len_list == 4):
                    style["marginLeft"] = "0px"  # avoid div to overlap

                style["height"] = "50%" if len_list == 4 else "100%"

                return style

        return {"display": "none"}

    return show_hide_graph_callback


# Resize pair div according to the number of charts displayed
def generate_size_graph_div_callback(pair):
    def size_graph_div_callback(charts_clicked):
        if charts_clicked is None:
            return ""

        charts_clicked = charts_clicked.split(",")  # [:4] max of 4 graph
        len_list = len(charts_clicked)
        if pair not in charts_clicked:
            return ""

        width = (
            "six columns"
            if len_list % 2 == 0
            else "twelve columns"
            if len_list == 1
            else "four columns"
        )
        return width

    return size_graph_div_callback


app.config.supress_callback_exceptions = True


# Loop through all currencies
for pair in currencies:

    # Callback for style of div for graphs
    app.callback(
        Output(pair + "graph_div", "style"), [Input("charts_clicked", "children")]
    )(generate_show_hide_graph_div_callback(pair))

    # Callback for class of div for graphs
    app.callback(
        Output(pair + "graph_div", "className"), [Input("charts_clicked", "children")]
    )(generate_size_graph_div_callback(pair))

    # Callback to update the actual graph
    app.callback(
        Output(pair + "chart", "figure"),
        [
            Input("i_tris", "n_intervals"),
            Input(pair + "dropdown_period", "value"),
            Input(pair + "chart_type", "value"),
            Input(pair + "studies", "values"),
            Input("charts_clicked", "children"),
        ],
        [
            State(pair + "ask", "children"),
            State(pair + "bid", "children"),
            State(pair + "chart", "figure"),
        ],
    )(generate_figure_callback(pair))

    # updates the ask and bid prices
    app.callback(
        Output(pair + "row", "children"),
        [Input("i_bis", "n_intervals")],
        [
            State(pair + "index", "children"),
            State(pair + "bid", "children"),
            State(pair + "ask", "children"),
        ],
    )(generate_ask_bid_row_callback(pair))

    # close graph by setting to 0 n_clicks property
    app.callback(
        Output(pair + "Button_chart", "n_clicks"),
        [Input(pair + "close", "n_clicks")],
        [State(pair + "Button_chart", "n_clicks")],
    )(generate_close_graph_callback())

    # show or hide graph menu
    app.callback(
        Output(pair + "menu", "className"),
        [Input(pair + "menu_button", "n_clicks")],
        [State(pair + "menu", "className")],
    )(generate_open_close_menu_callback())

    # stores in hidden div name of clicked tab name
    app.callback(
        Output(pair + "menu_tab", "children"),
        [
            Input(pair + "style_header", "n_clicks_timestamp"),
            Input(pair + "studies_header", "n_clicks_timestamp"),
        ],
    )(generate_active_menu_tab_callback())

    # hide/show STYLE tab content if clicked or not
    app.callback(
        Output(pair + "style_tab", "style"), [Input(pair + "menu_tab", "children")]
    )(generate_style_content_tab_callback())

    # hide/show MENU tab content if clicked or not
    app.callback(
        Output(pair + "studies_tab", "style"), [Input(pair + "menu_tab", "children")]
    )(generate_studies_content_tab_callback())

    #####
    # show modal
    app.callback(Output(pair + "modal", "style"), [Input(pair + "Buy", "n_clicks")])(
        generate_modal_open_callback()
    )

    # set modal value SL to O
    app.callback(Output(pair + "SL", "value"), [Input(pair + "Buy", "n_clicks")])(
        generate_clean_sl_callback()
    )

    # set modal value TP to O
    app.callback(Output(pair + "TP", "value"), [Input(pair + "Buy", "n_clicks")])(
        generate_clean_tp_callback()
    )

    # hide modal
    app.callback(
        Output(pair + "Buy", "n_clicks"),
        [
            Input(pair + "closeModal", "n_clicks"),
            Input(pair + "button_order", "n_clicks"),
        ],
    )(generate_modal_close_callback())

    # updates modal figure
    app.callback(
        Output(pair + "modal_graph", "figure"),
        [Input(pair + "index", "children"), Input(pair + "Buy", "n_clicks")],
        [State(pair + "modal_graph", "figure")],
    )(generate_modal_figure_callback(pair))

    # each pair saves its orders in hidden div
    app.callback(
        Output(pair + "orders", "children"),
        [Input(pair + "button_order", "n_clicks")],
        [
            State(pair + "volume", "value"),
            State(pair + "trade_type", "value"),
            State(pair + "SL", "value"),
            State(pair + "TP", "value"),
            State(pair + "orders", "children"),
            State(pair + "ask", "children"),
            State(pair + "bid", "children"),
        ],
    )(generate_order_button_callback(pair))

# updates hidden div with all the clicked charts
app.callback(
    Output("charts_clicked", "children"),
    [Input(pair + "Button_chart", "n_clicks") for pair in currencies],
    [State("charts_clicked", "children")],
)(generate_chart_button_callback())

# updates hidden orders div with all pairs orders
app.callback(
    Output("orders", "children"),
    [Input(pair + "orders", "children") for pair in currencies]
    + [Input(pair + "bid", "children") for pair in currencies]
    + [Input(pair + "ask", "children") for pair in currencies]
    + [Input("closable_orders", "value")],
    [State("orders", "children")],
)(generate_update_orders_div_callback())

# Orders Table
@app.callback(
    Output("orders_table", "children"),
    [Input("orders", "children"), Input("dropdown_positions", "value")],
)
def update_order_table(orders, position):
    headers = [
        "Order Id",
        "Time",
        "Type",
        "Volume",
        "Symbol",
        "TP",
        "SL",
        "Price",
        "Profit",
        "Status",
    ]

    # If tab is closed
    if position == "closed":
        headers += ["Close Time", "Close Price"]

    # If there are no orders
    if orders is None or orders is "[]":
        return html.Tr(children=[html.Th(title) for title in headers])

    rows = []
    list_order = json.loads(orders)
    for order in list_order:
        tr_childs = []
        for attr in order:
            if str(order["status"]) == position:
                tr_childs.append(html.Td(order[attr]))
        # Color row based on profitability of order
        if float(order["profit"]) >= 0:
            rows.append(html.Tr(className="profit", children=tr_childs))
        else:
            rows.append(html.Tr(className="no-profit", children=tr_childs))

    return [html.Tr([html.Th(title) for title in headers])] + rows


# Update Options in dropdown for Open and Close positions
@app.callback(Output("dropdown_positions", "options"), [Input("orders", "children")])
def update_positions_dropdown(orders):
    closeOrders = 0
    openOrders = 0
    if orders is not None:
        orders = json.loads(orders)
        for order in orders:
            if order["status"] == "closed":
                closeOrders += 1
            if order["status"] == "open":
                openOrders += 1
    return [
        {"label": "Open positions (" + str(openOrders) + ")", "value": "open"},
        {"label": "Closed positions (" + str(closeOrders) + ")", "value": "closed"},
    ]


# Updates close orders dropdown options
@app.callback(Output("closable_orders", "options"), [Input("orders", "children")])
def update_close_dropdown(orders):
    options = []
    if orders is not None:
        orders = json.loads(orders)
        for order in orders:
            if order["status"] == "open":
                options.append({"label": order["id"], "value": order["id"]})
    return options


# Callback to update Top Bar values
@app.callback(Output("top_bar", "children"), [Input("orders", "children")])
def update_top_bar(orders):
    if orders is None or orders is "[]":
        return get_top_bar()

    orders = json.loads(orders)
    open_pl = 0
    balance = 50000
    free_margin = 50000
    margin = 0

    for order in orders:
        if order["status"] == "open":
            open_pl += float(order["profit"])
            conversion_price = (
                1 if order["symbol"][:3] == "USD" else float(order["price"])
            )
            margin += (float(order["volume"]) * 100000) / (200 * conversion_price)
        else:
            balance += float(order["profit"])

    equity = balance - open_pl
    free_margin = equity - margin
    margin_level = "%" if margin == 0 else "%2.F" % ((equity / margin) * 100) + "%"
    equity = "%.2F" % equity
    balance = "%.2F" % balance
    open_pl = "%.2F" % open_pl
    free_margin = "%.2F" % free_margin
    margin = "%2.F" % margin

    return get_top_bar(balance, equity, margin, free_margin, margin_level, open_pl)


@app.callback(Output("live_clock", "children"), [Input("interval", "n_intervals")])
def update_time(n):
    return datetime.datetime.now().strftime("%H:%M:%S")


@app.callback(Output("news", "children"), [Input("i_news", "n_intervals")])
def update_news_div(n):
    return update_news()


if __name__ == "__main__":
    app.run_server(debug=True, threaded=True)
