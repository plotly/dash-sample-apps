import pandas as pd
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq

import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
from datetime import datetime, date
from data_preprocessing import data_preprocessing
import pickle

app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    external_stylesheets=[dbc.themes.SUPERHERO],
)
server = app.server
app.title = "Predictive Maintenance Dashboard"


def logo(app):
    title = html.H5(
        "PREDICTIVE MAINTENANCE DASHBOARD FOR WIND TURBINES",
        style={"marginTop": 5, "marginLeft": "10px"},
    )

    info_about_app = html.H6(
        "This Dashboard is focused on estimating the Remaining Useful Life (RUL) in wind turbines. RUL is defined "
        " as the time until the next fault.",
        style={"marginLeft": "10px"},
    )

    logo_image = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 50}
    )
    link = html.A(logo_image, href="https://plotly.com/dash/")

    return dbc.Row(
        [dbc.Col([dbc.Row([title]), dbc.Row([info_about_app])]), dbc.Col(link)]
    )


df, df_button, x_test, y_test = data_preprocessing()

predict_button = dbc.Card(
    className="mt-auto",
    children=[
        dbc.CardBody(
            [
                html.Div(
                    [
                        dbc.Button(
                            "Predict",
                            id="predict-button",
                            color="#fec036",
                            size="lg",
                            style={"color": "#fec036"},
                        ),
                    ]
                )
            ],
            style={
                "text-align": "center",
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
                "border-left": "1px solid rgb(216, 216, 216)",
                "border-right": "1px solid rgb(216, 216, 216)",
                "border-bottom": "1px solid rgb(216, 216, 216)",
            },
        )
    ],
)

get_new_information_button = dbc.Card(
    className="mt-auto",
    children=[
        dbc.CardBody(
            [
                html.Div(
                    [
                        dbc.Button(
                            "Get New Data",
                            id="get-new-info-button",
                            color="#fec036",
                            size="lg",
                            style={"color": "#fec036"},
                        ),
                    ]
                )
            ],
            style={
                "text-align": "center",
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
                "border-left": "1px solid rgb(216, 216, 216)",
                "border-right": "1px solid rgb(216, 216, 216)",
                "border-bottom": "1px solid rgb(216, 216, 216)",
            },
        )
    ],
)


graphs = dbc.Card(
    children=[
        dbc.CardBody(
            [
                html.Div(
                    [
                        dcc.Graph(
                            id="Main-Graph",
                            figure={
                                "layout": {
                                    "margin": {"t": 30, "r": 35, "b": 40, "l": 50},
                                    "xaxis": {
                                        "dtick": 5,
                                        "gridcolor": "#636363",
                                        "showline": False,
                                    },
                                    "yaxis": {"showgrid": False, "showline": False},
                                    "plot_bgcolor": "black",
                                    "paper_bgcolor": "black",
                                    "font": {"color": "gray"},
                                },
                            },
                            config={"displayModeBar": False},
                        ),
                        html.Pre(id="update-on-click-data"),
                    ],
                    style={"width": "98%", "display": "inline-block"},
                ),
                html.Div(
                    [
                        dcc.Dropdown(
                            id="feature-dropdown",
                            options=[
                                {"label": label, "value": label} for label in df.columns
                            ],
                            value="",
                            multi=False,
                            searchable=False,
                        )
                    ],
                    style={
                        "width": "33%",
                        "display": "inline-block",
                        "color": "black",
                    },
                ),
                html.Div(
                    [
                        dcc.DatePickerRange(
                            id="date-picker",
                            min_date_allowed=date(2014, 5, 1),  # need to change this
                            max_date_allowed=date(2015, 4, 9),
                            initial_visible_month=date(2014, 5, 1),
                            start_date_placeholder_text="Start Period",
                            end_date_placeholder_text="End Period",
                            calendar_orientation="vertical",
                        ),
                        html.Div(id="output-container-date-picker-range"),
                    ],
                    style={
                        "vertical-align": "top",
                        "position": "absolute",
                        "right": "3%",
                        "float": "right",
                        "display": "inline-block",
                        "color": "black",
                    },
                ),
            ],
            style={
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        )
    ]
)

rul_estimation_indicator = dbc.Card(
    children=[
        dbc.CardHeader(
            "System RUL Estimation (days)",
            style={
                "text-align": "center",
                "color": "white",
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
        dbc.CardBody(
            [
                daq.LEDDisplay(
                    id="rul-estimation-indicator-led",
                    size=24,
                    color="#fec036",
                    style={"color": "#black"},
                    backgroundColor="#2b2b2b",
                    value="0.0",
                )
            ],
            style={
                "text-align": "center",
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
    ]
)

info_box = dbc.Card(
    children=[
        dbc.CardBody(
            [
                html.Div(
                    dcc.Textarea(
                        id="Info-Textbox",
                        placeholder="This field is used to display information about a feature displayed "
                        "on the graph and estimated RUL. In order to estimate the RUL, use "
                        "the button 'Get New Data' and then, 'Predict'. The estimated RUL will be "
                        "printed.",
                        rows=8,
                        style={
                            "width": "100%",
                            "height": "100%",
                            "background-color": "black",
                            "color": "#fec036",
                            "placeholder": "#fec036",
                            "fontFamily": "Arial",
                            "fontSize": "16",
                            "display": "inline-block",
                        },
                    )
                )
            ],
            style={
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
    ],
)

blade_angle_display = dbc.Card(
    children=[
        dbc.CardHeader(
            "Blade Angle",
            style={
                "text-align": "center",
                "color": "white",
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
        dbc.CardBody(
            [
                html.Div(
                    daq.Gauge(
                        id="blade-angle-information-gauge",
                        min=min(df["WEC: ava. blade angle A"]),
                        max=max(
                            df["WEC: ava. blade angle A"]
                        ),  # This one should be the theoretical maximum
                        value=0,
                        showCurrentValue=True,
                        color="#fec036",
                        style={
                            "align": "center",
                            "display": "flex",
                            "marginTop": "5%",
                            "marginBottom": "-10%",
                        },
                    ),
                    className="m-auto",
                    style={
                        "display": "flex",
                        "backgroundColor": "black",
                        "border-radius": "1px",
                        "border-width": "5px",
                    },
                )
            ],
            className="d-flex",
            style={
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
    ],
    style={"height": "95%"},
)

active_power_display = dbc.Card(
    children=[
        dbc.CardHeader(
            "Active Power [kW]",
            style={
                "text-align": "center",
                "color": "white",
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
        dbc.CardBody(
            [
                html.Div(
                    daq.Gauge(
                        id="active-power-information-gauge",
                        min=min(df["WEC: ava. Power"]),
                        max=max(
                            df["WEC: ava. Power"]
                        ),  # This one should be the theoretical maximum
                        value=100,
                        showCurrentValue=True,
                        color="#fec036",
                        style={
                            "align": "center",
                            "display": "flex",
                            "marginTop": "5%",
                            "marginBottom": "-10%",
                        },
                    ),
                    className="m-auto",
                    style={
                        "display": "flex",
                        "backgroundColor": "black",
                        "border-radius": "1px",
                        "border-width": "5px",
                    },
                )
            ],
            className="d-flex",
            style={
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
    ],
    style={"height": "95%"},
)

active_power_from_wind_display = dbc.Card(
    children=[
        dbc.CardHeader(
            "Active Power Available from Wind [kW]",
            style={
                "display": "inline-block",
                "text-align": "center",
                "color": "white",
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
        dbc.CardBody(
            [
                html.Div(
                    daq.Gauge(
                        id="active-power-from-wind-information-gauge",
                        min=min(df["WEC: ava. available P from wind"]),
                        max=max(df["WEC: ava. available P from wind"]),
                        value=10,
                        showCurrentValue=True,
                        color="#fec036",
                        style={
                            "align": "center",
                            "display": "flex",
                            "marginTop": "5%",
                            "marginBottom": "-10%",
                        },
                    ),
                    className="m-auto",
                    style={
                        "display": "flex",
                        "backgroundColor": "black",
                        "border-radius": "1px",
                        "border-width": "5px",
                    },
                )
            ],
            className="d-flex",
            style={
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
    ],
    style={"height": "95%"},
)

wind_speed_information = dbc.Card(
    className="mt-auto",
    children=[
        dbc.CardHeader(
            "Wind Speed [m/s]",
            style={
                "text-align": "center",
                "color": "white",
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
        dbc.CardBody(
            [
                html.Div(
                    daq.Gauge(
                        id="wind-power-information-gauge",
                        min=min(df["WEC: ava. windspeed"]),
                        max=int(max(df["WEC: ava. windspeed"])),
                        value=0,
                        showCurrentValue=True,
                        color="#fec036",
                        style={
                            "align": "center",
                            "display": "flex",
                            "marginTop": "5%",
                            "marginBottom": "-10%",
                        },
                    ),
                    className="m-auto",
                    style={
                        "display": "flex",
                        "backgroundColor": "black",
                        "border-radius": "1px",
                        "border-width": "5px",
                    },
                )
            ],
            className="d-flex",
            style={
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
    ],
    style={"height": "95%"},
)

reactive_power_display = dbc.Card(
    className="mt-auto",
    children=[
        dbc.CardHeader(
            "Reactive Power [kVAR]",
            style={
                "text-align": "center",
                "color": "white",
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
        dbc.CardBody(
            [
                html.Div(
                    daq.Gauge(
                        id="reactive-power-information-gauge",
                        min=min(df["WEC: ava. reactive Power"]),
                        max=max(df["WEC: ava. reactive Power"]),
                        value=0,
                        showCurrentValue=True,
                        color="#fec036",
                        style={
                            "align": "center",
                            "display": "flex",
                            "marginTop": "5%",
                            "marginBottom": "-10%",
                        },
                    ),
                    className="m-auto",
                    style={
                        "display": "flex",
                        "backgroundColor": "black",
                        "border-radius": "1px",
                        "border-width": "5px",
                    },
                )
            ],
            className="d-flex",
            style={
                "backgroundColor": "black",
                "border-radius": "1px",
                "border-width": "5px",
                "border-top": "1px solid rgb(216, 216, 216)",
            },
        ),
    ],
    style={"height": "95%"},
)

gauge_size = "auto"
app.layout = dbc.Container(
    fluid=True,
    children=[
        logo(app),
        dbc.Row(
            [
                dbc.Col(graphs, xs=10, md=10, lg=10, width=10),
                dbc.Col(
                    [
                        dbc.Row(
                            dbc.Col(
                                rul_estimation_indicator, xs=12, md=12, lg=12, width=12
                            )
                        ),
                        dbc.Row(dbc.Col(info_box, xs=12, md=12, lg=12, width=12)),
                        dbc.Row(
                            dbc.Col(
                                get_new_information_button,
                                xs=12,
                                md=12,
                                lg=12,
                                width=12,
                            )
                        ),
                        dbc.Row(dbc.Col(predict_button, xs=12, md=12, lg=12, width=12)),
                    ]
                ),
            ],
            justify="start",
            style={"display": "flex", "marginBottom": "-3%"},
        ),
        dbc.Row(
            [
                dbc.Col(
                    active_power_display,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
                dbc.Col(
                    active_power_from_wind_display,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
                dbc.Col(
                    reactive_power_display,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
                dbc.Col(
                    wind_speed_information,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
                dbc.Col(
                    blade_angle_display,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
            ],
            style={"marginTop": "3%"},
        ),
    ],
)


def fig_update_layout(fig):
    fig.update_layout(
        xaxis=dict(
            showline=False,
            showgrid=False,
            showticklabels=True,
            zeroline=False,
            gridcolor="#636363",
            linecolor="rgb(204, 204, 204)",
            linewidth=2,
            tickfont=dict(family="Arial", size=12, color="white",),
            title=dict(font=dict(family="Arial", size=24, color="#fec036"),),
        ),
        yaxis=dict(
            showline=False,
            showgrid=False,
            showticklabels=True,
            zeroline=False,
            gridcolor="#636363",
            linecolor="rgb(204, 204, 204)",
            linewidth=2,
            tickfont=dict(family="Arial", size=12, color="white",),
            title=dict(font=dict(family="Arial", size=24, color="#fec036"),),
        ),
        autosize=True,
        margin=dict(autoexpand=True, l=50, b=40, r=35, t=30),
        showlegend=False,
        paper_bgcolor="black",
        plot_bgcolor="black",
        title=dict(
            font=dict(family="Arial", size=32, color="darkgray"),
            xanchor="center",
            yanchor="top",
            y=1,
            x=0.5,
        ),
    )
    return fig


@app.callback(
    [
        Output("Main-Graph", "figure"),
        Output("rul-estimation-indicator-led", "value"),
        Output("Info-Textbox", "value"),
    ],
    [
        Input("feature-dropdown", "value"),
        Input("date-picker", "start_date"),
        Input("date-picker", "end_date"),
        Input("get-new-info-button", "n_clicks"),
        Input("predict-button", "n_clicks"),
    ],
)
def update_graph(selected_column, start_date, end_date, n_get_new_info, n_pred):
    if n_pred is None:  # here is my work before prediction button is activated.
        value_rul = 0.0
        information_update = (
            "This field is used to display information about a feature displayed "
            "on the graph and estimated RUL. In order to estimate the RUL, use "
            "the button 'Get New Data' and then, 'Predict'. The estimated RUL will be "
            "printed."
        )
        if n_get_new_info is None:
            if selected_column in list(df):
                if start_date and end_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    end_date_object = datetime.strptime(end_date, "%Y-%m-%d")
                    mask = (df.index > start_date_object) & (
                        df.index <= end_date_object
                    )
                    df_within_dates = df.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                elif start_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    mask = df.index > start_date_object
                    df_within_dates = df.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                else:
                    fig = go.Figure(
                        data=[go.Scatter(x=df.index, y=df[selected_column])]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
            else:
                fig = go.Figure()
                fig = fig_update_layout(fig)
                return fig, value_rul, information_update
        else:
            if selected_column in list(df_button):
                if start_date and end_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    end_date_object = datetime.strptime(end_date, "%Y-%m-%d")
                    mask = (df_button.index > start_date_object) & (
                        df_button.index <= end_date_object
                    )
                    df_within_dates = df_button.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    _information_update = (
                        "New information is received for the last week and covers periods from "
                        + str(df_button.index[0])
                        + " to "
                        + str(df_button.index[-1])
                        + ". To predict"
                        " RUL, use 'Predict' button. To view data for the aforementioned period, choose"
                        " appropriate dates on the calendar."
                    )
                    return fig, value_rul, _information_update
                elif start_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    mask = df_button.index > start_date_object
                    df_within_dates = df_button.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    _information_update = (
                        "New information is received for the last week and covers periods from "
                        + str(df_button.index[0])
                        + " to "
                        + str(df_button.index[-1])
                        + ". To predict"
                        " RUL, use 'Predict' button. To view data for the aforementioned period, choose"
                        " appropriate dates on the calendar."
                    )
                    return fig, value_rul, _information_update
                else:
                    fig = go.Figure(
                        data=[
                            go.Scatter(x=df_button.index, y=df_button[selected_column])
                        ]
                    )
                    fig = fig_update_layout(fig)
                    _information_update = (
                        "New information is received for the last week and covers periods from "
                        + str(df_button.index[0])
                        + " to "
                        + str(df_button.index[-1])
                        + ". To predict"
                        " RUL, use 'Predict' button. To view data for the aforementioned period, choose"
                        " appropriate dates on the calendar."
                    )
                    return fig, value_rul, _information_update
            else:
                fig = go.Figure()
                fig = fig_update_layout(fig)
                _information_update = (
                    "New information is received for the last week and covers periods from "
                    + str(df_button.index[0])
                    + " to "
                    + str(df_button.index[-1])
                    + ". To predict"
                    " RUL, use 'Predict' button. To view data for the aforementioned period, choose"
                    " appropriate dates on the calendar."
                )
                return fig, value_rul, _information_update
    else:
        if n_get_new_info is None:
            value_rul = 0.0
            information_update = "To predict RUL, please use 'Get New Data' button."
            if selected_column in list(df):
                if start_date and end_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    end_date_object = datetime.strptime(end_date, "%Y-%m-%d")
                    mask = (df.index > start_date_object) & (
                        df.index <= end_date_object
                    )
                    df_within_dates = df.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                elif start_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    mask = df.index > start_date_object
                    df_within_dates = df.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                else:
                    fig = go.Figure(
                        data=[go.Scatter(x=df.index, y=df[selected_column])]
                    )
                    return fig, value_rul, information_update
            else:
                fig = go.Figure()
                fig = fig_update_layout(fig)
                return fig, value_rul, information_update
        else:
            model = pickle.load(open("assets/xgb_reg.pkl", "rb"))
            y_pred = model.predict(x_test)
            df_out = pd.DataFrame()
            df_out["pred"] = y_pred
            value_rul = round(max(df_out["pred"]))
            information_update = "RUL is estimated based on the readings from the last week: " "from " + str(
                x_test.index[0]
            ) + " to " + str(
                x_test.index[-1]
            )
            if selected_column in list(df_button):
                if start_date and end_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    end_date_object = datetime.strptime(end_date, "%Y-%m-%d")
                    mask = (df_button.index > start_date_object) & (
                        df_button.index <= end_date_object
                    )
                    df_within_dates = df_button.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                elif start_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    mask = df_button.index > start_date_object
                    df_within_dates = df_button.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                else:
                    fig = go.Figure(
                        data=[
                            go.Scatter(x=df_button.index, y=df_button[selected_column])
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
            else:
                fig = go.Figure()
                fig = fig_update_layout(fig)
                return fig, value_rul, information_update


@app.callback(
    [
        Output("active-power-information-gauge", "value"),
        Output("active-power-from-wind-information-gauge", "value"),
        Output("wind-power-information-gauge", "value"),
        Output("reactive-power-information-gauge", "value"),
        Output("blade-angle-information-gauge", "value"),
    ],
    Input("Main-Graph", "clickData"),
)
def display_click_data(clickData):
    if clickData:
        data_time = clickData["points"][0]["x"]
        value_active_power = df["WEC: ava. Power"].loc[df.index == data_time].values[0]
        value_active_power_wind = (
            df["WEC: ava. available P from wind"].loc[df.index == data_time].values[0]
        )
        value_reactive_power = (
            df["WEC: ava. reactive Power"].loc[df.index == data_time].values[0]
        )
        value_wind_speed = (
            df["WEC: ava. windspeed"].loc[df.index == data_time].values[0]
        )
        value_blade_angle = (
            df["WEC: ava. blade angle A"].loc[df.index == data_time].values[0]
        )
        return (
            value_active_power,
            value_active_power_wind,
            value_wind_speed,
            value_reactive_power,
            value_blade_angle,
        )
    else:
        value_active_power = 0
        value_active_power_wind = 0
        value_reactive_power = 0
        value_wind_speed = 0
        value_blade_angle = 0
        return (
            value_active_power,
            value_active_power_wind,
            value_wind_speed,
            value_reactive_power,
            value_blade_angle,
        )


if __name__ == "__main__":
    app.run_server(debug=True, use_reloader=True)
