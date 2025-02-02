from dash import Dash, html, dcc, Input, Output, State, callback, callback_context
import dash_bootstrap_components as dbc
import dash_daq as daq
from utils.helper_functions import df
from datetime import datetime, date
from constants import gauge_size


def header(
    app, header_color, header, subheader=None, header_background_color="transparent"
):
    left_headers = html.Div(
        [
            html.Div(header, className="header-title"),
            html.Div(subheader, className="subheader-title"),
        ],
        style={"color": header_color},
    )

    logo = html.Img(src=app.get_asset_url("images/plotly-logo-light-theme.png"))
    logo_link = html.A(logo, href="https://plotly.com/get-demo/", target="_blank")
    demo_link = html.A(
        "LEARN MORE",
        href="https://plotly.com/dash/",
        target="_blank",
        className="demo-button",
    )
    right_logos = html.Div([demo_link, logo_link], className="header-logos")

    return html.Div(
        [left_headers, right_logos],
        className="header",
        style={"background-color": header_background_color},
    )


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

    logo_image = html.Img(src=app.get_asset_url("dash-logo.png"), style={"height": 50})

    link_btns = html.Div(
        style={"float": "right"},
        children=[
            html.A(
                dbc.Button(
                    "Enterprise Demo",
                    color="primary",
                    className="mr-1",
                ),
                href="https://plotly.com/get-demo/",
                target="_blank",
            ),
            html.A(
                dbc.Button("Source Code", color="secondary", className="mr-1"),
                href="https://github.com/plotly/dash-sample-apps/tree/main/apps/dash-turbine-maintenance",
                target="_blank",
            ),
            html.A(
                logo_image,
                href="https://plotly.com/dash/",
                style={"margin-left": "15px"},
            ),
        ],
    )

    return dbc.Row(
        [dbc.Col([dbc.Row([title]), dbc.Row([info_about_app])]), dbc.Col([link_btns])]
    )


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
