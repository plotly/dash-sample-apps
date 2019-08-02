import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_core_components as dcc
import dash_daq as daq

import plotly.graph_objs as go
from time import sleep
import os
import numpy as np

import fgen_afg3021 as fgen
import osc_tds350 as osc


app = dash.Dash()

app.config["suppress_callback_exceptions"] = True

server = app.server

tabs = [{"label": "Run #{}".format(i), "value": i} for i in range(1, 2)]

tab = 1

runs = {}

fgen.open_port()


def knobs(cur_input, cur_tab):
    return html.Div(
        [
            daq.Knob(
                value=cur_input[cur_tab]["frequency_input"],
                id="frequency-input",
                label="Frequency (Hz)",
                labelPosition="bottom",
                size=70,
                color=theme["primary"],
                max=2500000,
                min=1e5,
                style={"background": "transparent"},
                className="four columns",
            ),
            daq.Knob(
                value=cur_input[cur_tab]["amplitude_input"],
                id="amplitude-input",
                label="Amplitude (mV)",
                labelPosition="bottom",
                size=70,
                scale={"labelInterval": 10},
                color=theme["primary"],
                max=10,
                min=0,
                className="four columns",
            ),
            daq.Knob(
                value=cur_input[cur_tab]["offset_input"],
                id="offset-input",
                label="Offset (mV)",
                labelPosition="bottom",
                size=72,
                scale={"labelInterval": 10},
                color=theme["primary"],
                max=10,
                min=0,
                className="four columns",
            ),
        ],
        className="knobs",
    )


def led_displays(cur_input, cur_tab):
    return html.Div(
        [
            daq.LEDDisplay(
                id="frequency-display",
                size=10,
                value=cur_input[cur_tab]["frequency_input"],
                label="Frequency (Hz)",
                labelPosition="bottom",
                color=theme["primary"],
                className="four columns",
            ),
            daq.LEDDisplay(
                id="amplitude-display",
                size=10,
                value=cur_input[cur_tab]["amplitude_input"],
                label="Amplitude (mV)",
                labelPosition="bottom",
                color=theme["primary"],
                className="four columns",
            ),
            daq.LEDDisplay(
                id="offset-display",
                size=10,
                value=cur_input[cur_tab]["amplitude_input"],
                label="Offset (mV)",
                labelPosition="bottom",
                color=theme["primary"],
                className="four columns",
            ),
        ],
        className="led-displays",
    )


def radioitem(cur_input, cur_tab):
    return html.Div(
        className="radio-items",
        children=[
            dcc.RadioItems(
                id="function-type",
                options=[
                    {"label": "Sine", "value": "SIN"},
                    {"label": "Square", "value": "SQUARE"},
                    {"label": "Ramp", "value": "RAMP"},
                ],
                value=cur_input[cur_tab]["function_type"],
                labelStyle={"display": "inline-block"},
            )
        ],
    )


def power_setting_div(cur_inputs, cur_tab):
    if cur_inputs is None or len(cur_inputs) == 0:
        cur_inputs = init_input
    return html.Div(
        className="row power-settings-tab",
        children=[
            # Title
            html.Div(
                className="Title",
                children=html.H3(
                    "Power", id="power-title", style={"color": theme["primary"]}
                ),
            ),
            # Power Controllers
            html.Div(
                [
                    html.Div(
                        [
                            daq.PowerButton(
                                id="function-generator",
                                on=cur_inputs[cur_tab]["function_generator"],
                                label="Function Generator",
                                labelPosition="bottom",
                                color=theme["primary"],
                            )
                        ],
                        className="six columns",
                    ),
                    html.Div(
                        [
                            daq.PowerButton(
                                id="oscilloscope",
                                on=cur_inputs[cur_tab]["oscilloscope"],
                                label="Oscilloscope",
                                labelPosition="bottom",
                                color=theme["primary"],
                            )
                        ],
                        className="six columns",
                    ),
                ]
            ),
        ],
    )


def function_setting_div(cur_input, cur_tab):
    if cur_input is None or len(cur_input) == 0:
        cur_input = init_input
    return html.Div(
        className="row power-settings-tab",
        children=[
            html.Div(
                className="Title",
                style={"color": theme["primary"]},
                children=html.H3("Function", id="function-title"),
            ),
            html.Div(
                children=[
                    # Knobs
                    knobs(cur_input, cur_tab),
                    # LED Displays
                    led_displays(cur_input, cur_tab),
                    # # RadioItems
                    radioitem(cur_input, cur_tab),
                ]
            ),
        ],
    )


# Main App
app.layout = html.Div(
    id="main-page",
    children=[
        # Header
        html.Div(
            id="header",
            className="banner row",
            children=[
                # Logo and Title
                html.Div(
                    className="banner-logo-and-title",
                    children=[
                        html.Img(
                            src=app.get_asset_url("dash-logo-white.png"),
                            className="logo",
                        ),
                        html.H2(
                            "Dash DAQ: Function Generator & Oscilloscope Control Panel"
                        ),
                    ],
                ),
                # Toggle
                html.Div(
                    className="row toggleDiv",
                    children=daq.ToggleSwitch(
                        id="toggleTheme", size=40, value=False, color="#00418e"
                    ),
                ),
            ],
        ),
        html.Div(
            className="row",
            children=[
                # LEFT PANEL - SETTINGS
                html.Div(
                    className="five columns",
                    children=[
                        html.Div(
                            id="left-panel",
                            children=[
                                html.Div(
                                    id="dark-theme-components",
                                    className="left-panel-controls",
                                    children=DarkThemeProvider(
                                        theme=theme,
                                        children=[
                                            power_setting_div(None, "1"),
                                            function_setting_div(None, "1"),
                                        ],
                                    ),
                                ),
                                html.Div(
                                    className="left-panel-color-picker",
                                    children=[
                                        daq.ColorPicker(
                                            id="color-picker",
                                            label=" ",
                                            value=dict(hex="#0054A6"),
                                            size=164,
                                        )
                                    ],
                                ),
                            ],
                        )
                    ],
                ),
                # RIGHT PANEL - OSCILLATIONS
                html.Div(
                    className="seven columns",
                    children=[
                        html.Div(
                            id="right-panel",
                            className="right-panel",
                            children=[
                                html.Div(
                                    id="card-right-panel-info",
                                    className="light-card",
                                    children=[
                                        html.Div(
                                            [html.H3("Graph", id="graph-title")],
                                            style={"color": theme["primary"]},
                                            className="Title",
                                        ),
                                        dcc.Tabs(
                                            id="tabs",
                                            children=[
                                                dcc.Tab(label="Run #1", value="1")
                                            ],
                                            value="1",
                                            className="oscillator-tabs",
                                        ),
                                        html.Div(
                                            className="row oscope-info",
                                            children=[
                                                html.Div(
                                                    id="div-graph-info",
                                                    className="graph-param",
                                                    children=html.Div(
                                                        id="graph-info", children="-"
                                                    ),
                                                ),
                                                html.Button(
                                                    "+",
                                                    id="new-tab",
                                                    n_clicks=0,
                                                    type="submit",
                                                ),
                                            ],
                                        ),
                                    ],
                                ),
                                html.Div(
                                    id="card-graph",
                                    className="light-card",
                                    children=dcc.Graph(id="oscope-graph", figure={}),
                                ),
                            ],
                        )
                    ],
                ),
            ],
        ),
        dcc.Store(id="control-inputs", data={}),
        dcc.Interval(id="update-oscope", interval=2000, n_intervals=0),
    ],
)

# Callback to update theme layout
@app.callback(
    Output("dark-theme-components", "children"),
    [Input("toggleTheme", "value"), Input("color-picker", "value")],
    [State("control-inputs", "data"), State("tabs", "value")],
)
def turn_dark(turn_dark, color_pick, cur_inputs, cur_tab_value):
    theme.update(dark=turn_dark)

    if color_pick is not None:
        theme.update(primary=color_pick["hex"])

    return DarkThemeProvider(
        theme=theme,
        children=[
            power_setting_div(cur_inputs, cur_tab_value),
            function_setting_div(cur_inputs, cur_tab_value),
        ],
    )


# Callback updating backgrounds
@app.callback(
    [
        Output("main-page", "className"),
        Output("left-panel", "className"),
        Output("card-right-panel-info", "className"),
        Output("card-graph", "className"),
    ],
    [Input("toggleTheme", "value")],
)
def update_background(turn_dark):

    if turn_dark:
        return ["dark-main-page", "dark-card", "dark-card", "dark-card"]
    else:
        return ["light-main-page", "light-card", "light-card", "light-card"]


@app.callback(Output("graph-info", "style"), [Input("color-picker", "value")])
def color_info(color):
    return {"textAlign": "center", "border": "2px solid " + color["hex"]}


# Update colors upon color-picker changes
@app.callback(
    [
        Output("power-title", "style"),
        Output("function-title", "style"),
        Output("graph-title", "style"),
        Output("div-graph-info", "style"),
        Output("tabs", "color"),
    ],
    [Input("color-picker", "value")],
)
def color_update(color):
    return (
        list({"color": color["hex"]} for _ in range(3))
        + [{"border": ("1px solid " + color["hex"]), "color": "inherit"}]
        + [{"border": color["hex"]}]
    )


# Callbacks for knob inputs
@app.callback(Output("frequency-display", "value"), [Input("frequency-input", "value")])
def update_frequency_display(value):
    fgen.set_frequency(value)
    return value


@app.callback(Output("amplitude-display", "value"), [Input("amplitude-input", "value")])
def update_amplitude_display(value):
    fgen.set_amplitude(value)
    return value


@app.callback(Output("offset-display", "value"), [Input("offset-input", "value")])
def update_offset_display(value):
    fgen.set_offset(value)
    return value


@app.callback(Output("function-type", "value"), [Input("function_type", "value")])
def update_fgen_wave(value):
    fgen.set_wave(value)
    return value


# Callbacks graph and graph info
@app.callback(
    Output("graph-info", "children"),
    [Input("oscope-graph", "figure"), Input("tabs", "value")],
)
def update_info(_, value):
    if "" + str(value) in runs:
        return runs["" + str(value)][1]
    return "-"


@app.callback(
    Output("oscope-graph", "figure"),
    [Input("update-oscope", "n_intervals"), Input("tabs", "value")],
)
def update_output(_, value):
    global tab
    time = np.linspace(-0.000045, 0.000045, 1e3)
    zero = dict(
        data=[dict(x=time, y=[0] * len(time), marker={"color": "#2a3f5f"})],
        layout=go.Layout(
            xaxis={
                "title": "s",
                "color": "#506784",
                "titlefont": dict(family="Dosis", size=15),
            },
            yaxis={
                "title": "Voltage (mV)",
                "color": "#506784",
                "titlefont": dict(family="Dosis", size=15),
            },
            margin={"l": 40, "b": 40, "t": 0, "r": 50},
            plot_bgcolor="#F3F6FA",
        ),
    )

    if tab is not value:
        if "" + str(value) in runs:
            tab = value
            return runs["" + str(value)][0]
        tab = value
        return zero

    else:
        figure = {
            "data": osc.get_data(),
            "layout": go.Layout(
                xaxis={
                    "title": "s",
                    "color": "#506784",
                    "titlefont": dict(family="Dosis", size=15),
                },
                yaxis={
                    "title": "Voltage (mV)",
                    "color": "#506784",
                    "titlefont": dict(family="Dosis", size=15),
                    "autorange": False,
                    "range": [-10, 10],
                },
                margin={"l": 40, "b": 40, "t": 0, "r": 50},
                plot_bgcolor="#F3F6FA",
            ),
        }

        runs["" + str(value)] = (
            figure,
            str(fgen.get_wave())
            + " | "
            + str(fgen.get_frequency())
            + "Hz"
            + " | "
            + str(fgen.get_amplitude())
            + "mV"
            + " | "
            + str(fgen.get_offset())
            + "mV",
        )

        # wait to update the runs variable
        sleep(0.10)

        return figure


# Callback for adding tabs
@app.callback(
    [Output("tabs", "children"), Output("tabs", "value")],
    [Input("new-tab", "n_clicks")],
    [State("control-inputs", "data"), State("color-picker", "value")],
)
def update_total_tab_number(n_clicks, cur_inputs, color):
    return (
        list(
            dcc.Tab(
                label="Run #{}".format(i),
                value="{}".format(i),
                selected_style={
                    "color": color["hex"],
                    "backgroundColor": "transparent",
                    "border": "none",
                },
                style={
                    "color": "black",
                    "backgroundColor": "transparent",
                    "border": "none",
                },
            )
            for i in range(1, len(cur_inputs) + 2)
        ),
        str(len(cur_inputs) + 1),
    )


if __name__ == "__main__":
    app.run_server(port=8000, debug=True)
