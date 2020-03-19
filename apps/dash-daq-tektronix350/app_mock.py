import dash
from dash.dependencies import Input, Output, State
import dash_daq as daq
from dash_daq import DarkThemeProvider
import dash_html_components as html
import dash_core_components as dcc

import numpy as np
from scipy import signal

app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)
app.config["suppress_callback_exceptions"] = True
server = app.server

axis_color = {"dark": "#EBF0F8", "light": "#506784"}
marker_color = {"dark": "#f2f5fa", "light": "#2a3f5f"}

theme = {
    "dark": False,
    "primary": "#447EFF",
    "secondary": "#D3D3D3",
    "detail": "#D3D3D3",
}

init_input = {
    "1": {
        "function_generator": True,
        "oscilloscope": True,
        "frequency_input": 1e6,
        "amplitude_input": 1,
        "offset_input": 0,
        "function_type": "SIN",
    }
}


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
                style={"backgroundColor": "transparent"},
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
                size=70,
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
                className="power-controllers",
                children=[
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
                ],
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
                                    style={"height": 705},
                                    children=DarkThemeProvider(
                                        theme=theme,
                                        children=[
                                            power_setting_div(None, "1"),
                                            function_setting_div(None, "1"),
                                        ],
                                    ),
                                )
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
    ],
)


@app.callback(
    [
        Output("oscilloscope", "on"),
        Output("function-generator", "on"),
        Output("frequency-input", "value"),
        Output("amplitude-input", "value"),
        Output("offset-input", "value"),
        Output("function-type", "value"),
    ],
    [Input("tabs", "value")],
    [
        State("control-inputs", "data"),
        State("oscilloscope", "on"),
        State("function-generator", "on"),
    ],
)
def update_controls(tab_index: str, cur_inputs, osci_on, func_gen):
    if tab_index not in cur_inputs:
        return osci_on, func_gen, 1000000, 1, 0, "SIN"

    td = cur_inputs[tab_index]
    return (
        td["oscilloscope"],
        td["function_generator"],
        td["frequency_input"],
        td["amplitude_input"],
        td["offset_input"],
        td["function_type"],
    )


# update control values
@app.callback(
    Output("control-inputs", "data"),
    [
        Input("oscilloscope", "on"),
        Input("function-generator", "on"),
        Input("frequency-input", "value"),
        Input("amplitude-input", "value"),
        Input("offset-input", "value"),
        Input("function-type", "value"),
    ],
    [State("tabs", "value"), State("control-inputs", "data")],
)
def update_control_values(
    osc_on, fnct_on, frequency, amplitude, offset, wave, sel_tab, cur_inputs
):
    cur_inputs.update(
        {
            sel_tab: dict(
                oscilloscope=osc_on,
                function_generator=fnct_on,
                frequency_input=frequency,
                amplitude_input=amplitude,
                offset_input=offset,
                function_type=wave,
            )
        }
    )
    return cur_inputs


# new tab created not saved to store unless control inputs changes
@app.callback(
    [Output("oscope-graph", "figure"), Output("graph-info", "children")],
    [Input("control-inputs", "data"), Input("toggleTheme", "value")],
    [State("tabs", "value")],
)
def generate_graph(cur_inputs, theme_value, tab_index: str):
    theme_select = "dark" if theme_value else "light"
    axis = axis_color[theme_select]
    marker = marker_color[theme_select]
    time = np.linspace(-0.000045, 0.000045, 1000)

    base_figure = dict(
        data=[dict(x=time, y=[0] * len(time), marker={"color": marker})],
        layout=dict(
            xaxis=dict(
                title="s",
                color=axis,
                titlefont=dict(family="Asap", size=13),
                gridcolor="#61626370",
                showgrid=True,
            ),
            yaxis=dict(
                title="Voltage (mV)",
                color=axis,
                range=[-10, 10],
                titlefont=dict(family="Asap", size=13),
                gridcolor="#61626370",
                showgrid=True,
            ),
            margin={"l": 0, "b": 0, "t": 0, "r": 0},
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
        ),
    )
    if tab_index not in cur_inputs:
        return base_figure, "-"

    tab_data = cur_inputs[tab_index]

    if not tab_data["oscilloscope"]:
        base_figure.update(data=[])
        base_figure["layout"]["xaxis"].update(
            showticklabels=False, showline=False, zeroline=False
        )
        base_figure["layout"]["yaxis"].update(
            showticklabels=False, showline=False, zeroline=False
        )
        return base_figure, "-"

    if not tab_data["function_generator"]:
        return base_figure, "-"

    if tab_data["function_type"] == "SIN":
        y = [
            float(tab_data["offset_input"])
            + (
                float(tab_data["amplitude_input"])
                * np.sin(
                    np.radians(2.0 * np.pi * float(tab_data["frequency_input"]) * t)
                )
            )
            for t in time
        ]

    elif tab_data["function_type"] == "SQUARE":
        y = [
            float(tab_data["offset_input"])
            + float(tab_data["amplitude_input"])
            * (signal.square(2.0 * np.pi * float(tab_data["frequency_input"]) / 10 * t))
            for t in time
        ]

    elif tab_data["function_type"] == "RAMP":
        y = float(tab_data["amplitude_input"]) * (
            np.abs(
                signal.sawtooth(
                    2 * np.pi * float(tab_data["frequency_input"]) / 10 * time
                )
            )
        )
        y = float(tab_data["offset_input"]) + 2 * y - float(tab_data["amplitude_input"])
    else:
        return base_figure, "-"

    base_figure["data"][0].update(y=y)

    info = (
        f'{tab_data["function_type"]}|{tab_data["frequency_input"]}Hz|'
        f'{tab_data["amplitude_input"]} mV | {tab_data["offset_input"]} mV'
    )

    return base_figure, info


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


# Callbacks for knob inputs
@app.callback(Output("frequency-display", "value"), [Input("frequency-input", "value")])
def update_frequency_display(value):
    return value


@app.callback(Output("amplitude-display", "value"), [Input("amplitude-input", "value")])
def update_amplitude_display(value):
    return value


@app.callback(Output("offset-display", "value"), [Input("offset-input", "value")])
def update_offset_display(value):
    return value


# Callback for adding tabs
@app.callback(
    [Output("tabs", "children"), Output("tabs", "value")],
    [Input("new-tab", "n_clicks")],
    [State("control-inputs", "data")],
)
def update_total_tab_number(n_clicks, cur_inputs):
    return (
        list(
            dcc.Tab(
                label="Run #{}".format(i),
                value="{}".format(i),
                selected_style={
                    "color": "black",
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
    app.run_server(debug=True)
