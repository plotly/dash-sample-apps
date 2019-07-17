import dash
from dash.dependencies import Input, Output, State
import dash_daq as daq
from dash_daq import DarkThemeProvider
import dash_html_components as html
import dash_core_components as dcc
from dash.exceptions import PreventUpdate

import numpy as np
from scipy import signal

app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)
app.config["suppress_callback_exceptions"] = True
server = app.server

font_color = {"dark": "#ffffff", "light": "#222"}
background_color = {"dark": "#2a3f5f", "light": "#ffffff"}
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
        className="knobs"
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
        className="led-displays"
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
                labelStyle={"display": "inline-block"}
            )
        ]
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
                        style={"margin-bottom": "15px"},
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
                        style={"margin-bottom": "15px"},
                    ),
                ],
                style={"margin": "15px 0"},
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
            )
        ],
    )


app.layout = html.Div(
    id="main-page",
    children=[

        # Header
        html.Div(
            id="header",
            className="banner row",
            children=[
                # Toggle
                html.Div(
                    className="row toggleDiv",
                    children=[
                        daq.ToggleSwitch(
                            id="toggleTheme",
                            size=30,
                            value=False,
                        )
                    ],
                ),
                # Title and Image
                html.Div(
                    className="row",
                    children=[
                        html.H2("Dash DAQ: Function Generator & Oscilloscope Control Panel"),
                        html.Img(
                            src="https://s3-us-west-1.amazonaws.com/plotly-tutorials/"
                            + "excel/dash-daq/dash-daq-logo-by-plotly-stripe+copy.png",
                            className="logo"
                        ),
                        html.H6("Dash DAQ: Function Generator & Oscilloscope Control Panel"),
                    ]
                )
            ],
        ),

        html.Div(
            className="row",
            children=[
                    # LEFT PANEL - SETTINGS
                    html.Div(
                        className="five columns left-panel",
                        children=[
                            html.Div(
                                id="dark-theme-components",
                                children=DarkThemeProvider(
                                    theme=theme,
                                    children=[
                                        power_setting_div(None, "1"),
                                        function_setting_div(None, "1"),
                                    ],
                                ),
                            ),
                            daq.ColorPicker(
                                id="color-picker",
                                label="Color Picker",
                                value=dict(hex="#6682C0"),
                                size=164,
                                style={
                                    "marginTop": "20px",
                                    "backgroundColor": "inherit",
                                },
                            ),
                        ],
                    ),

                    # RIGHT PANEL - OSCILLATIONS
                    html.Div(
                        className="seven columns right-panel",
                        children=[
                            html.Div(
                                [html.H3("Graph", id="graph-title")],
                                style={"color": theme["primary"]},
                                className="Title",
                            ),
                            dcc.Tabs(
                                id="tabs",
                                children=[dcc.Tab(label="Run #1", value="1")],
                                value="1",
                                className="oscillator-tabs",
                                colors={
                                    "border": "#d6d6d6",
                                    "primary": "#6682C0",
                                    "background": "#f2f2f2",
                                },
                            ),
                            html.Div(
                                className="row oscope-info",
                                children=[
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Div(
                                                        id="graph-info",
                                                        children="-",
                                                        style={
                                                            "border": "1px solid"
                                                            + theme["primary"]
                                                        },
                                                    )
                                                ],
                                                className="row graph-param",
                                            )
                                        ],
                                        className="six columns",
                                    ),
                                    html.Button(
                                        "+",
                                        id="new-tab",
                                        n_clicks=0,
                                        type="submit",
                                        style={
                                            "height": "20px",
                                            "width": "20px",
                                            "padding": "2px",
                                            "lineHeight": "10px",
                                            "float": "right",
                                            "color": "inherit",
                                        },
                                    ),
                                ],
                            ),
                            html.Hr(),
                            dcc.Graph(id="oscope-graph", figure={}),
                        ],
                    ),
                ]
            
        ),
        dcc.Store(id="control-inputs", data={}),  # {tabs_number: {value1:x, value2:x}
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
            xaxis=dict(title="s", color=axis, titlefont=dict(family="Dosis", size=13)),
            yaxis=dict(
                title="Voltage (mV)",
                color=axis,
                range=[-10, 10],
                titlefont=dict(family="Dosis", size=13),
            ),
            margin={"l": 40, "b": 40, "t": 20, "r": 50},
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


# Update colors upon color-picker changes
@app.callback(
    [
        Output("power-title", "style"),
        Output("function-title", "style"),
        Output("graph-title", "style"),
        Output("graph-info", "style"),
        Output("tabs", "color"),
        Output("header", "style"),
    ],
    [Input("color-picker", "value")],
)
def color_update(color):
    return (
        list({"color": color["hex"]} for _ in range(3))
        + [{"border": ("1px solid " + color["hex"]), "color": "inherit"}]
        + [{"border": color["hex"]}]
        + [{"backgroundColor": color["hex"]}]
    )


@app.callback(
    Output("new-tab", "style"),
    [Input("toggleTheme", "value")],
    [State("new-tab", "style")],
)
def update_click_btn_color(turn_dark, cur_style):
    if turn_dark:
        cur_style.update(backgroundColor="#EBF0F8")
    raise PreventUpdate


# Callback updating backgrounds
@app.callback(Output("main-page", "style"), [Input("toggleTheme", "value")])
def update_background(turn_dark):
    if turn_dark:
        return {
            "backgroundColor": background_color["dark"],
            "color": font_color["dark"],
        }
    else:
        return {
            "backgroundColor": background_color["light"],
            "color": font_color["light"],
        }


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
            dcc.Tab(label="Run #{}".format(i), value="{}".format(i))
            for i in range(1, len(cur_inputs) + 2)
        ),
        str(len(cur_inputs) + 1),
    )


if __name__ == "__main__":
    app.run_server(debug=True)
