# -*- coding: utf-8 -*-

import os
import sys
import numpy
import time
from threading import Lock
from textwrap import dedent

import dash
import dash_html_components as html
import dash_core_components as dcc

import plotly.graph_objs as go

import dash_daq as daq
from dash.dependencies import Input, Output, State

import DashOceanOpticsSpectrometer as doos
from DashOceanOpticsSpectrometer import Control

DEMO = False

# lock for modifying information about spectrometer
spec_lock = Lock()
# lock for communicating with spectrometer
comm_lock = Lock()

# initialize spec
spec = doos.DashOceanOpticsSpectrometer(spec_lock, comm_lock)

# demo or actual
if ("DASH_PATH_ROUTING" in os.environ) or (
    len(sys.argv) == 2 and sys.argv[1] == "demo"
):
    spec = doos.DemoSpectrometer(spec_lock, comm_lock)
    DEMO = True
else:
    spec = doos.PhysicalSpectrometer(spec_lock, comm_lock)

spec.assign_spec()


app = dash.Dash()
server = app.server

colors = {
    "background": "#bbbbbb",
    "primary": "#efefef",
    "secondary": "#efefef",
    "tertiary": "#dfdfdf",
    "grid-colour": "#eeeeee",
    "accent": "#2222ff",
}

controls = []

# integration time, microseconds
int_time = Control(
    "integration-time",
    "int. time (μs)",
    "NumericInput",
    {
        "id": "integration-time-input",
        "max": spec.int_time_max(),
        "min": spec.int_time_min(),
        "size": 150,
        "value": spec.int_time_min(),
    },
)
controls.append(int_time)

# scans to average over
nscans_avg = Control(
    "nscans-to-average",
    "number of scans",
    "NumericInput",
    {"id": "nscans-to-average-input", "max": 100, "min": 1, "size": 150, "value": 1},
)
controls.append(nscans_avg)

# strobe
strobe_enable = Control(
    "continuous-strobe-toggle",
    "strobe",
    "BooleanSwitch",
    {"id": "continuous-strobe-toggle-input", "color": colors["accent"], "on": False},
)
controls.append(strobe_enable)

# strobe period
strobe_period = Control(
    "continuous-strobe-period",
    "strobe pd. (μs)",
    "NumericInput",
    {
        "id": "continuous-strobe-period-input",
        "max": 100,
        "min": 1,
        "size": 150,
        "value": 1,
    },
)
controls.append(strobe_period)

# light sources
light_sources = Control(
    "light-source",
    "light source",
    "Dropdown",
    {
        "id": "light-source-input",
        "options": spec.light_sources(),
        "placeholder": "select light source",
        "value": "l2" if DEMO else "",
    },
)
controls.append(light_sources)


page_layout = [
    html.Div(
        [
            html.Div(
                [
                    html.Img(
                        src=app.get_asset_url("dash-daq-logo.png"), className="logo"
                    ),
                    html.Div(
                        [
                            html.Label("Int. Time (μs)"),
                            dcc.Input(
                                id="integration-timetimes",
                                type="number",
                                max=spec.int_time_max(),
                                min=spec.int_time_min(),
                                size="150",
                                value=spec.int_time_min(),
                                className="control__dropdowns",
                            ),
                        ],
                        className="control",
                    ),
                    html.Div(
                        [
                            html.Label("Number of Scans"),
                            dcc.Input(
                                id="number-of-snansss",
                                type="number",
                                max="100",
                                min="1",
                                size="150",
                                value="1",
                                className="control__dropdowns",
                            ),
                        ],
                        className="control",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    daq.BooleanSwitch(
                                        label="Strobe",
                                        id="my-daq-booleanswitch",
                                        color="#565656",
                                        on=True,
                                    )
                                ]
                            ),
                            html.Div(
                                [
                                    html.Label("Strobe Pd. (μs)"),
                                    dcc.Input(
                                        id="number-of-snansssss",
                                        type="number",
                                        max="100",
                                        min="1",
                                        size="150",
                                        value="1",
                                        className="control__dropdowns",
                                    ),
                                ]
                            ),
                        ],
                        className="strobe",
                    ),
                    html.Div(
                        [
                            html.Label("Light Source"),
                            dcc.Dropdown(
                                id="light-source-inputt",
                                options=spec.light_sources(),
                                placeholder="Select light source",
                                value="l2" if DEMO else "",
                            ),
                        ],
                        className="control",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Label("Light Intensity"),
                                    daq.Knob(
                                        id="light-intensity-knob",
                                        size=110,
                                        color="#565656",
                                        value=0,
                                    ),
                                ],
                                id="light-intensity-knob-container",
                            ),
                        ],
                        className="control",
                    ),
                    html.Div(
                        id="controls",
                        title="All of the spectrometer parameters that can be changed.",
                        children=[ctrl.create_ctrl_div(True) for ctrl in controls],
                    ),
                    html.Div(
                        [
                            html.Div([html.Label("Autoscale Plot"),]),
                            html.Div(
                                [
                                    daq.BooleanSwitch(
                                        id="my-daq-booleanswitchsss",
                                        color="#565656",
                                        on=True,
                                    ),
                                ],
                            ),
                        ],
                        className="control autoscale",
                    ),
                    html.Div(
                        [html.Button("update", id="submit-buttossn")],
                        className="control",
                    ),
                    html.Div(
                        id="submit-status",
                        title="Contains information about the success or failure of your \
                commands.",
                        children=[""],
                    ),
                ],
                className="one-third column left__section",
            ),
            html.Div(
                [
                    html.Div(
                        id="graph-container",
                        children=[
                            html.Div(
                                children=[
                                    html.Div(
                                        id="graph-title", children=["ocean optics"]
                                    ),
                                    dcc.Graph(id="spec-readings", animate=True),
                                    dcc.Interval(
                                        id="spec-reading-interval",
                                        interval=1 * 1000,
                                        n_intervals=0,
                                    ),
                                ]
                            )
                        ],
                    ),
                ],
                className="two-thirds column right__section",
            ),
        ]
    ),
    html.Div(
        id="page",
        children=[
            # banner
            # html.Div(
            #     id="logo",
            #     title="Dash DAQ by Plotly",
            #     style={
            #         "position": "absolute",
            #         "left": "10px",
            #         "top": "10px",
            #         "zIndex": 100,
            #     },
            #     children=[
            #         html.A(
            #             html.Img(
            #                 src="https://s3-us-west-1.amazonaws.com/plotly-tutorials/excel/dash-daq/dash-daq-logo-by-plotly-stripe+copy.png",
            #                 style={"height": "65px"},
            #             ),
            #             href="http://www.dashdaq.io",
            #         )
            #     ],
            # ),
            # plot
            # power button
            html.Div(
                id="power-button-container",
                title="Turn the power on to begin viewing the data and controlling \
        the spectrometer.",
                children=[
                    daq.PowerButton(
                        id="power-button",
                        size=50,
                        color=colors["accent"],
                        on=True if DEMO else False,
                    )
                ],
            ),
            # status box
            html.Div(
                id="status-box",
                children=[
                    # light intensity
                    # html.Div(
                    #     className="status-box-title", children=["light intensity"]
                    # ),
                    # html.Div(
                    #     id="light-intensity-knob-container",
                    #     title="Controls the intensity of the light source, if any.",
                    #     children=[
                    #         daq.Knob(
                    #             id="light-intensity-knob",
                    #             size=110,
                    #             color=colors["accent"],
                    #             value=0,
                    #         )
                    #     ],
                    # ),
                    # autoscale
                    html.Div(className="status-box-title", children=["autoscale plot"]),
                    html.Div(
                        id="autoscale-switch-container",
                        title="Controls whether the plot automatically resizes \
                to fit the spectra.",
                        children=[
                            daq.BooleanSwitch(
                                id="autoscale-switch", on=True, color=colors["accent"]
                            )
                        ],
                    ),
                    # submit button
                    html.Div(
                        id="submit-button-container",
                        title="Sends all of the control values below the graph \
                to the spectrometer.",
                        children=[
                            html.Button(
                                "update",
                                id="submit-button",
                                n_clicks=0,
                                n_clicks_timestamp=0,
                            )
                        ],
                    ),
                    # displays whether the parameters were successfully changed
                ],
            ),
            # all controls
            # html.Div(
            #     id="controls",
            #     title="All of the spectrometer parameters that can be changed.",
            #     children=[ctrl.create_ctrl_div(True) for ctrl in controls],
            # ),
            # about the app
            html.Div(
                id="infobox",
                children=[
                    html.Div("about this app", id="infobox-title"),
                    dcc.Markdown(
                        dedent(
                            """
            This app was created to act as an interface for an Ocean Optics \
            spectrometer. The options above are used to control various \
            properties of the instrument; the integration time, the number of \
            scans to average over, the strobe and strobe period, and the \
            light source.

            Clicking \"Update\" after putting in the desired settings will \
            result in them being sent to the device. A status message \
            will appear below the button indicating which commands, if any, \
            were unsuccessful; below the unsuccessful commands, a list of \
            successful commands can be found.

            (Note that the box containing the status information is \
            scrollable.)


            The dial labelled \"light intensity\" will affect the current \
            selected light source, if any. The switch labelled \"autoscale \
            plot\" will change the axis limits of the plot to fit all of the \
            data. Please note that the animations and speed of the graph will \
            improve if this feature is turned off, and that it will not be \
            possible to zoom in on any portion of the plot if it is turned \
            on.
            """
                        )
                    ),
                ],
            ),
        ],
    ),
]

app.layout = html.Div(id="main", children=page_layout)


############################
# Callbacks
############################

# disable/enable the update button depending on whether options have changed
@app.callback(
    Output("submit-button", "style"),
    [Input(ctrl.component_attr["id"], ctrl.val_string()) for ctrl in controls]
    + [Input("submit-button", "n_clicks_timestamp")],
)
def update_button_disable_enable(*args):
    now = time.time() * 1000
    disabled = {
        "color": colors["accent"],
        "backgroundColor": colors["background"],
        "cursor": "not-allowed",
    }
    enabled = {
        "color": colors["background"],
        "backgroundColor": colors["accent"],
        "cursor": "pointer",
    }

    # if the button was recently clicked (less than a second ago), then
    # it's safe to say that the callback was triggered by the button; so
    # we have to "disable" it
    if int(now) - int(args[-1]) < 500 and int(args[-1]) > 0:
        return disabled
    else:
        return enabled


# spec model
@app.callback(Output("graph-title", "children"), [Input("power-button", "on")])
def update_spec_model(_):
    return "ocean optics %s" % spec.model()


# disable/enable controls
@app.callback(Output("controls", "children"), [Input("power-button", "on")])
def disable_enable_controls(pwr_on):
    return [ctrl.create_ctrl_div(not pwr_on) for ctrl in controls]


# keep light intensity from resetting, update the value,
# or disable in the event of no light sources
@app.callback(
    Output("light-intensity-knob-container", "children"),
    [Input("light-intensity-knob", "value")],
    state=[State("light-source-input", "value")] + [State("power-button", "on")],
)
def preserve_set_light_intensity(intensity, ls, pwr):
    if ls != "" and ls is not None:
        spec.send_light_intensity(ls, intensity)
    disable = not (pwr and ls != "" and ls is not None)
    return [
        daq.Knob(
            id="light-intensity-knob",
            size=110,
            color=colors["accent"],
            scale={"interval": "1", "labelInterval": "1"},
            disabled=disable,
            value=intensity,
        )
    ]


# send user-selected options to spectrometer
@app.callback(
    Output("submit-status", "children"),
    [Input("submit-button", "n_clicks")],
    state=[State(ctrl.component_attr["id"], ctrl.val_string()) for ctrl in controls]
    + [State("power-button", "on")],
)
def update_spec_params(n_clicks, *args):

    # don't return anything if the device is off
    if not args[-1]:
        return [
            'Press the power button to the top-right of the app, then \
            press the "update" button above to apply your options to \
            the spectrometer.'
        ]

    # dictionary of commands; component id and associated value
    commands = {controls[i].component_attr["id"]: args[i] for i in range(len(controls))}

    failed, succeeded = spec.send_control_values(commands)

    summary = []

    if len(failed) > 0:
        summary.append(
            "The following parameters were not \
        successfully updated: "
        )
        summary.append(html.Br())
        summary.append(html.Br())

        for f in failed:
            # get the name as opposed to the id of each control
            # for readability
            [ctrlName] = [c.ctrl_name for c in controls if c.component_attr["id"] == f]
            summary.append(ctrlName.upper() + ": " + failed[f])
            summary.append(html.Br())

        summary.append(html.Br())
        summary.append(html.Hr())
        summary.append(html.Br())
        summary.append(html.Br())

    if len(succeeded) > 0:
        summary.append("The following parameters were successfully updated: ")
        summary.append(html.Br())
        summary.append(html.Br())

        for s in succeeded:
            [ctrlName] = [c.ctrl_name for c in controls if c.component_attr["id"] == s]
            summary.append(ctrlName.upper() + ": " + succeeded[s])
            summary.append(html.Br())

    return html.Div(summary)


# update the plot
@app.callback(
    Output("spec-readings", "figure"),
    state=[State("power-button", "on"), State("autoscale-switch", "on")],
    inputs=[Input("spec-reading-interval", "n_intervals")],
)
def update_plot(on, auto_range, _):

    traces = []
    wavelengths = []
    intensities = []

    x_axis = {
        "title": "Wavelength (nm)",
        "titlefont": {"family": "Helvetica, sans-serif", "color": colors["secondary"]},
        "tickfont": {"color": colors["tertiary"]},
        "dtick": 100,
        "color": colors["secondary"],
        "gridcolor": colors["grid-colour"],
    }
    y_axis = {
        "title": "Intensity (AU)",
        "titlefont": {"family": "Helvetica, sans-serif", "color": colors["secondary"]},
        "tickfont": {"color": colors["tertiary"]},
        "color": colors["secondary"],
        "gridcolor": colors["grid-colour"],
    }

    if on:
        spectrum = spec.get_spectrum()
        wavelengths = spectrum[0]
        intensities = spectrum[1]
    else:
        wavelengths = numpy.linspace(400, 900, 5000)
        intensities = [0 for wl in wavelengths]

    if on:
        if auto_range:
            x_axis["range"] = [min(wavelengths), max(wavelengths)]
            y_axis["range"] = [min(intensities), max(intensities)]
    traces.append(
        go.Scatter(
            x=wavelengths,
            y=intensities,
            name="Spectrometer readings",
            mode="lines",
            line={"width": 1, "color": colors["accent"]},
        )
    )

    layout = go.Layout(
        height=600,
        font={"family": "Helvetica Neue, sans-serif", "size": 12},
        margin={"t": 20},
        titlefont={
            "family": "Helvetica, sans-serif",
            "color": colors["primary"],
            "size": 26,
        },
        xaxis=x_axis,
        yaxis=y_axis,
        paper_bgcolor=colors["background"],
        plot_bgcolor=colors["background"],
    )

    return {"data": traces, "layout": layout}


if __name__ == "__main__":
    app.run_server(debug=True)
