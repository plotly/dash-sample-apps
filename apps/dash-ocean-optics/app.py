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

# # demo or actual
if ("DASH_PATH_ROUTING" in os.environ) or (
    len(sys.argv) == 2 and sys.argv[1] == "demo"
):
    spec = doos.DemoSpectrometer(spec_lock, comm_lock)
    DEMO = True
else:
    spec = doos.PhysicalSpectrometer(spec_lock, comm_lock)

spec.assign_spec()

app = dash.Dash(__name__)
server = app.server

colors = {
    "background": "#bbbbbb",
    "primary": "#efefef",
    "secondary": "#efefef",
    "tertiary": "#dfdfdf",
    "grid-colour": "#cccccc",
    "accent": "#2222ff",
}

controls = []

# integration time, microseconds
int_time = Control(
    "integration-time",
    "Integration Time (μs)",
    "NumericInput",
    {
        "id": "integration-time-input",
        "className": "control__dropdowns",
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
    "Number of Scans",
    "NumericInput",
    {
        "id": "nscans-to-average-input",
        "className": "control__dropdowns",
        "max": 100,
        "min": 1,
        "size": 150,
        "value": 1,
    },
)
controls.append(nscans_avg)

# strobe
strobe_enable = Control(
    "continuous-strobe-toggle",
    "Strobe",
    "BooleanSwitch",
    {
        "id": "continuous-strobe-toggle-input",
        "className": "control__dropdowns",
        "color": colors["accent"],
        "on": False,
    },
)
controls.append(strobe_enable)

# strobe period
strobe_period = Control(
    "continuous-strobe-period",
    "Strobe pd. (μs)",
    "NumericInput",
    {
        "id": "continuous-strobe-period-input",
        "className": "control__dropdowns",
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
    "Light source",
    "Dropdown",
    {
        "id": "light-source-input",
        "className": "light-source-dropdown",
        "options": spec.light_sources(),
        "placeholder": "select light source",
        "value": "Lamp 2" if DEMO else "",
    },
)
controls.append(light_sources)

base_intro = """This app was created to act as an interface for an Ocean Optics \
            spectrometer. Use controlling elements to control various \
            properties of the instrument; the integration time, the number of \
            scans to average over, the strobe and strobe period, and the \
            light source.
"""

extend_intro = """This app was created to act as an interface for an Ocean Optics 
spectrometer. Use controlling elements to control various 
properties of the instrument; the integration time, the number of 
scans to average over, the strobe and strobe period, and the
light source.

Clicking "Update" after putting in the desired settings will 
result in your parameter settings being sent to the device. A status message 
will appear below the button indicating which commands, if any, 
were unsuccessful; below the unsuccessful commands, a list of
successful commands can be found.
           
The dial labelled "Light intensity" will affect the current 
selected light source, if any. The switch labelled autoscale 
plot will change the axis limits of the plot to fit all of the 
data. Please note that the animations and speed of the graph will 
improve if this autoscale is turned off, and that it will not be 
possible to zoom in on any portion of the plot if it is turned 
on.
"""

page_layout = [
    html.Div(
        [
            html.Div(
                [
                    html.Img(
                        src=app.get_asset_url("dash-daq-logo.png"), className="logo"
                    ),
                    html.Div(
                        id="controls",
                        title="All of the spectrometer parameters that can be changed.",
                        children=[ctrl.create_ctrl_div(True) for ctrl in controls],
                    ),
                    html.Div(html.Label("Light Intensity"), className="control"),
                    html.Div(
                        [
                            html.Div(
                                [
                                    daq.Knob(
                                        id="light-intensity-knob",
                                        size=110,
                                        color="#565656",
                                        value=0,
                                    )
                                ],
                                id="light-intensity-knob-container",
                            )
                        ],
                        className="control",
                    ),
                    html.Div(
                        [
                            html.Div([html.Label("Autoscale Plot")]),
                            html.Div(
                                [
                                    daq.BooleanSwitch(
                                        id="autoscale-switch", color="#565656", on=True
                                    )
                                ],
                                title="Controls whether the plot automatically resizes to fit the spectra.",
                            ),
                        ],
                        className="control autoscale",
                    ),
                    html.Div(
                        [
                            html.Button(
                                "update",
                                id="submit-button",
                                n_clicks=0,
                                n_clicks_timestamp=0,
                            )
                        ],
                        title="Click to send all of the control values to the spectrometer.",
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
                                    html.H6(
                                        id="graph-title", children=["Ocean Optics"]
                                    ),
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
                                    dcc.Markdown(
                                        dedent(base_intro), id="graph-title-intro"
                                    ),
                                    html.Button(
                                        "Learn More", id="learn-more-btn", n_clicks=0
                                    ),
                                    dcc.Graph(
                                        id="spec-readings",
                                        animate=True,
                                        figure=dict(
                                            data=[],
                                            layout=dict(
                                                height=600,
                                                paper_bgcolor="rgba(0,0,0,0)",
                                                plot_bgcolor="rgba(0,0,0,0)",
                                            ),
                                        ),
                                    ),
                                    dcc.Interval(
                                        id="spec-reading-interval",
                                        interval=2
                                        * 1000,  # change from 1 sec to 2 seconds
                                        n_intervals=0,
                                    ),
                                ]
                            )
                        ],
                    )
                ],
                className="two-thirds column right__section",
            ),
        ]
    )
]

app.layout = html.Div(id="main", children=page_layout)


############################
# Callbacks
############################


@app.callback(
    [Output("graph-title-intro", "children"), Output("learn-more-btn", "children")],
    [Input("learn-more-btn", "n_clicks")],
)
def display_info_box(btn_click):
    if (btn_click % 2) == 1:
        return dedent(extend_intro), "Close"
    else:
        return dedent(base_intro), "Learn More"


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
    return "Ocean Optics %s" % spec.model()


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
        margin=dict(l=40, r=40, t=40, b=40, pad=10),
        titlefont={
            "family": "Helvetica, sans-serif",
            "color": colors["primary"],
            "size": 26,
        },
        xaxis=x_axis,
        yaxis=y_axis,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    return {"data": traces, "layout": layout}


if __name__ == "__main__":
    app.run_server(debug=True)
