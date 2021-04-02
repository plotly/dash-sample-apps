import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import numpy as np


import skrf as rf
from skrf.time import detect_span


external_stylesheets = [
    "https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css",
    "https://codepen.io/chriddyp/pen/bWLwgP.css",
]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server


""" FYI  The device under test was made with this

air = rf.air
air.npoints =1001
dut = air.shunt_delay_load(.2,180)**\
    air.line(100)**\
    air.shunt_delay_load(.4,100)**\
    air.line(4000)**\
    air.shunt(air.load(.7))**\
    air.line(1000)**\
    air.impedance_mismatch(1,3)
    #**air.short()
dut.add_noise_polar(.01,3e-5)
"""
dut = rf.Network("dut.s2p")


# pull out some frequency info to set slider bounds on center/span
freq = dut.frequency
dt = freq.t_ns[1] - freq.t_ns[0]
t_max = freq.t_ns.max()


######## ## the APP
app.layout = html.Div(
    className="page",
    children=[
        html.Div(
            className="sub_page",
            children=[
                # html.Div(className='col-2'),
                html.Div(
                    children=[
                        html.H3(
                            className="product",
                            children=[
                                "Time Gating",
                                html.A(
                                    href="http://www.plotly.com",
                                    children=[
                                        html.Img(
                                            src=app.get_asset_url("dash.png"),
                                            style={"height": "50px", "float": "right"},
                                        )
                                    ],
                                ),
                                html.A(
                                    href="http://www.scikit-rf.org",
                                    children=[
                                        html.Img(
                                            src=app.get_asset_url(
                                                "powered_by_skrf.png"
                                            ),
                                            style={"height": "50px", "float": "right"},
                                        ),
                                    ],
                                ),
                            ],
                        ),
                        html.Div(
                            className="row",
                            children=[
                                html.Div(
                                    className="col-4",
                                    children=[
                                        html.Span(
                                            children=[
                                                html.Label(
                                                    "S-parameter: ", className="inline"
                                                ),
                                                dcc.RadioItems(
                                                    id="sparam-radio",
                                                    options=[
                                                        {
                                                            "label": "S11",
                                                            "value": "s11",
                                                        },
                                                        {
                                                            "label": "S21",
                                                            "value": "s21",
                                                        },
                                                        {
                                                            "label": "S12",
                                                            "value": "s12",
                                                        },
                                                        {
                                                            "label": "S22",
                                                            "value": "s22",
                                                        },
                                                    ],
                                                    value="s11",
                                                    labelStyle={
                                                        "display": "inline-block",
                                                        "margin": "6px",
                                                    },
                                                ),
                                            ]
                                        ),
                                        html.Label("Window Type: "),
                                        dcc.Dropdown(
                                            id="window-dropdown",
                                            options=[
                                                {"label": "kaiser", "value": "kaiser"},
                                                {"label": "boxcar", "value": "boxcar"},
                                                {
                                                    "label": "hamming",
                                                    "value": "hamming",
                                                },
                                            ],
                                            value="hamming",
                                        ),
                                        html.Label("Mode:  ", className="inline"),
                                        dcc.RadioItems(
                                            className="inline",
                                            id="mode-radio",
                                            options=[
                                                {"label": "pass", "value": "bandpass"},
                                                {"label": "stop", "value": "bandstop"},
                                            ],
                                            value="bandpass",
                                            labelStyle={
                                                "display": "inline-block",
                                                "margin": "6px",
                                            },
                                        ),
                                    ],
                                ),
                                html.Div(
                                    className="col",
                                    children=[
                                        html.Button(
                                            id="autogate-button", children="AutoGate"
                                        ),
                                        html.Label("Center Frequency(ns)"),
                                        dcc.Slider(
                                            id="center-slider",
                                            min=-t_max,
                                            max=t_max,
                                            value=0,
                                            step=dt * 10,
                                            marks={
                                                str(int(x)): str(int(x))
                                                for x in np.linspace(-t_max, t_max, 13)
                                            },
                                        ),
                                        html.Label("Span Frequency(ns)"),
                                        dcc.Slider(
                                            id="span-slider",
                                            min=0,
                                            max=t_max,
                                            value=10,
                                            step=dt * 10,
                                            marks={
                                                str(int(x)): str(int(x))
                                                for x in np.linspace(-t_max, t_max, 13)
                                            },
                                        ),
                                    ],
                                ),
                            ],
                        ),
                        html.Div(
                            [
                                dcc.Graph(
                                    id="time-graph", figure={"layout": {"height": 300}}
                                ),
                                dcc.Graph(
                                    id="freq-graph", figure={"layout": {"height": 400}}
                                ),
                            ]
                        ),
                    ]
                ),
            ],
        )
    ],
)


@app.callback(
    Output("time-graph", "figure"),
    Output("freq-graph", "figure"),
    Input("sparam-radio", "value"),
    Input("center-slider", "value"),
    Input("span-slider", "value"),
    Input("window-dropdown", "value"),
    Input("mode-radio", "value"),
)
def update_figure(sparam, center, span, window, mode):
    """
    Update figures given a change in s-parameter of interest, or
    gate parameters
    """

    if window == "kaiser":
        window = ("kaiser", 6)
    oneport = dut.__getattribute__(sparam)
    gated = oneport.time_gate(center=center, span=span, window=window, mode=mode)

    # time
    time_fig = go.Figure()
    time_fig.add_trace(
        go.Line(
            x=freq.t_ns, y=oneport.s_time_db.flatten(), line_color="#204a87", name="OG"
        )
    )

    time_fig.add_trace(go.Line(x=freq.t_ns, y=gated.s_time_db.flatten(), name="gated"))
    dat = oneport.s_time_db.flatten()
    span_height = (dat.min() * 1.5, dat.max() * 0.5)
    time_fig.update_layout(
        shapes=[
            dict(
                type="rect",
                xref="x",
                yref="y",
                x0=center - span / 2,
                y0=span_height[0],
                x1=center + span / 2,
                y1=span_height[1],
                fillcolor="lightgray",
                opacity=0.8,
                line_width=0,
                layer="below",
            )
        ]
    )

    time_fig.update_layout(
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        margin=dict(l=20, r=20, t=40, b=20),
        title="Time Domain",
        xaxis_title="Time (ns)",
        yaxis_title="Magnitude (dB)",
        yaxis_range=(-80, 0),
    )

    # frequency
    freq_fig = go.Figure()
    freq_fig.add_trace(
        go.Line(
            x=freq.f_scaled, y=oneport.s_db.flatten(), name="OG", line_color="#204a87"
        )
    )

    freq_fig.add_trace(go.Line(x=freq.f_scaled, y=gated.s_db.flatten(), name="gated"))

    freq_fig.update_layout(
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        margin=dict(l=20, r=20, t=40, b=20),
        title="Frequency",
        xaxis_title="Frequency (%s)" % dut.frequency.unit,
        yaxis_title="Magnitude (dB)",
        yaxis_range=(-60, 10),
    )
    return time_fig, freq_fig


@app.callback(
    Output("center-slider", "value"),
    Output("span-slider", "value"),
    Input("autogate-button", "n_clicks"),
    Input("sparam-radio", "value"),
)
def autogate(autogate, sparam):
    """
    Autogate chooses center/span values based on basic peak search
    """
    oneport = dut.__getattribute__(sparam)
    n = oneport.s_time_mag.argmax()
    center = oneport.frequency.t_ns[n]
    span = detect_span(oneport)
    return center, span


if __name__ == "__main__":
    app.run_server(debug=True)
