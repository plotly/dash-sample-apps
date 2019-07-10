# In[]:
# Import required libraries
import numpy as np

import plotly.graph_objs as go
import dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State

import dash_daq as daq

from dash_daq_drivers import keithley_instruments

# Instance of a Keithley2400
iv_generator = keithley_instruments.KT2400('COM3', mock_mode=True)

# Define the app
app = dash.Dash('')
server = app.server
app.config.suppress_callback_exceptions = False

# Load css file
external_css = ["https://codepen.io/bachibouzouk/pen/ZRjdZN.css"]
for css in external_css:
    app.css.append_css({"external_url": css})


class UsefulVariables:
    """Class to store information useful to callbacks"""
    def __init__(self):
        self.n_clicks = 0
        self.n_clicks_clear_graph = 0
        self.n_refresh = 0
        self.source = 'V'
        self.is_source_being_changed = False
        self.mode = 'single'
        self.sourced_values = []
        self.measured_values = []

    def change_n_clicks(self, nclicks):
        self.n_clicks = nclicks

    def change_n_clicks_clear_graph(self, nclicks):
        self.n_clicks_clear_graph = nclicks

    def reset_n_clicks(self):
        self.n_clicks = 0
        self.n_clicks_clear_graph = 0

    def change_n_refresh(self, nrefresh):
        self.n_refresh = nrefresh

    def reset_interval(self):
        self.n_refresh = 0

    def clear_graph(self):
        self.sourced_values = []
        self.measured_values = []

    def sorted_values(self):
        """ Sort the data so the are ascending according to the source """
        data_array = np.vstack(
            [
                local_vars.sourced_values,
                local_vars.measured_values
            ]
        )
        data_array = data_array[:, data_array[0, :].argsort()]

        return data_array


local_vars = UsefulVariables()

# font and background colors associated with each themes
bkg_color = {'dark': '#2a3f5f', 'light': '#F3F6FA'}
grid_color = {'dark': 'white', 'light': '#C8D4E3'}
text_color = {'dark': 'white', 'light': '#506784'}


def get_source_labels(source='V'):
    """labels for source/measure elements"""
    if source == 'V':
        # we source voltage and measure current
        source_label = 'Voltage'
        measure_label = 'Current'
    elif source == 'I':
        # we source current and measure voltage
        source_label = 'Current'
        measure_label = 'Voltage'

    return source_label, measure_label


def get_source_units(source='V'):
    """units for source/measure elements"""
    if source == 'V':
        # we source voltage and measure current
        source_unit = 'V'
        measure_unit = 'A'
    elif source == 'I':
        # we source current and measure voltage
        source_unit = 'A'
        measure_unit = 'V'

    return source_unit, measure_unit


h_style = {
    'display': 'flex',
    'flex-direction': 'row',
    'alignItems': 'center',
    'justifyContent': 'space-between',
    'margin': '5px'
}


# Create controls using a function
def generate_main_layout(
        theme='light',
        src_type='V',
        mode_val='single',
        fig=None,
):
    """generate the layout of the app"""

    source_label, measure_label = get_source_labels(src_type)
    source_unit, measure_unit = get_source_units(src_type)

    if mode_val == 'single':
        single_style = {
            'display': 'flex',
            'flex-direction': 'column',
            'alignItems': 'center'
        }
        sweep_style = {'display': 'none'}

        label_btn = 'Single measure'
    else:
        single_style = {'display': 'none'}
        sweep_style = {
            'display': 'flex',
            'flex-direction': 'column',
            'alignItems': 'center'
        }

        label_btn = 'Start sweep'

    # As the trigger-measure btn will have its n_clicks reset by the reloading
    # of the layout we need to reset this one as well
    local_vars.reset_n_clicks()

    # Doesn't clear the data of the graph
    if fig is None:
        data = []
    else:
        data = fig['data']

    html_layout = [
        html.Div(
            className='row',
            children=[
                # graph to trace out the result(s) of the measurement(s)
                html.Div(
                    className="eight columns",
                    children=[
                        dcc.Graph(
                            id='IV_graph',
                            figure={
                                'data': data,
                                'layout': dict(
                                    paper_bgcolor=bkg_color[theme],
                                    plot_bgcolor=bkg_color[theme],
                                    font=dict(
                                        color=text_color[theme],
                                        size=15,
                                    ),
                                    xaxis={
                                        'color': grid_color[theme],
                                        'gridcolor': grid_color[theme]
                                    },
                                    yaxis={
                                        'color': grid_color[theme],
                                        'gridcolor': grid_color[theme]
                                    }
                                )
                            }
                        )
                    ]
                ),
                # controls and options for the IV tracer
                html.Div(
                    className="two columns",
                    id='IV-options',
                    children=[
                        html.H4(
                            'Sourcing',
                            title='Choose whether you want to source voltage '
                                  'and measure current or source current and '
                                  'measure voltage'
                        ),
                        dcc.RadioItems(
                            id='source-choice',
                            options=[
                                {'label': 'Voltage', 'value': 'V'},
                                {'label': 'Current', 'value': 'I'}
                            ],
                            value=src_type
                        ),
                        html.Br(),
                        html.H4(
                            'Measure mode',
                            title='Choose if you want to do single measurement'
                                  ' or to start a sweep'
                        ),
                        dcc.RadioItems(
                            id='mode-choice',
                            options=[
                                {'label': 'Single measure', 'value': 'single'},
                                {'label': 'Sweep', 'value': 'sweep'}
                            ],
                            value=mode_val
                        ),
                        html.Br(),
                        html.Div(
                            daq.StopButton(
                                id='clear-graph_btn',
                                buttonText='Clear graph',
                                size=150
                            ),
                            style={
                                'alignItems': 'center',
                                'display': 'flex',
                                'flex-direction': 'row'
                            }
                        ),
                        html.Br(),
                        daq.Indicator(
                            id='clear-graph_ind',
                            value=False,
                            style={'display': 'none'}
                        )
                    ]
                )
            ]
        ),
        html.Div(
            id='measure_controls',
            className='row',
            children=[
                # Sourcing controls
                html.Div(
                    id='source-div',
                    className="three columns",
                    children=[
                        # To perform single measures adjusting the source with
                        # a knob
                        html.Div(
                            id='single_div',
                            children=[
                                daq.Knob(
                                    id='source-knob',
                                    value=0.00,
                                    min=0,
                                    max=10,
                                    label='%s (%s)' % (
                                        source_label,
                                        source_unit
                                    )
                                ),
                                daq.LEDDisplay(
                                    id="source-knob-display",
                                    label='Knob readout',
                                    value=0.00
                                )
                            ],
                            style=single_style
                        ),
                        # To perfom automatic sweeps of the source
                        html.Div(
                            id='sweep_div',
                            children=[
                                html.Div(
                                    id='sweep-title',
                                    children=html.H4(
                                        "%s sweep:" % source_label
                                    )
                                ),
                                html.Div(
                                    [
                                        'Start',
                                        html.Br(),
                                        daq.PrecisionInput(
                                            id='sweep-start',
                                            precision=4,
                                            label=' %s' % source_unit,
                                            labelPosition='right',
                                            value=0,
                                            style={'margin': '5px'}
                                        ),

                                    ],
                                    title='The lowest value of the sweep',
                                    style=h_style
                                ),
                                html.Div(
                                    [
                                        'Stop',
                                        daq.PrecisionInput(
                                            id='sweep-stop',
                                            precision=4,
                                            label=' %s' % source_unit,
                                            labelPosition='right',
                                            value=10,
                                            style={'margin': '5px'}
                                        )
                                    ],
                                    title='The highest value of the sweep',
                                    style=h_style
                                ),
                                html.Div(
                                    [
                                        'Step',
                                        daq.PrecisionInput(
                                            id='sweep-step',
                                            precision=4,
                                            label=' %s' % source_unit,
                                            labelPosition='right',
                                            value=0.2,
                                            style={'margin': '5px'}
                                        )
                                    ],
                                    title='The increment of the sweep',
                                    style=h_style
                                ),
                                html.Div(
                                    [
                                        'Time of a step',
                                        daq.NumericInput(
                                            id='sweep-dt',
                                            value=0.2,
                                            min=0.01,
                                            style={'margin': '5px'}
                                        ),
                                        's'
                                    ],
                                    title='The time spent on each increment',
                                    style=h_style
                                ),
                                html.Div(
                                    [
                                        daq.Indicator(
                                            id='sweep-status',
                                            label='Sweep active',
                                            value=False
                                        )
                                    ],
                                    title='Indicates if the sweep is running',
                                    style=h_style
                                )
                            ],
                            style=sweep_style
                        )
                    ]
                ),
                # measure button and indicator
                html.Div(
                    id='trigger-div',
                    className="two columns",
                    children=[
                        daq.StopButton(
                            id='trigger-measure_btn',
                            buttonText=label_btn,
                            size=150
                        ),
                        daq.Indicator(
                            id='measure-triggered',
                            value=False,
                            label='Measure active'
                        ),
                    ]
                ),
                # Display the sourced and measured values
                html.Div(
                    id='measure-div',
                    className="five columns",
                    children=[
                        daq.LEDDisplay(
                            id="source-display",
                            label='Applied %s (%s)' % (
                                source_label,
                                source_unit
                            ),
                            value="0.0000"
                        ),
                        daq.LEDDisplay(
                            id="measure-display",
                            label='Measured %s (%s)' % (
                                measure_label,
                                measure_unit
                            ),
                            value="0.0000"
                        )
                    ]
                )

            ],
            style={
                'width': '100%',
                'flexDirection': 'column',
                'alignItems': 'center',
                'justifyContent': 'space-between'
            }
        ),
        html.Div(
            children=[
                html.Div(
                    children=dcc.Markdown('''
**What is this mock app about?**

This is an app to show the graphic elements of Dash DAQ used to create an
interface for an IV curve tracer using a Keithley 2400 SourceMeter. This mock
demo does not actually connect to a physical instrument the values displayed
are generated from an IV curve model for demonstration purposes.

**How to use the mock app**

First choose if you want to source (apply) current or voltage, using the radio
item located on the right of the graph area. Then choose if you want to operate
in a single measurement mode or in a sweep mode.

***Single measurement mode***

Adjust the value of the source with the knob at the bottom of the graph area
and click on the `SINGLE MEASURE` button, the measured value will be displayed.
Repetition of this procedure for different source values will reveal the full
IV curve.

***Sweep mode***

Set the sweep parameters `start`, `stop` and `step` as well as the time
spent on each step, then click on the button `START SWEEP`, the result of the
sweep will be displayed on the graph.

The data is never erased unless the button `CLEAR GRAPH is pressed` or if the
source type is changed.

You can purchase the Dash DAQ components at [
dashdaq.io](https://www.dashdaq.io/)
                    '''),
                    style={
                        'max-width': '600px',
                        'margin': '15px auto 300 px auto',
                        'padding': '40px',
                        'alignItems': 'left',
                        'box-shadow': '10px 10px 5px rgba(0, 0, 0, 0.2)',
                        'border': '1px solid #DFE8F3',
                        'color': text_color[theme],
                        'background': bkg_color[theme]
                    }
                )
            ]
        )
    ]

    if theme == 'dark':
        return daq.DarkThemeProvider(children=html_layout)
    elif theme == 'light':
        return html_layout


root_layout = html.Div(
    id='main_page',
    children=[
        dcc.Location(id='url', refresh=False),
        dcc.Interval(id='refresher', interval=1000000),
        html.Div(
            id='header',
            className='banner',
            children=[
                html.H2('Dash DAQ: IV curve tracer (mock app)'),
                daq.ToggleSwitch(
                    id='toggleTheme',
                    label='Dark/Light layout',
                    size=30,
                    style={'display': 'none'}
                ),
                html.Img(
                    src='https://s3-us-west-1.amazonaws.com/plotly'
                        '-tutorials/excel/dash-daq/dash-daq-logo'
                        '-by-plotly-stripe.png',
                    style={
                        'height': '100',
                        'float': 'right',
                    }
                )
            ],
            style={
                'width': '100%',
                'display': 'flex',
                'flexDirection': 'row',
                'alignItems': 'center',
                'justifyContent': 'space-between',
                'background': '#A2B1C6',
                'color': '#506784'
            }
        ),
        html.Div(
            id='page-content',
            children=generate_main_layout(),
            # className='ten columns',
            style={
                'width': '100%',
            }
        )
    ]
)


# In[]:
# Create app layout
app.layout = root_layout


# In[]:
# Create callbacks
# ======= Dark/light themes callbacks =======
@app.callback(
    Output('page-content', 'children'),
    [
        Input('toggleTheme', 'value')
    ],
    [
        State('source-choice', 'value'),
        State('mode-choice', 'value'),
        State('IV_graph', 'figure')
    ]
)
def page_layout(value, src_type, mode_val, fig):
    """update the theme of the daq components"""

    if value:
        return generate_main_layout('dark', src_type, mode_val, fig)
    else:
        return generate_main_layout('light', src_type, mode_val, fig)


@app.callback(
    Output('page-content', 'style'),
    [Input('toggleTheme', 'value')],
    [State('page-content', 'style')]
)
def page_style(value, style_dict):
    """update the theme of the app"""
    if value:
        theme = 'dark'
    else:
        theme = 'light'

    style_dict['color'] = text_color[theme]
    style_dict['background'] = bkg_color[theme]
    return style_dict


# ======= Callbacks for changing labels =======
@app.callback(
    Output('source-knob', 'label'),
    [
        Input('source-choice', 'value'),
        Input('mode-choice', 'value')
    ],
)
def source_knob_label(src_type, _):
    """update label upon modification of Radio Items"""
    source_label, measure_label = get_source_labels(src_type)
    return source_label


@app.callback(
    Output('source-knob-display', 'label'),
    [
        Input('source-choice', 'value'),
        Input('mode-choice', 'value')
    ],
)
def source_knob_display_label(scr_type, _):
    """update label upon modification of Radio Items"""
    source_label, measure_label = get_source_labels(scr_type)
    source_unit, measure_unit = get_source_units(scr_type)
    return 'Value : %s (%s)' % (source_label, source_unit)


@app.callback(
    Output('sweep-start', 'label'),
    [
        Input('source-choice', 'value'),
        Input('mode-choice', 'value')
    ],
)
def sweep_start_label(src_type, _):
    """update label upon modification of Radio Items"""
    source_unit, measure_unit = get_source_units(src_type)
    return '(%s)' % source_unit


@app.callback(
    Output('sweep-stop', 'label'),
    [
        Input('source-choice', 'value'),
        Input('mode-choice', 'value')
    ],
)
def sweep_stop_label(src_type, _):
    """update label upon modification of Radio Items"""
    source_unit, measure_unit = get_source_units(src_type)
    return '(%s)' % source_unit


@app.callback(
    Output('sweep-step', 'label'),
    [
        Input('source-choice', 'value'),
        Input('mode-choice', 'value')
    ],
)
def sweep_step_label(src_type, _):
    """update label upon modification of Radio Items"""
    source_unit, measure_unit = get_source_units(src_type)
    return '(%s)' % source_unit


@app.callback(
    Output('sweep-title', 'children'),
    [
        Input('source-choice', 'value'),
        Input('mode-choice', 'value')
    ],
)
def sweep_title_label(src_type, _):
    """update label upon modification of Radio Items"""
    source_label, measure_label = get_source_labels(src_type)
    return html.H4("%s sweep:" % source_label)


@app.callback(
    Output('source-display', 'label'),
    [
        Input('source-choice', 'value'),
        Input('mode-choice', 'value')
    ],
)
def source_display_label(src_type, _):
    """update label upon modification of Radio Items"""
    source_label, measure_label = get_source_labels(src_type)
    source_unit, measure_unit = get_source_units(src_type)
    return 'Applied %s (%s)' % (source_label, source_unit)


@app.callback(
    Output('measure-display', 'label'),
    [
        Input('source-choice', 'value'),
        Input('mode-choice', 'value')
    ],
)
def measure_display_label(src_type, _):
    """update label upon modification of Radio Items"""
    source_label, measure_label = get_source_labels(src_type)
    source_unit, measure_unit = get_source_units(src_type)
    return 'Measured %s (%s)' % (measure_label, measure_unit)


@app.callback(
    Output('trigger-measure_btn', 'buttonText'),
    [
        Input('mode-choice', 'value')
    ],
)
def trigger_measure_label(mode_val):
    """update the measure button upon choosing single or sweep"""
    if mode_val == 'single':
        return 'Single measure'
    else:
        return 'Start sweep'


# ======= Callbacks to change elements in the layout =======
@app.callback(
    Output('single_div', 'style'),
    [
        Input('mode-choice', 'value')
    ],
)
def single_div_toggle_style(mode_val):
    """toggle the layout for single measure"""
    if mode_val == 'single':
        return {
            'display': 'flex',
            'flex-direction': 'column',
            'alignItems': 'center'
        }
    else:
        return {'display': 'none'}


@app.callback(
    Output('sweep_div', 'style'),
    [
        Input('mode-choice', 'value')
    ],
)
def sweep_div_toggle_style(mode_val):
    """toggle the layout for sweep"""
    if mode_val == 'single':
        return {'display': 'none'}
    else:
        return {
            'display': 'flex',
            'flex-direction': 'column',
            'alignItems': 'center'
        }


# ======= Applied/measured values display =======
@app.callback(
    Output('source-knob', 'value'),
    [
        Input('source-choice', 'value')
    ],
    [
        State('source-knob', 'value'),
    ]
)
def source_change(src_type, src_val):
    """modification upon source-change
    change the source type in local_vars
    reset the knob to zero
    reset the measured values on the graph
    """
    if src_type == local_vars.source:
        local_vars.is_source_being_changed = False
        return src_val
    else:
        local_vars.is_source_being_changed = True
        local_vars.source = src_type
        return 0.00


# ======= Interval callbacks =======
@app.callback(
    Output('refresher', 'interval'),
    [
        Input('sweep-status', 'value')
    ],
    [
        State('mode-choice', 'value'),
        State('sweep-dt', 'value')
    ]
)
def interval_toggle(swp_on, mode_val, dt):
    """change the interval to high frequency for sweep"""
    if dt <= 0:
        # Precaution against the user
        dt = 0.5
    if mode_val == 'single':
        return 1000000
    else:
        if swp_on:
            return dt * 1000
        else:
            return 1000000


@app.callback(
    Output('refresher', 'n_intervals'),
    [
        Input('trigger-measure_btn', 'n_clicks'),
        Input('mode-choice', 'value')
    ],
    [
        State('sweep-status', 'value'),
        State('refresher', 'n_intervals')
    ]
)
def reset_interval(_, mode_val, swp_on, n_interval):
    """reset the n_interval of the dcc.Interval once a sweep is done"""
    if mode_val == 'single':
        local_vars.reset_interval()
        return 0
    else:
        if swp_on:
            return n_interval
        else:
            local_vars.reset_interval()
            return 0


@app.callback(
    Output('sweep-status', 'value'),
    [
        Input('trigger-measure_btn', 'n_clicks'),
        Input('source-display', 'value')
    ],
    [
        State('measure-triggered', 'value'),
        State('sweep-status', 'value'),
        State('sweep-stop', 'value'),
        State('sweep-step', 'value'),
        State('mode-choice', 'value')
    ]
)
def sweep_activation_toggle(
        _,
        sourced_val,
        meas_triggered,
        swp_on,
        swp_stop,
        swp_step,
        mode_val
):
    """decide whether to turn on or off the sweep
    when single mode is selected, it is off by default
    when sweep mode is selected, it enables the sweep if is wasn't on
    otherwise it stops the sweep once the sourced value gets higher or equal
    than the sweep limit minus the sweep step
    """
    if mode_val == 'single':
        return False
    else:
        if swp_on:
            # The condition of continuation is to source lower than the sweep
            # limit minus one sweep step
            answer = float(sourced_val) <= float(swp_stop)-float(swp_step)
            return answer
        else:
            if not meas_triggered:
                # The 'trigger-measure_btn' wasn't pressed yet
                return False
            else:
                # Initiate a sweep
                return True


# ======= Measurements callbacks =======
@app.callback(
    Output('source-knob-display', 'value'),
    [
        Input('source-knob', 'value')
    ]
)
def set_source_knob_display(knob_val):
    """"set the value of the knob on a LED display"""
    return knob_val


@app.callback(
    Output('measure-triggered', 'value'),
    [
        Input('trigger-measure_btn', 'n_clicks'),
        Input('mode-choice', 'value')
    ]
)
def update_trigger_measure(nclick, _):
    """ Controls if a measure can be made or not
    The indicator 'measure-triggered' can be set to True only by a click
    on the 'trigger-measure_btn' button or by the 'refresher' interval
    """

    if nclick is None:
        nclick = 0

    if int(nclick) != local_vars.n_clicks:
        # It was triggered by a click on the trigger-measure_btn button
        local_vars.change_n_clicks(int(nclick))
        return True
    else:
        # It was triggered by a change of the mode
        return False


@app.callback(
    Output('source-display', 'value'),
    [
        Input('refresher', 'n_intervals'),
        Input('measure-triggered', 'value'),
    ],
    [
        State('source-knob', 'value'),
        State('source-display', 'value'),
        State('sweep-start', 'value'),
        State('sweep-stop', 'value'),
        State('sweep-step', 'value'),
        State('mode-choice', 'value'),
        State('sweep-status', 'value')
    ]
)
def set_source_display(
        n_interval,
        meas_triggered,
        knob_val,
        old_source_display_val,
        swp_start,
        swp_stop,
        swp_step,
        mode_val,
        swp_on
):
    """"set the source value to the instrument"""

    # Default answer
    answer = old_source_display_val

    if mode_val == 'single':
        answer = knob_val
    else:
        if meas_triggered:
            if swp_on:
                answer = float(swp_start) \
                         + (int(n_interval) - 1) * float(swp_step)
                if answer > float(swp_stop):
                    answer = old_source_display_val

    return float("%.4f" % answer)


@app.callback(
    Output('measure-display', 'value'),
    [
        Input('source-display', 'value')
    ],
    [
        State('measure-triggered', 'value'),
        State('measure-display', 'value'),
        State('source-choice', 'value'),
        State('mode-choice', 'value'),
        State('sweep-status', 'value')
    ]
)
def update_measure_display(
        src_val,
        meas_triggered,
        meas_old_val,
        src_type,
        mode_val,
        swp_on
):
    """"read the measured value from the instrument
    check if a measure should be made
    initiate a measure of the KT2400
    read the measure value and return it
    by default it simply return the value previously available
    """
    source_value = float(src_val)
    measured_value = meas_old_val

    if mode_val == 'single':
        if meas_triggered:
            # Save the sourced value
            local_vars.sourced_values.append(source_value)
            # Initiate a measurement
            measured_value = iv_generator.source_and_measure(src_type, src_val)
            # Save the measured value
            local_vars.measured_values.append(measured_value)
    else:
        if meas_triggered and swp_on:
            # Save the sourced value
            local_vars.sourced_values.append(source_value)
            # Initiate a measurement
            measured_value = iv_generator.source_and_measure(src_type, src_val)
            # Save the measured value
            local_vars.measured_values.append(measured_value)

    return measured_value


# ======= Graph related callbacks =======
@app.callback(
    Output('clear-graph_ind', 'value'),
    [
        Input('source-knob', 'value'),
        Input('clear-graph_btn', 'n_clicks'),
        Input('measure-triggered', 'value')
    ],
    []
)
def clear_graph_click(src_val, nclick, meas_triggered):
    """clear the data on the graph
    Uses the callback of the knob value triggered by source-choice change
    or the click on the clear-graph_btn
    everytime a measure is initiated, this value is reset to False, this
    is why we need the input of measure_triggered
    """
    if nclick is None:
        nclick = 0

    if local_vars.is_source_being_changed:
        # The callback was triggered by a source change
        local_vars.is_source_being_changed = False
        local_vars.clear_graph()
        return True
    else:
        if int(nclick) != local_vars.n_clicks_clear_graph:
            # It was triggered by a click on the clear-graph_btn button
            local_vars.change_n_clicks_clear_graph(int(nclick))
            # Reset the data
            local_vars.clear_graph()
            return True
        else:
            return False


@app.callback(
    Output('IV_graph', 'figure'),
    [
        Input('measure-display', 'value'),
        Input('clear-graph_ind', 'value')
    ],
    [
        State('toggleTheme', 'value'),
        State('measure-triggered', 'value'),
        State('IV_graph', 'figure'),
        State('source-choice', 'value'),
        State('mode-choice', 'value'),
        State('sweep-status', 'value')
    ]
)
def update_graph(
        measured_val,
        clear_graph,  # Had to do this because of the lack of multiple Outputs
        theme,
        meas_triggered,
        graph_data,
        src_type,
        mode_val,
        swp_on
):
    """"update the IV graph"""
    if theme:
        theme = 'dark'
    else:
        theme = 'light'

    # Labels for sourced and measured quantities
    source_label, measure_label = get_source_labels(src_type)
    source_unit, measure_unit = get_source_units(src_type)

    if mode_val == 'single':
        if meas_triggered:
            # The change to the graph was triggered by a measure

            # Sort the stored data so the are ascending in x
            data_array = local_vars.sorted_values()

            xdata = data_array[0, :]
            ydata = data_array[1, :]

            data_for_graph = [
                go.Scatter(
                    x=xdata,
                    y=ydata,
                    mode='lines+markers',
                    name='IV curve',
                    line={
                        'color': '#EF553B',
                        'width': 2
                    }
                )
            ]

            return {
                'data': data_for_graph,
                'layout': dict(
                    xaxis={
                        'title': 'Applied %s (%s)' % (
                            source_label, source_unit
                        ),
                        'color': text_color[theme],
                        'gridcolor': grid_color[theme]
                    },
                    yaxis={
                        'title': 'Measured %s (%s)' % (
                            measure_label,
                            measure_unit
                        ),
                        'gridcolor': grid_color[theme]
                    },
                    font=dict(
                        color=text_color[theme],
                        size=15,
                    ),
                    margin={'l': 100, 'b': 100, 't': 50, 'r': 20, 'pad': 0},
                    plot_bgcolor=bkg_color[theme],
                    paper_bgcolor=bkg_color[theme]
                )
            }
        else:
            if clear_graph:
                graph_data['data'] = local_vars.sorted_values()
            return graph_data
    else:
        if swp_on:
            # The change to the graph was triggered by a measure

            # Sort the stored data so the are ascending in x
            data_array = local_vars.sorted_values()

            xdata = data_array[0, :]
            ydata = data_array[1, :]

            data_for_graph = [
                go.Scatter(
                    x=xdata,
                    y=ydata,
                    mode='lines+markers',
                    name='IV curve',
                    line={
                        'color': '#EF553B',
                        'width': 2
                    }
                )
            ]

            return {
                'data': data_for_graph,
                'layout': dict(
                    xaxis={
                        'title': 'Applied %s (%s)' % (
                            source_label, source_unit
                        ),
                        'color': text_color[theme],
                        'gridcolor': grid_color[theme]
                    },
                    yaxis={
                        'title': 'Measured %s (%s)' % (
                            measure_label,
                            measure_unit
                        ),
                        'gridcolor': grid_color[theme]
                    },
                    font=dict(
                        color=text_color[theme],
                        size=15,
                    ),
                    margin={'l': 100, 'b': 100, 't': 50, 'r': 20, 'pad': 0},
                    plot_bgcolor=bkg_color[theme],
                    paper_bgcolor=bkg_color[theme]
                )
            }
        else:
            if clear_graph:
                graph_data['data'] = local_vars.sorted_values()
            return graph_data


# In[]:
# Main
if __name__ == '__main__':
    app.run_server(debug=True)
