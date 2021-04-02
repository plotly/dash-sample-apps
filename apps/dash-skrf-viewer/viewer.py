# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:59:14 2019

@author: ander906
"""
import base64
import io
import os
import json
import re

import dash
import dash_table
import numpy as np
import matplotlib.pyplot as plt
import skrf as rf
import plotly.graph_objs as go
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from flask_caching import Cache
from uuid import uuid4



external_stylesheets = [
    'https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css',
    'https://codepen.io/chriddyp/pen/bWLwgP.css'
    ]

app = dash.Dash(__name__,external_stylesheets=external_stylesheets )
server = app.server
app.config.suppress_callback_exceptions = True



try:
    cache = Cache(app.server, config={
        'CACHE_TYPE': 'redis',
        'CACHE_REDIS_URL': os.environ.get('REDIS_URL', ''),
        'CACHE_THRESHOLD': 20,
        'CACHE_DEFAULT_TIMEOUT': 1200
    })
except RuntimeError:
    cache = Cache(app.server, config={
        'CACHE_TYPE': 'filesystem',

        'CACHE_DIR': os.path.join(os.getcwd(), 'cache'),
        'CACHE_THRESHOLD': 20,
        'CACHE_DEFAULT_TIMEOUT': 1200
    })
write_snp = False

col_header_test = ["m{}".format(i) for i in range(4)]

app.layout = html.Div(className='page',children=[
    html.Div(className='sub_page',children=[
    #html.Div(className='col-2'),
    html.Div( children=[
    html.H3(className='product',children=[
        'Network Viewer',
        html.A(href='http://www.plotly.com',  children=[
        html.Img(src=app.get_asset_url('dash.png'),
                 style={'height':'50px','float':'right'})
        ]),
        html.A(href='http://www.scikit-rf.org',  children=[
            html.Img(src=app.get_asset_url('powered_by_skrf.png'),
                    style={'height':'50px','float':'right'}),
            ]),

    ]),
    
    html.Div([
        html.Div(
            [
                
                html.Div([dcc.Upload(
            id='upload-data',
                children=html.Div([
                    html.A('Upload File(s)',style={'font-size':'16px'}),
                ]),
                style={
                    'width': '80%',
                    'height': '40px',
                    'lineHeight': '40px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                # Allow multiple files to be uploaded
                multiple=True
            ),
            ]),
                dcc.Loading(id="loading-upload",
                            children=[html.Div(id='data-uploading')],
                            type="default"),
                html.Div(
                    dash_table.DataTable(
                        id='uploaded-data-table',
                        style_cell={
                            'whiteSpace': 'normal',
                            'textAlign': 'left'
                        },
                        style_table={
                            'table-layout': 'fixed',
                            'width': '100%',
                        },
                        css=[{
                            'selector': '.dash-cell div.dash-cell-value',
                            'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'
                        }],
                        columns=[{"name": "Loaded Networks",
                                  "id": "data"}],
                        editable=False,
                        row_selectable="multi",
                        selected_rows=[],
                    ), style={'overflowX': 'auto'}
                ),
                
            html.Div([
                    html.Button('Plot', id='button',style={'width':'100%'}),
                    html.Details(id='parameter-control-div', open=False, children=[
                        html.Summary('Configure',style={'font-size':'16px'}),
                        html.Label('Parameters',style={'font-size':'14px','border-bottom':'1px solid #000'}),
                        dcc.RadioItems(
                            id='parm-select',className='inline',
                            options=[
                                {'label': 'S', 'value': 'S'},
                                {'label': 'Y', 'value': 'Y'},
                                {'label': 'Z', 'value': 'Z'}
                                # {'label': 'ABCD Parameters', 'value': 'A'}
                            ],
                            value='S'),
                        html.Label('Axes',style={'font-size':'14px','border-bottom':'1px solid #000'}),
                        dcc.RadioItems(
                            id='axes-select',
                            options=[
                                {'label': 'Magnitude / Phase', 'value': 'MAG'},
                                {'label': 'Real / Imaginary', 'value': 'RI'},
                                #{'label': 'Magnitude / Time Domain', 'value': 'TIME'},
                                {'label': 'Smith', 'value': 'SMITH'}
                            ],
                            value='MAG')]),
                    html.Details(id='port-table-div', children=[
                        html.Summary('Ports',style={'font-size':'16px'}),
                        html.Div(
                            dash_table.DataTable(
                                id='port-table',
                                css=[{
                                    'selector': '.dash-cell div.dash-cell-value',
                                    'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'
                                }],
                                columns=[{"name": "Parameters",
                                          "id": "Parameters"}],
                                editable=False,
                                row_selectable="multi",
                                selected_rows=[],
                            ),
                        ),
                        
                    ]
                                 ),
                    html.Div(id='plot-options'),
                    # html.Hr(),
                ]),    
            ], className='three columns', style={'word-wrap': 'break-word'},
        ),
        
        html.Div([
 
            html.Div([
                html.Div([

                    html.Div(id='output-plot-text'),
                    dcc.Loading(id="loading-plot",
                        children=[
                            html.Div(
                                id='output-plot1',
                                children=[dcc.Graph(id='graph1', figure={'data': []}, style={'display': 'none'})]
                        ),
                            html.Div(
                                id='output-plot2',
                                children=[dcc.Graph(id='graph2', figure={'data': []}, style={'display': 'none'})]
                        ), ],
                        type="default")
                ], className="twelve columns"),
            ], className="row"),
            html.Div(id='currently-plotted',
                     style={'display': 'none'})
        ], id='tabs-content', className='nine columns'),
        ], className='row'
    ),
    html.Div(id='data-uploading1', style={'display': 'none'}),
    html.Div(id='uuid-hidden-div',
             children=str(uuid4()),
             style={'display': 'none'})
    ]),
    ])])




############################################
# Todo: Remove these functions once skrf v15.0 released \
#  (functions implemented natively in https://github.com/scikit-rf/scikit-rf/pull/279).

class TouchstoneEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, complex):
            return np.real(obj), np.imag(obj)  # split into [real, im]
        if isinstance(obj, rf.Frequency):
            return obj.f.tolist()
        return json.JSONEncoder.default(self, obj)

def to_json(network):
    return json.dumps(network, cls=TouchstoneEncoder)

def from_json(obj):
    ntwk = rf.Network(f_unit='hz')
    ntwk.name = obj['name']
    ntwk.comments = obj['comments']
    ntwk.port_names = obj['port_names']
    ntwk.z0 = np.array(obj['_z0'])[..., 0] + np.array(obj['_z0'])[..., 1] * 1j  # recreate complex numbers
    ntwk.s = np.array(obj['_s'])[..., 0] + np.array(obj['_s'])[..., 1] * 1j
    ntwk.f = np.array(obj['_frequency'])
    ntwk.variables = obj['variables']
    return ntwk


################################################

def load_touchstone(content_string: str, filename: str) -> rf.Network:
    """
    Loads s-parameter data into a skrf network object from an uploaded encoded str.

    Parameters
    ----------
    content_string : str
        Encoded string containing touchstone data
    filename : str
        The filename of the touchstone file.
    Returns
    -------
    d : skrf.Network
        A skrf Network object holding s-parameter data from the touchstone file.

    """
    class dataIO(io.BytesIO):
        """Class used to trick skrf into thinking it is reading a file object."""
        _filename: str

        def __init__(self, data, filename):
            super(dataIO, self).__init__(data)
            self.name = filename

        @property
        def name(self):
            """Filename of SnP data."""
            return self._name

        @name.setter
        def name(self, value):
            assert isinstance(value, str)
            self._name = value

    data = dataIO(content_string, filename)
    d = rf.Network()
    d.read_touchstone(data)
    return d


@app.callback([Output('port-table-div', 'open'),
               Output('port-table', 'data')],
              # [Input('tabs-example', 'value')],
              [Input('button', 'n_clicks')],  # "derived_virtual_data"
              [State('uploaded-data-table', "derived_virtual_data"),
               State('uploaded-data-table', "derived_virtual_selected_rows")])
def update_port_table(_, data, selected_rows):
    if not selected_rows:
        return False, []
    else:
        maxports = 0
        for r in selected_rows:
            d = data[r]
            try:
                nport = int(re.search(r'(\d*)-Port Network', d['data']).group(1))
            except AttributeError:
                print("Number of ports not found in print(skrf.Network()).")
                return False, []
            if maxports < nport:
                maxports = nport
        ports = []
        for i in range(maxports):
            for j in range(maxports):
                ports.append({"Parameters": '{}{}'.format(i + 1, j + 1)})
        return True, ports


@app.callback([Output('uploaded-data-table', 'data'),
               Output('data-uploading', 'children')],
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified'),
               State('uuid-hidden-div', 'children')])
def update_s_output(list_of_contents, list_of_names, list_of_dates, uuid):
    if list_of_contents is not None:
        ch = []
        ports = []
        d = {}
        d2 = []
        content_type = []
        content_string = []
        for i, c in enumerate(list_of_contents):
            ct, cs = c.split(',')
            content_type.append(ct)
            content_string.append(cs)
        for c, n in zip(content_string, list_of_names):
            decoded = base64.b64decode(c)
            try:
                data = load_touchstone(decoded, n)
            except Exception as e:
                print(e)
                return html.Div(['There was an error processing this file.']), None
            d2.append({'data': data.__str__()})
            if write_snp:
                d[n] = data.write_touchstone(return_string=True)
            else:
                d[n] = data.__dict__

        cache.set(uuid, json.dumps(d, cls=TouchstoneEncoder))
        return d2, ''
    else:
        raise PreventUpdate


@app.callback(Output('plot-options', 'children'),
              [Input('currently-plotted', 'children')],
              [State('graph1', "figure"),
               State('graph2', "figure")])
def generate_options(plotted_axes_format, fig1, fig2):
    if plotted_axes_format == 'MAG':
        return html.Div()
    elif plotted_axes_format == 'RI':
        return html.Div()
    elif plotted_axes_format == 'TIME':
        fig1_x = fig1['data'][0]['x']
        fig2_x = fig2['data'][0]['x']
        return html.Div(
            html.Details([
                html.Summary('Timegating Options'),
                html.Button(id='timegate-button', children='Time Gate Visible Traces'),
                html.H6("Frequency Gate"),
                dcc.RangeSlider(
                    id='freq-gate-slider',
                    min=min(fig1_x),
                    max=max(fig1_x),
                    step=(max(fig1_x) - min(fig1_x)) / len(fig1_x),
                    value=[min(fig1_x), max(fig1_x)]
                ),
                html.Div(id='freq-gate-slider-value'),
                html.Div('\n'),
                html.H6("Time Gate"),
                dcc.RangeSlider(
                    id='time-gate-slider',
                    min=min(fig2_x),
                    max=max(fig2_x),
                    step=(max(fig2_x) - min(fig2_x)) / len(fig2_x),
                    value=[min(fig2_x), max(fig2_x)]
                ),
                html.Div(id='time-gate-slider-value'),
                html.Div('\n'),
            ])
        )
    elif plotted_axes_format == 'SMITH':
        return html.Div()
    else:
        return html.Div()


for slider in ['freq-gate-slider', 'time-gate-slider']:
    @app.callback(
        Output(f'{slider}-value', 'children'),
        [Input(slider, "value")]
    )
    def update_slider_values(value):
        return f'{value[0]:0.2g} to {value[1]:0.2g}.'


@app.callback(
    [Output('graph1', "figure"),
     Output('graph2', "figure")],
    [Input('timegate-button', 'n_clicks')],
    [State('graph1', "figure"),
     State('graph2', "figure"),
     State('freq-gate-slider', "value"),
     State('time-gate-slider', 'value')
     ])
def time_gate_plot(_, fig1, fig2, fgate, tgate):
    ftracelist = []
    ttracelist = []
    for d in fig1['data']:
        ntwk = rf.Network(f=np.asfarray(d['x']), f_unit='hz', name=d['name'],
                          s=np.asfarray(d['y']), z0=50)
        # print(ntwk)
        # print(ntwk.time_gate(start=tgate[0], stop=tgate[1]).s_mag)
        ftracelist.append(
            go.Scatter(x=ntwk[f"{fgate[0]}-{fgate[1]}"].f,
                       y=ntwk[f"{fgate[0]}-{fgate[1]}"].time_gate(
                           start=tgate[0] * 1e9,
                           stop=tgate[1] * 1e9).s_mag.flatten(),
                       mode='lines',
                       name=f"{d['name']} Gated ({fgate[0]:0.2g}-{fgate[1]:0.2g} Hz," +
                            f" {tgate[0]:0.2g}-{tgate[1]:0.2g} sec)"
                       )
        )
        # [f"{fgate[0]}-{fgate[1]}"]
        ttracelist.append(
            go.Scatter(x=ntwk[f"{fgate[0]}-{fgate[1]}"].frequency.t,
                       y=ntwk[f"{fgate[0]}-{fgate[1]}"].time_gate(
                           start=tgate[0] * 1e9,
                           stop=tgate[1] * 1e9).s_time_mag.flatten(),
                       mode='lines',
                       name=f"{d['name']} Gated ({fgate[0]:0.2g}-{fgate[1]:0.2g} Hz)"
                       )
        )
    [fig1['data'].append(t) for t in ftracelist]
    [fig2['data'].append(t) for t in ttracelist]
    return fig1, fig2


@app.callback(
    [Output('output-plot-text', "children"),
     Output('output-plot1', "children"),
     Output('output-plot2', "children"),
     Output('currently-plotted', 'children')],
    [Input('button', 'n_clicks')],
    #    Input('nsmooth_input', 'n_submit'), Input('nsmooth_input', 'n_blur')],
    [State('parm-select', 'value'),
     State('axes-select', 'value'),
     State('currently-plotted', 'children'),
     State('uploaded-data-table', "derived_virtual_selected_rows"),
     State('uploaded-data-table', "derived_virtual_data"),
     State('port-table', "derived_virtual_selected_rows"),
     State('port-table', "derived_virtual_data"),
     State('uuid-hidden-div', 'children'),
     ])
@cache.memoize()
def update_graph(_, parm, axes_format, plotted_axes_format, selected_ntwk_rows,
                 selected_ntwk_data,
                 selected_rows, selected_data, uuid):
    plotted_axes_format_output = 'None'
    if not selected_ntwk_rows:
        return (html.Div(children=[
                    html.Div(style={'height':'20px'}),
                    '<- Upload a Touchstone File to Start. Then Select Data to Plot',
                    html.Div(style={'height':'150px'}),
                    html.Div(children = ['Originally written by',html.Br(),'Jackson Anderson,',
                                        html.A('ander906@purdue.edu',
                                            href='mailto:ander906@purdue.edu')],)
                ], style={'text-align':'center'}),
                dcc.Graph(id='graph1', figure={'data': []}, style={'display': 'none'}),
                dcc.Graph(id='graph2', figure={'data': []}, style={'display': 'none'}),
                plotted_axes_format_output)
    elif not selected_rows:
        return (html.Div(children='Please Select Ports to Plot.'),
                dcc.Graph(id='graph1', figure={'data': []}, style={'display': 'none'}),
                dcc.Graph(id='graph2', figure={'data': []}, style={'display': 'none'}),
                plotted_axes_format_output)
    else:

        if axes_format == plotted_axes_format:
            plotted_axes_format_output = dash.no_update
        else:
            plotted_axes_format_output = axes_format
        try:
            s_data = json.loads(cache.get(uuid))
        except TypeError:
            return (html.Div(children='Data could not be retrieved from cache. '
                                      'Your session may have timed out.'
                                      'Please re-upload your data and try again.'),
                    dcc.Graph(id='graph1', figure={'data': []}, style={'display': 'none'}),
                    dcc.Graph(id='graph2', figure={'data': []}, style={'display': 'none'}),
                    plotted_axes_format_output)
        data = {}
        for r in selected_ntwk_rows:
            m = re.search(r"Network: '(.*)', ", selected_ntwk_data[r]['data'])
            if m:
                for key, val in s_data.items():
                    if re.search(m.group(1), key):
                        data[key] = val
        # json_data = s_data

        traces1 = []
        traces2 = []
        layout1 = []
        layout2 = []

        for key, val in data.items():
            if write_snp:
                ntwk = load_touchstone(val.encode(), key)
            else:
                ntwk = from_json(val)
            if parm == 'S':
                pass
            elif parm == 'Y':
                ntwk.s = ntwk.y
            elif parm == 'Z':
                ntwk.s = ntwk.z
            elif parm == 'A':
                ntwk.s = ntwk.a
            for i in range(len(ntwk.s[0, :, 0])):
                for j in range(len(ntwk.s[0, 0, :])):
                    for k in selected_rows:
                        if selected_data[k]['Parameters'] == f'{i+1}{j+1}':
                            yvals1 = []
                            yvals2 = []
                            if axes_format == "MAG":
                                yvals1.append(ntwk.s_mag[:, i, j])
                                yvals2.append(ntwk.s_deg[:, i, j])
                            elif axes_format == "RI":
                                yvals1.append(ntwk.s_re[:, i, j])
                                yvals2.append(ntwk.s_im[:, i, j])
                            elif axes_format == "TIME":
                                yvals1.append(ntwk.s_mag[:, i, j])
                                yvals2.append(ntwk.s_time_mag[:, i, j])
                            elif axes_format == "SMITH":
                                # mpl_to_plotly and plotly don't support Smith plots, so data is plotted
                                # in mpl and read in as image. traces1 stores y data, traces2 stores
                                # trace labels for figure
                                traces1.append(ntwk.s[:, i, j])
                                traces2.append(f'{parm}{i + 1}{j + 1}')
                                continue  # skip normal plotly output

                            traces1.append(
                                go.Scattergl(x=ntwk.f, y=yvals1[0], mode='lines',
                                             name='{}{}{} {}'.format(parm, i + 1, j + 1, key)
                                             )
                            )
                            if axes_format == "TIME":
                                traces2.append(
                                    go.Scattergl(x=ntwk.frequency.t, y=yvals2[0],
                                                 mode='lines',
                                                 name='{}{}{} {}'.format(parm, i + 1,
                                                                       j + 1, key)
                                                 )
                                )
                            else:
                                traces2.append(
                                    go.Scattergl(x=ntwk.f, y=yvals2[0],
                                                 mode='lines',
                                                 name='{}{}{} {}'.format(parm, i + 1,
                                                                       j + 1, key)
                                                 )
                                )
                        else:
                            continue

# Define Layouts
        if axes_format == "MAG":
            layout1.append(go.Layout(
                xaxis={'title': 'Frequency [Hz]',
                       'exponentformat': 'SI'},
                yaxis={'type': 'log',
                       'title': '{} Parameter Magnitude'.format(parm)},
                margin={'l': 60, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 1, 'y': 1},
                hovermode='closest'
            ))
            layout2.append(go.Layout(
                xaxis={'title': 'Frequency [Hz]',
                       'exponentformat': 'SI'},
                yaxis={'type': 'linear',
                       'title': '{} Parameter Phase'.format(parm)},
                margin={'l': 60, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 1, 'y': 1},
                hovermode='closest'
            ))
        elif axes_format == "RI":
            layout1.append(go.Layout(
                xaxis={'title': 'Frequency [Hz]',
                       'exponentformat': 'SI'},
                yaxis={'type': 'linear',
                       'title': '{} Parameter Real'.format(parm)},
                margin={'l': 60, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 1, 'y': 1},
                hovermode='closest'
            ))
            layout2.append(go.Layout(
                xaxis={'title': 'Frequency [Hz]',
                       'exponentformat': 'SI'},
                yaxis={'type': 'linear',
                       'title': '{} Parameter Imaginary'.format(parm)},
                margin={'l': 60, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 1, 'y': 1},
                hovermode='closest'
            ))
        elif axes_format == "TIME":
            layout1.append(go.Layout(
                xaxis={'title': 'Frequency [Hz]',
                       'exponentformat': 'SI'},
                yaxis={'type': 'log',
                       'title': '{} Parameter Magnitude'.format(parm)},
                margin={'l': 60, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 1, 'y': 1},
                hovermode='closest'
            ))
            layout2.append(go.Layout(
                xaxis={'title': 'Time [s]',
                       'exponentformat': 'SI'},
                yaxis={'type': 'log',
                       'title': '{} Parameter, Time Domain'.format(parm)},
                margin={'l': 60, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 1, 'y': 1},
                hovermode='closest'
            ))
        elif axes_format == "SMITH":
            # See https://github.com/4QuantOSS/DashIntro/blob/master/notebooks/Tutorial.ipynb
            # for example of encoding mpl figure to image for dash
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            legend_entry = []
            for i, s in enumerate(traces1):
                if parm == 'Y':
                    rf.plotting.plot_smith(s, ax=ax1, label=traces2[i], chart_type='y')
                else:
                    rf.plotting.plot_smith(s, ax=ax1, label=traces2[i])
            out_img = io.BytesIO()
            fig.savefig(out_img, format='png')
            fig.clf()
            plt.close('all')
            out_img.seek(0)  # rewind file
            encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
            return (
                html.Div([
                    # skip normal plotly graph return and return image of plt.figure() instead
                    html.Div('Interactive Smith Chart coming soon...'),
                    html.Img(src="data:image/png;base64,{}".format(encoded))
                ]),
                dcc.Graph(id='graph1', figure={'data': []}, style={'display': 'none'}),
                dcc.Graph(id='graph2', figure={'data': []}, style={'display': 'none'}),
                plotted_axes_format_output
            )
        return (
            "",
            dcc.Graph(id='graph1', figure={'data': traces1, 'layout': layout1[0]}),
            dcc.Graph(id='graph2', figure={'data': traces2, 'layout': layout2[0]}),
            plotted_axes_format_output
        )



if __name__ == '__main__':
    app.run_server(debug=True)
