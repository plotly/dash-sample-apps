# -*- utf-8 -*-

import dash_bootstrap_components as dbc
from dash import (Input, Output, State, html, Dash, dcc, dash_table,
                  get_asset_url, ALL, MATCH, ClientsideFunction,
                  callback_context)
from dash.exceptions import PreventUpdate
#==========================================

#==========================================
import time
import numpy as np
import pandas as pd

#==========================================
from utils.sampling import sampling, df_split, data_split
from utils.modeling import modeling
from utils.charting import prediction_plot, confusion_matrix_plot, roc_curve_plot
#==========================================
from utils.handle_func import *
from utils.split_components import *

#==========================================

#==========================================

#==========================================

#==========================================

#==========================================

#==========================================

#==========================================

app = Dash(
    __name__,
    title='SVM',
    update_title='Eating...',
    #external_scripts=external_scripts,
    external_stylesheets=[dbc.themes.FLATLY])

server = app.server

app.config.suppress_callback_exceptions = True
#==========================================
#==========================================

#=============layout=======================

app.layout = html.Div([
    dbc.Navbar(
        dbc.Container([
            html.A(
                dbc.Row([
                    dbc.Col(
                        html.Img(src=get_asset_url('logo.png'),
                                 height="30px")),
                    dbc.Col(
                        dbc.NavbarBrand("Support Vector Machines",
                                        className="ms-2")),
                ],
                        align="center",
                        className="g-0"),
                href="/",
                style={"textDecoration": "none"},
            ),
        ])),
    dbc.Row([
        dbc.Col(
            dbc.Container([
                html.Br(),
                dbc.Row(dbc.Col(fig_0 := dcc.Graph(), )),
                html.Br(),
                dbc.Row([
                    dbc.Col(fig_1 := dcc.Graph(), width=6, align='center'),
                    dbc.Col(fig_2 := dcc.Graph(), width=6, align='center'),
                ]),
                html.Br(),
                dbc.Row(
                    alert := dbc.Alert(is_open=False,
                                       dismissable=True,
                                       duration=2000,
                                       style={
                                           'padding-left': '45px',
                                           'margin-top': '55px'
                                       }), )
            ], ),
            width=8,
        ),
        dbc.Col(
            dbc.Container([
                html.Br(), uploader_btn,
                html.Br(), threshold,
                html.Br(), kernel,
                html.Br(), formula,
                html.Br(), cost,
                html.Br(), degree,
                html.Br(), gamma,
                html.Br(),
                html.Br(), shrinking,
                html.Br()
            ]),
            width=4,
        )
    ]),
    dbc.Row([
        svm_params := html.Div(style={'display': 'none'}),
        dcc.Store(id='datasets_params', storage_type='memory')
    ])
])

#==========================================


#================callbacks=================
@app.callback(
    [
        Output({
            'type': 'uploader_parameter',
            'index': 'x'
        }, 'options'),
        Output({
            'type': 'uploader_parameter',
            'index': 'y'
        }, 'options'),
        Output({
            'type': 'uploader_parameter',
            'index': 'c'
        }, 'options')
    ],
    [Input({
        'type': 'uploader_parameter',
        'index': 'uploader'
    }, 'contents')],
    [State({
        'type': 'uploader_parameter',
        'index': 'uploader'
    }, 'filename')])
def update_output(data, filename):
    if data is None:
        raise PreventUpdate
    else:
        header = parse_contents(data, filename, header=True)
        return 3 * [header]


@app.callback(Output({
    'type': 'tabs-table',
    'id': 't2'
}, 'data'), [Input({
    'type': 'uploader_parameter',
    'index': ALL
}, 'value')], [
    State({
        'type': 'uploader_parameter',
        'index': ALL
    }, 'id'),
    State({
        'type': 'uploader_parameter',
        'index': 'uploader'
    }, 'contents'),
    State({
        'type': 'uploader_parameter',
        'index': 'uploader'
    }, 'filename')
])
def update_output(xyc, idx, uploaded_data, filename):
    data_params = {j['index']: xyc[i] for i, j in enumerate(idx)}
    if not all([v for k, v in data_params.items() if k not in ['uploader']
                ]) or not uploaded_data:
        raise PreventUpdate
    else:

        df0 = parse_contents(uploaded_data,
                             filename,
                             header=False,
                             usecols=[
                                 v for k, v in data_params.items()
                                 if k not in ['uploader', 'test_size']
                             ])

        df0.rename(
            columns={v: k
                     for k, v in data_params.items() if k != 'test_size'},
            inplace=1)

        data_params.update({'df': df0})

        data = df_split(**data_params)

        split_data = data[0]

        df1 = pd.DataFrame(split_data[0], columns=['x', 'y'])
        df1['c'] = split_data[2]
        df1['s'] = 'TRAIN'
        df2 = pd.DataFrame(split_data[1], columns=['x', 'y'])
        df2['c'] = split_data[3]
        df2['s'] = 'TEST'
        df = pd.concat([df1, df2])

        return df.to_dict('records')


@app.callback(Output({
    'type': 'tabs-table',
    'id': 't3'
}, 'data'), [
    Input({
        'type': 'canvas_parameter',
        'index': 'test_size'
    }, 'value'),
    Input({
        'type': 'canvas_parameter',
        'index': 'canvas'
    }, 'json_data')
])
def canvas_output(test_size, canvas_data):
    if not canvas_data:
        raise PreventUpdate

    X, y = handle_json(canvas_data)

    params = {'test_size': test_size, 'X': X, 'y': y}

    data = data_split(**params)

    split_data = data[0]

    df1 = pd.DataFrame(split_data[0], columns=['x', 'y'])
    df1['c'] = split_data[2]
    df1['s'] = 'TRAIN'
    df2 = pd.DataFrame(split_data[1], columns=['x', 'y'])
    df2['c'] = split_data[3]
    df2['s'] = 'TEST'
    df = pd.concat([df1, df2])

    return df.to_dict('records')


@app.callback(Output({
    'type': 'tabs-table',
    'id': 't1'
}, 'data'), [Input({
    'type': 'dataset_parameter',
    'index': ALL
}, 'value')], [State({
    'type': 'dataset_parameter',
    'index': ALL
}, 'id')])
def generate_data(value, idx):

    data_params = {j['index']: value[i] for i, j in enumerate(idx)}

    data = sampling(**data_params)

    split_data = data[0]

    df1 = pd.DataFrame(split_data[0], columns=['x', 'y'])
    df1['c'] = split_data[2]
    df1['s'] = 'TRAIN'
    df2 = pd.DataFrame(split_data[1], columns=['x', 'y'])
    df2['c'] = split_data[3]
    df2['s'] = 'TEST'
    df = pd.concat([df1, df2])
    return df.to_dict('records')


@app.callback(Output(svm_params, 'children'),
              [Input({
                  'type': 'svm_parameter',
                  'index': ALL
              }, 'value')],
              [State({
                  'type': 'svm_parameter',
                  'index': ALL
              }, 'id')])
def params_update(value, idx):
    df = pd.DataFrame({'index': [i['index'] for i in idx], 'value': value})
    return dash_table.DataTable(df.to_dict('records'), [{
        "name": i,
        "id": i
    } for i in df.columns])


@app.callback([
    Output(fig_0, 'figure'),
    Output(fig_1, 'figure'),
    Output(fig_2, 'figure'),
    Output(alert, 'children'),
    Output(alert, 'is_open')
], [
    Input(save_btn, 'n_clicks'),
    Input({
        'type': 'svm_parameter',
        'index': ALL
    }, 'value')
], [
    State({
        'type': 'svm_parameter',
        'index': ALL
    }, 'id'),
    State({
        'type': 'dataset_parameter',
        'index': ALL
    }, 'id'),
    State({
        'type': 'dataset_parameter',
        'index': ALL
    }, 'value'),
    State({
        'type': 'uploader_parameter',
        'index': ALL
    }, 'id'),
    State({
        'type': 'uploader_parameter',
        'index': ALL
    }, 'value'),
    State({
        'type': 'uploader_parameter',
        'index': ALL
    }, 'contents'),
    State({
        'type': 'uploader_parameter',
        'index': ALL
    }, 'filename'),
    State({
        'type': 'canvas_parameter',
        'index': ALL
    }, 'json_data'),
    State({
        'type': 'canvas_parameter',
        'index': ALL
    }, 'value'),
    State(tabs, 'active_tab'),
    State('datasets_params', 'data')
])
def params_update(n_clicks, value, idx, data_1_idx, data_1_value, data_2_idx,
                  data_2_value, uploaded_data, filename, canvas_data,
                  canvas_params, at, store):
    t1 = time.perf_counter()

    if store and callback_context.triggered[0]["prop_id"].split(
            ".")[0] != save_btn.id:
        [
            data_1_idx, data_1_value, data_2_idx, data_2_value, uploaded_data,
            filename, canvas_data, canvas_params, at
        ] = store

    if at == 'tab-0':
        data_1_params = {
            j['index']: data_1_value[i]
            for i, j in enumerate(data_1_idx)
        }

        data = sampling(**data_1_params)

    elif at == 'tab-1':
        data_2_params = {
            j['index']: data_2_value[i]
            for i, j in enumerate(data_2_idx)
        }

        df0 = parse_contents(list(filter(None, uploaded_data))[0],
                             list(filter(None, filename))[0],
                             header=False,
                             usecols=[
                                 v for k, v in data_2_params.items()
                                 if k not in ['uploader', 'test_size']
                             ])

        df0.rename(columns={
            v: k
            for k, v in data_2_params.items() if k != 'test_size'
        },
                   inplace=1)

        data_2_params.update({'df': df0})

        data = df_split(**data_2_params)

    elif at == 'tab-2':
        X, y = handle_json(list(filter(None, canvas_data))[0], )

        split_params = {
            'test_size': list(filter(None, canvas_params))[0],
            'X': X,
            'y': y
        }

        data = data_split(**split_params)

    if not data:
        raise PreventUpdate

    params = {j['index']: value[i] for i, j in enumerate(idx)}

    params.update({'data': data})
    params.update({'cost': 10**params['cost_power'] * params['cost_coef']})
    params.update({'gamma': 10**params['gamma_power'] * params['gamma_coef']})

    model = modeling(**params)
    params.update({'model': model})

    fig_0 = prediction_plot(**params)
    fig_1 = roc_curve_plot(**params)
    fig_2 = confusion_matrix_plot(**params)

    t2 = time.perf_counter()

    alert_info = 'Takes {:.3} seconds'.format(t2 - t1)

    return fig_0, fig_1, fig_2, alert_info, True


@app.callback(Output(offcanvas, 'is_open'),
              [Input(offcanvas_btn, 'n_clicks'),
               Input(save_btn, 'n_clicks')], [State(offcanvas, 'is_open')])
def open_offcanvas(n, sn, is_open):
    if n or sn:
        return not is_open
    return is_open


@app.callback(Output(save_btn, 'disabled'),
              [Input({
                  'type': 'tabs-table',
                  'id': ALL
              }, 'is_loading')])
def btn_disabled(tabs_table):
    if any(tabs_table):
        return True
    else:
        return False


@app.callback(Output(offcanvas_content, "children"),
              [Input(tabs, "active_tab")])
def switch_tab(at):
    if at == 'tab-0':
        return tab_1_content
    elif at == 'tab-1':
        return tab_2_content
    elif at == 'tab-2':
        return tab_3_content
    return html.P("Something wrong...")


@app.callback(Output({
    'type': 'svm_parameter',
    'index': 'degree'
}, 'disabled'), [Input({
    'type': 'svm_parameter',
    'index': 'kernel'
}, 'value')])
def disable_param_degree(kernel):
    return kernel != 'poly'


@app.callback(Output(latex_formula, 'children'),
              [Input({
                  'type': 'svm_parameter',
                  'index': 'kernel'
              }, 'value')])
def disable_param_degree(kernel):
    ltx = {
        'rbf': r'$K(x, z) = exp(-\gamma||x-z||^2)$',
        'linear': r'$K(x, z) = x \bullet z$',
        'poly': r'$K(x,z) = (\gamma x \bullet z+r)^d$',
        'sigmoid': r'$K(x,z) = tanh(\gamma x \bullet z+r)$'
    }
    return dcc.Markdown(ltx[kernel], mathjax=True)


@app.callback([
    Output({
        'type': 'svm_parameter',
        'index': 'gamma_power'
    }, 'disabled'),
    Output({
        'type': 'svm_parameter',
        'index': 'gamma_coef'
    }, 'disabled')
], [Input({
    'type': 'svm_parameter',
    'index': 'kernel'
}, 'value')])
def disable_param_gamma(kernel):
    _ = kernel not in ['rbf', 'poly', 'sigmoid']
    return _, _


@app.callback(
    Output({
        'type': 'svm_parameter',
        'index': 'cost_coef'
    }, 'marks'),
    [Input({
        'type': 'svm_parameter',
        'index': 'cost_power'
    }, 'value')])
def update_slider_svm_parameter_C_coef(power):
    scale = 10**power
    return {i: str(round(i * scale, 8)) for i in range(1, 10, 2)}


@app.callback(Output({
    'type': 'svm_parameter',
    'index': 'gamma_coef'
}, 'marks'), Input({
    'type': 'svm_parameter',
    'index': 'gamma_power'
}, 'value'))
def scale_param_gamma(power):
    scale = 10**power
    return {i: str(round(i * scale, 8)) for i in range(1, 10, 2)}


@app.callback(Output({
    'type': 'svm_parameter',
    'index': 'threshold'
}, 'value'), [Input(threshold_btn, 'n_clicks')], [State(fig_0, 'figure')])
def reset_threshold(n_clicks, fig):
    if n_clicks:
        Z = np.array(fig['data'][0]['z'])
        value = -Z.min() / (Z.max() - Z.min())
    else:
        value = 0.5
    return value


#==========================================

app.clientside_callback(
    ClientsideFunction(namespace='clientside', function_name='canvas_toggle'),
    Output({
        'type': 'canvas_parameter',
        'index': 'canvas'
    }, 'lineColor'),
    [Input({
        'type': 'canvas_parameter',
        'index': 'toggle'
    }, 'n_clicks')],
    [State({
        'type': 'canvas_parameter',
        'index': 'canvas'
    }, 'lineColor')])

app.clientside_callback(
    ClientsideFunction(namespace='clientside',
                       function_name='datasets_params_store'),
    Output('datasets_params', 'data'),
    [Input(save_btn, 'n_clicks')],
    [
        # data_1_idx, data_1_value, data_2_idx,data_2_value, uploaded_data, filename, canvas_data, canvas_params
        State({
            'type': 'dataset_parameter',
            'index': ALL
        }, 'id'),
        State({
            'type': 'dataset_parameter',
            'index': ALL
        }, 'value'),
        State({
            'type': 'uploader_parameter',
            'index': ALL
        }, 'id'),
        State({
            'type': 'uploader_parameter',
            'index': ALL
        }, 'value'),
        State({
            'type': 'uploader_parameter',
            'index': ALL
        }, 'contents'),
        State({
            'type': 'uploader_parameter',
            'index': ALL
        }, 'filename'),
        State({
            'type': 'canvas_parameter',
            'index': ALL
        }, 'json_data'),
        State({
            'type': 'canvas_parameter',
            'index': ALL
        }, 'value'),
        State(tabs, 'active_tab')
    ])
#==========================================
#==========================================
#==========================================
#==========================================

if __name__ == "__main__":
    app.run_server(debug=True)
