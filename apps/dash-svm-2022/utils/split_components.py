import dash_bootstrap_components as dbc
import dash_daq as daq
import dash_canvas
from dash import html, dcc, dash_table

dataset = html.Div([
    html.Strong('Dataset'),
    dbc.RadioItems(id={
        'type': 'dataset_parameter',
        'index': 'dataset'
    },
                   className="btn-group",
                   inputClassName="btn-check",
                   labelClassName="btn btn-outline-primary",
                   labelCheckedClassName="active",
                   options=[
                       {
                           "label": "Moons",
                           "value": 'moons'
                       },
                       {
                           "label": "Linearly Separable",
                           "value": 'LS'
                       },
                       {
                           "label": "Circles",
                           "value": 'circles'
                       },
                   ],
                   value='LS')
])

sample_size = html.Div([
    html.Strong('Sample Size'),
    html.Br(),
    daq.Slider(
        id={
            'type': 'dataset_parameter',
            'index': 'sample_size'
        },
        min=100,
        max=500,
        value=100,
        step=100,
        marks={i * 100: i * 100
               for i in range(6)},
    )
])

noise_level = html.Div([
    html.Strong('Noise Level'),
    html.Br(),
    daq.Slider(id={
        'type': 'dataset_parameter',
        'index': 'noise'
    },
               max=1,
               value=0.4,
               step=0.1,
               marks={i / 5: i / 5
                      for i in range(1, 5)})
],
                       style={'margin-bottom': '15px'})

test_size = html.Div([
    html.Strong('Test Size'),
    html.Br(),
    daq.Slider(id={
        'type': 'dataset_parameter',
        'index': 'test_size'
    },
               max=0.5,
               value=0.25,
               step=0.05,
               marks={i / 10: i / 10
                      for i in range(1, 5)})
],
                     style={'margin-bottom': '15px'})

#==========================================


def reuse_table(i):
    return dash_table.DataTable(
        id={
            'type': 'tabs-table',
            'id': i
        },
        columns=[{
            "name":
            i,
            "id":
            i,
            "type":
            'numeric' if i != 'c' else 'text',
            'format':
            dash_table.Format.Format(precision=4,
                                     scheme=dash_table.Format.Scheme.fixed)
            if i not in ['s', 'c'] else None
        } for i in ['s', 'x', 'y', 'c']],
        #export_format ='csv',
        fixed_rows={'headers': True},
        style_table={
            'height': '300px',
            'overflow': 'auto'
        },
        page_action='none',
        style_cell={
            'width': '{}%'.format(100 / 4),
            'textOverflow': 'ellipsis',
            'overflow': 'hidden'
        })


#==========================================
#==========================================
#==========================================

tab_1_content = html.Div([
    dbc.Card(
        dbc.CardBody([
            dbc.Row(dataset),
            html.Br(),
            dbc.Row([
                dbc.Col(sample_size, width=4),
                dbc.Col(noise_level, width=4),
                dbc.Col(test_size, width=4)
            ], )
        ]),
        className="mt-3",
    ),
    dbc.Card(
        dbc.CardBody([dcc.Loading(reuse_table('t1'))]),
        className="mt-3",
    )
])

tab_2_content = html.Div([
    dbc.Card(
        dbc.CardBody([
            dcc.Upload(
                id={
                    'type': 'uploader_parameter',
                    'index': 'uploader'
                },
                children=html.Div(
                    ['Drag and Drop or ',
                     html.A('Select Files')]),
                style={
                    #'width': '72%',
                    'height': '45px',
                    'lineHeight': '45px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                })
        ]),
        className="mt-3",
    ),
    dbc.Card(dcc.Loading([
        dbc.CardBody(
            dbc.Row([
                dbc.Col(
                    html.Div([
                        html.Strong('x'),
                        html.Br(),
                        dcc.Dropdown(id={
                            'type': 'uploader_parameter',
                            'index': 'x'
                        })
                    ],
                             style={'margin-bottom': '15px'})),
                dbc.Col(
                    html.Div([
                        html.Strong('y'),
                        html.Br(),
                        dcc.Dropdown(id={
                            'type': 'uploader_parameter',
                            'index': 'y'
                        })
                    ],
                             style={'margin-bottom': '15px'})),
                dbc.Col(
                    html.Div([
                        html.Strong('c'),
                        html.Br(),
                        dcc.Dropdown(id={
                            'type': 'uploader_parameter',
                            'index': 'c'
                        })
                    ],
                             style={'margin-bottom': '15px'})),
                dbc.Col(
                    html.Div([
                        html.Strong('Test Size'),
                        html.Br(),
                        daq.Slider(id={
                            'type': 'uploader_parameter',
                            'index': 'test_size'
                        },
                                   max=0.5,
                                   value=0.25,
                                   step=0.05,
                                   marks={i / 10: i / 10
                                          for i in range(1, 5)})
                    ],
                             style={'margin-bottom': '15px'}))
            ]))
    ]),
             className="mt-3"),
    dbc.Card(dbc.CardBody([dcc.Loading(reuse_table('t2'))]), className="mt-3")
])

tab_3_content = html.Div([
    dbc.Card(
        dbc.CardBody([
            dbc.Row([
                dbc.Col(
                    [
                        html.Div(
                            [
                                dash_canvas.DashCanvas(  # This will be switched to the annotation component in the future.
                                    id={
                                        'type': 'canvas_parameter',
                                        'index': 'canvas'
                                    },
                                    filename=
                                    '/assets/canvas_bg.png',  #get_asset_url('bg.png'),
                                    lineWidth=5,
                                    goButtonTitle='Generate',
                                    lineColor='#509188',
                                    #width=canvas_width,
                                    hide_buttons=[
                                        "zoom", "pan", "line", "pencil",
                                        "rectangle", "select"
                                    ])
                            ],
                            className='canvas_container')
                    ],
                    width=8),
                dbc.Col([
                    dbc.Row([
                        html.Div(dbc.Button('Toggle',
                                            id={
                                                'type': 'canvas_parameter',
                                                'index': 'toggle'
                                            }),
                                 style={'margin-top': '45px'})
                    ]),
                    html.Br(),
                    dbc.Row([
                        html.Div([
                            html.Strong('Test Size'),
                            html.Br(),
                            daq.Slider(id={
                                'type': 'canvas_parameter',
                                'index': 'test_size'
                            },
                                       max=0.5,
                                       value=0.25,
                                       step=0.05,
                                       marks={
                                           i / 10: i / 10
                                           for i in range(1, 5)
                                       })
                        ])
                    ])
                ],
                        width=4)
            ])
        ]),
        className="mt-3",
    ),
    dbc.Card(dbc.CardBody([dcc.Loading(reuse_table('t3'))]), className="mt-3")
])

#==========================================
uploader_btn = html.Div([
    offcanvas_btn := dbc.Button("SELECT DATA",
                                outline=True,
                                color="primary",
                                size='lg'), offcanvas :=
    dbc.Offcanvas([
        tabs := dbc.Tabs(
            [
                tab_1 := dbc.Tab(label="Scikit-learn Datasets"), tab_2 :=
                dbc.Tab(label="Upload Data"), tab_3 :=
                dbc.Tab(label="Hand Painted")
            ],
            active_tab="tab-0",
        ), offcanvas_content := html.Div(),
        dbc.Card(
            html.Div([
                dbc.CardBody([save_btn := dbc.Button("SAVE", color="success")])
            ],
                     className="d-grid gap-2 mx-auto"),
            className="mt-3",
        )
    ],
                  placement='end',
                  is_open=False,
                  title='DATA UPLOAD',
                  style={'width': '85%'})
])

#==========================================
threshold = html.Div([
    html.Strong('Threshold'),
    html.Br(),
    daq.Knob(id={
        'type': 'svm_parameter',
        'index': 'threshold'
    },
             min=0,
             max=1,
             value=0.5,
             size=100), threshold_btn := dbc.Button("RESET THRESHOLD")
],
                     style={'margin-bottom': '15px'})

#==========================================
kernel = html.Div([
    html.Strong('Kernel'),
    html.Br(),
    dcc.Dropdown(id={
        'type': 'svm_parameter',
        'index': 'kernel'
    },
                 options={
                     'rbf': 'Radial basis function (RBF)',
                     'linear': 'Linear',
                     'poly': 'Polynomial',
                     'sigmoid': 'Sigmoid'
                 },
                 value='rbf',
                 style={'width': '75%'})
],
                  style={'margin-bottom': '15px'})

formula = html.Div(latex_formula := html.P())

cost = html.Div(
    [
        html.Strong('Cost (C)'),
        html.Br(),
        daq.Slider(id={
            'type': 'svm_parameter',
            'index': 'cost_power'
        },
                   min=-2,
                   max=4,
                   value=0,
                   marks={i: 10**i
                          for i in range(-2, 5)}),
        html.Br(),
        daq.Slider(
            id={
                'type': 'svm_parameter',
                'index': 'cost_coef'
            },
            min=1,
            max=9,
            value=1,
            step=1,
            handleLabel={
                #"showCurrentValue": True,
                "label": "COST"
            })
    ],
    style={'margin-bottom': '15px'})

degree = html.Div([
    html.Strong('Degree'),
    html.Br(),
    daq.Slider(id={
        'type': 'svm_parameter',
        'index': 'degree'
    },
               min=2,
               max=10,
               value=2,
               step=1,
               marks={i: i
                      for i in range(2, 9, 2)})
],
                  style={'margin-bottom': '15px'})

gamma = html.Div(
    [
        html.Strong('Gamma'),
        html.Br(),
        daq.Slider(
            id={
                'type': 'svm_parameter',
                'index': 'gamma_power'
            },
            min=-5,
            max=0,
            value=-1,
            marks={i: 10**i
                   for i in range(-5, 1)},
        ),
        html.Br(),
        daq.Slider(
            id={
                'type': 'svm_parameter',
                'index': 'gamma_coef'
            },
            min=1,
            max=9,
            value=5,
            step=1,
            handleLabel={
                #"showCurrentValue": True,
                "label": "GAMMA",
                "style": {
                    "height": "15px"
                }
            })
    ],
    style={'margin-bottom': '15px'})

shrinking = html.Div([
    html.Strong('Shrinking'),
    dbc.RadioItems(id={
        'type': 'svm_parameter',
        'index': 'shrinking'
    },
                   className="btn-group",
                   inputClassName="btn-check",
                   labelClassName="btn btn-outline-primary",
                   labelCheckedClassName="active",
                   options=[
                       {
                           "label": "Disable",
                           "value": False
                       },
                       {
                           "label": "Enable",
                           "value": True
                       },
                   ],
                   value=True)
])
