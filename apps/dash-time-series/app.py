import pandas as pd
import numpy as np
import warnings
import base64
import io
from pmdarima.arima import auto_arima
from statsmodels.tsa.stattools import acf
from statsmodels.tsa.statespace.sarimax import SARIMAX
import scipy.stats as stats
import plotly.graph_objects as go
import dash
import dash_table as dt
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
pd.options.mode.chained_assignment = None
warnings.filterwarnings('ignore')

external_stylesheets = ['https://fonts.googleapis.com/css?family=Open+Sans:300,400,700',
                        'https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
meta_tags=[{'name': 'viewport', 'content': 'width=device-width'}])

app.css.config.serve_locally = True
app.scripts.config.serve_locally = True
app.config.suppress_callback_exceptions = True

server = app.server

app.title = 'Time Series Analysis'

app.layout = html.Div(children=[

    ####################################################################################################################
    # header
    ####################################################################################################################

    html.Div(children=[

        # title
        html.H5(children=['Time Series Analysis'], style={'display': 'inline-block',
        'vertical-align': 'middle', 'margin-left': '2vw'}),

        # logo
        html.Img(src=app.get_asset_url('dash-logo-new.png'), style={'display': 'inline-block',
        'vertical-align': 'middle', 'margin-left': '68vw', 'width': '10%'}),

    ], className='row', style={'background-color': '#1A3399', 'color': 'white', 'width': '100%'}),

    # first column
    html.Div(children=[

        ################################################################################################################
        # data
        ################################################################################################################

        # title
        html.Label('Data', style={'margin': '1vw 0.25vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'color': '#2547A1'}),

        # tooltip icon
        html.Div(id='data_tooltip', children=['?'], style={'margin': '1vw 0vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'width': '0.8vw', 'height': '0.8vw', 'line-height': '0.8vw', 'font-size': '0.6vw',
        'font-weight': 'bold', 'color': 'white', 'background-color': '#2F5CAB', 'border-radius': '50%',
        'text-align': 'center'}),

        # tooltip text
        dbc.Tooltip(children=['The input file should contain the dates in the first column and the values in the '
        'second column, without column headers. The accepted file formats are CSV, XLS, XLSX, and TXT. If "Sample" '
        'is selected, the results will be generated using a sample time series.'], target='data_tooltip',
        placement='right', style={'line-height': '0.7vw', 'font-size': '0.6vw', 'color': 'white', 'text-align':
        'justify', 'background-color': '#2F5CAB', 'border-radius': '5px', 'padding': '0.25vw 0.25vw 0.25vw 0.25vw',
        'height': '3.5vw', 'width': '16vw'}),

        # checkbox
        dcc.Checklist(id='sample_data', options=[{'label': 'Sample', 'value': 'true'}], value=['true'],
        style={'margin': '1vw 0vw 0.25vw 7vw', 'display': 'inline-block', 'vertical-align': 'middle',
        'color': '#3770B4', 'font-size': '90%'}),

        # upload
        dcc.Upload(id='uploaded_file', children=html.Div([

            html.P('Drag and Drop or ', style={'display': 'inline'}),
            html.A('Select File', style={'display': 'inline', 'text-decoration': 'underline'}),

        ]), className='custom-upload', multiple=False),

        ################################################################################################################
        # transformation
        ################################################################################################################

        # title
        html.Label('Transformation', style={'margin': '1vw 0.25vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'color': '#2547A1'}),

        # tooltip icon
        html.Div(id='transformation_tooltip', children=['?'], style={'margin': '1vw 0vw 0.25vw 0vw',
        'display': 'inline-block', 'vertical-align': 'middle', 'width': '0.8vw', 'height': '0.8vw',
        'line-height': '0.8vw', 'font-size': '0.6vw', 'font-weight': 'bold', 'color': 'white',
        'background-color': '#2F5CAB', 'border-radius': '50%', 'text-align': 'center'}),

        # tooltip text
        dbc.Tooltip(children=['The data must be positive for any of the following transformations to be applied.'],
        target='transformation_tooltip', placement='right', style={'line-height': '0.7vw', 'font-size': '0.6vw',
        'color': 'white', 'text-align': 'justify', 'background-color': '#2F5CAB', 'border-radius': '5px',
        'padding': '0.25vw 0.25vw 0.25vw 0.25vw', 'height': '2vw', 'width': '10vw'}),

        # radio items
        dcc.RadioItems(id='transformation',
                       options=[{'label': 'Logarithm', 'value': 'logarithm'},
                                {'label': 'Square Root', 'value': 'square-root'},
                                {'label': 'Box-Cox', 'value': 'box-cox'},
                                {'label': 'None', 'value': 'none'}],
                       value='none',
                       labelStyle={'font-size': '90%'}),

        ################################################################################################################
        # differencing
        ################################################################################################################

        # title
        html.Label('Differencing Order', style={'margin': '1vw 0.25vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'color': '#2547A1'}),

        # tooltip icon
        html.Div(id='differencing_tooltip', children=['?'], style={'margin': '1vw 0vw 0.25vw 0vw',
        'display': 'inline-block', 'vertical-align': 'middle', 'width': '0.8vw', 'height': '0.8vw',
        'line-height': '0.8vw', 'font-size': '0.6vw', 'font-weight': 'bold', 'color': 'white',
        'background-color': '#2F5CAB', 'border-radius': '50%', 'text-align': 'center'}),

        # tooltip text
        dbc.Tooltip(children=['Must be positive. If zero is selected, no differencing is applied.'],
        target='differencing_tooltip', placement='right', style={'line-height': '0.7vw', 'font-size': '0.6vw',
        'color': 'white', 'text-align': 'justify', 'background-color': '#2F5CAB', 'border-radius': '5px',
        'padding': '0.25vw 0.25vw 0.25vw 0.25vw', 'height': '1.5vw', 'width': '10vw'}),

        # numeric input
        dcc.Input(id='differencing',
                  type='number',
                  min=0,
                  value=0,
                  style={'height': '2vw', 'width': '69.5%', 'font-size': '90%', 'display': 'block'}),

        ################################################################################################################
        # order
        ################################################################################################################

        # title
        html.Label('ARIMA Order', style={'margin': '1vw 0.25vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'color': '#2547A1'}),

        # tooltip icon
        html.Div(id='arima_tooltip', children=['?'], style={'margin': '1vw 0vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'width': '0.8vw', 'height': '0.8vw', 'line-height': '0.8vw', 'font-size': '0.6vw',
        'font-weight': 'bold', 'color': 'white', 'background-color': '#2F5CAB', 'border-radius': '50%',
        'text-align': 'center'}),

        # tooltip text
        dbc.Tooltip(children=['If "Auto" is selected, the optimal ARIMA order will be automatically determined, and '
        'the supplied input values of p, d and q will be ignored.'], target='arima_tooltip', placement='right',
        style={'line-height': '0.7vw', 'font-size': '0.6vw', 'color': 'white', 'text-align': 'justify',
        'background-color': '#2F5CAB', 'border-radius': '5px', 'padding': '0.25vw 0.25vw 0.25vw 0.25vw',
        'height': '2vw', 'width': '14vw'}),

        # checkbox
        dcc.Checklist(id='auto_arima_order', options=[{'label': 'Auto', 'value': 'true'}], value=['true'],
        style={'margin': '1vw 0vw 0.25vw 4.5vw', 'display': 'inline-block', 'vertical-align': 'middle',
        'color':'#3770B4', 'font-size': '90%'}),

        # p
        html.Div(children=[

            html.Label('p', style={'font-size': '90%', 'width': '0.75vw', 'margin': '1vw 0.5vw 0vw 0vw',
            'display': 'inline-block'}),

            dcc.Input(id='p',
                      type='number',
                      min=0,
                      value=0,
                      style={'height': '2vw', 'width': '63%', 'font-size': '90%', 'display': 'inline-block'}),

        ], className='row'),

        # d
        html.Div(children=[

            html.Label('d', style={'font-size': '90%', 'width': '0.75vw', 'margin': '1vw 0.5vw 0vw 0vw',
            'display': 'inline-block'}),

            dcc.Input(id='d',
                      type='number',
                      min=0,
                      value=0,
                      style={'height': '2vw', 'width': '63%', 'font-size': '90%', 'display': 'inline-block'}),

        ], className='row'),

        # q
        html.Div(children=[

            html.Label('q', style={'font-size': '90%', 'width': '0.75vw', 'margin': '1vw 0.5vw 0vw 0vw',
            'display': 'inline-block'}),

            dcc.Input(id='q',
                      type='number',
                      min=0,
                      value=0,
                      style={'height': '2vw', 'width': '63%', 'font-size': '90%', 'display': 'inline-block'}),

        ], className='row'),

        ################################################################################################################
        # seasonal order
        ################################################################################################################

        # title
        html.Label('Seasonal Order', style={'margin': '1vw 0.25vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'color': '#2547A1'}),

        # tooltip icon
        html.Div(id='seasonal_tooltip', children=['?'], style={'margin': '1vw 0vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'width': '0.8vw', 'height': '0.8vw', 'line-height': '0.8vw', 'font-size': '0.6vw',
        'font-weight': 'bold', 'color': 'white', 'background-color': '#2F5CAB', 'border-radius': '50%',
        'text-align': 'center'}),

        # tooltip text
        dbc.Tooltip(children=['If "Auto" is selected, the optimal seasonal order will be automatically determined, '
        'and the supplied input values of P, D and Q will be ignored. The seasonal periodicity s should instead '
        'always be supplied by the user. If s=1, the time series is assumed to be non-seasonal, and the seasonality '
        'parameters are not estimated.'],target='seasonal_tooltip', placement='right', style={'line-height': '0.7vw',
        'font-size': '0.6vw', 'color': 'white', 'text-align': 'justify', 'background-color': '#2F5CAB',
        'border-radius': '5px', 'padding': '0.25vw 0.25vw 0.25vw 0.25vw', 'height': '4vw', 'width': '17vw'}),

        # checkbox
        dcc.Checklist(id='auto_seasonal_order', options=[{'label': 'Auto', 'value': 'true'}], value=[],
        style={'margin': '1vw 0vw 0.25vw 3.5vw', 'display': 'inline-block', 'vertical-align': 'middle',
        'color':'#3770B4', 'font-size': '90%'}),

        # P
        html.Div(children=[

            html.Label('P', style={'font-size': '90%', 'width': '0.75vw', 'margin': '1vw 0.5vw 0vw 0vw',
            'display': 'inline-block'}),

            dcc.Input(id='P',
                      type='number',
                      min=0,
                      value=0,
                      style={'height': '2vw', 'width': '63%', 'font-size': '90%', 'display': 'inline-block'}),

        ], className='row'),

        # D
        html.Div(children=[

            html.Label('D', style={'font-size': '90%', 'width': '0.75vw', 'margin': '1vw 0.5vw 0vw 0vw',
            'display': 'inline-block'}),

            dcc.Input(id='D',
                      type='number',
                      min=0,
                      value=0,
                      style={'height': '2vw', 'width': '63%', 'font-size': '90%', 'display': 'inline-block'}),

        ], className='row'),

        # Q
        html.Div(children=[

            html.Label('Q', style={'font-size': '90%', 'width': '0.75vw', 'margin': '1vw 0.5vw 0vw 0vw',
            'display': 'inline-block'}),

            dcc.Input(id='Q',
                      type='number',
                      min=0,
                      value=0,
                      style={'height': '2vw', 'width': '63%', 'font-size': '90%', 'display': 'inline-block'}),

        ], className='row'),

        # s
        html.Div(children=[

            html.Label('s', style={'font-size': '90%', 'width': '0.75vw', 'margin': '1vw 0.5vw 0vw 0vw',
            'display': 'inline-block'}),

            dcc.Input(id='s',
                      type='number',
                      min=1,
                      value=1,
                      style={'height': '2vw', 'width': '63%', 'font-size': '90%', 'display': 'inline-block'}),

        ], className='row'),

        ################################################################################################################
        # trend
        ################################################################################################################

        # title
        html.Label('Trend', style={'margin': '1vw 0.25vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'color': '#2547A1'}),

        # tooltip icon
        html.Div(id='trend_tooltip', children=['?'], style={'margin': '1vw 0vw 0.25vw 0vw', 'display': 'inline-block',
        'vertical-align': 'middle', 'width': '0.8vw', 'height': '0.8vw', 'line-height': '0.8vw', 'font-size': '0.6vw',
        'font-weight': 'bold', 'color': 'white', 'background-color': '#2F5CAB', 'border-radius': '50%',
        'text-align': 'center'}),

        # tooltip text
        dbc.Tooltip(children=['If "Trend=Constant" the model will include an intercept and no linear time trend.',
        html.Br(), 'If "Trend=Linear" the model will include a linear time trend and no intercept.', html.Br(),
        'If "Trend=Both" the model will include both an intercept and a linear time trend.', html.Br(),
        'If "Trend=None" the model will include no intercept and no linear time trend.'], target='trend_tooltip',
        placement='right', style={'line-height': '0.7vw', 'font-size': '0.6vw', 'color': 'white', 'text-align':
        'justify', 'background-color': '#2F5CAB', 'border-radius': '5px', 'padding': '0.25vw 0.25vw 0.25vw 0.25vw'}),

        # dropdown
        dcc.Dropdown(id='trend',
                     options=[{'label': 'Constant', 'value': 'c'},
                              {'label': 'Linear', 'value': 't'},
                              {'label': 'Both', 'value': 'ct'},
                              {'label': 'None', 'value': 'n'}],
                     multi=False,
                     clearable=False,
                     value='c',
                     style={'width': '83.5%'}),

    ], style={'display': 'inline-block', 'vertical-align': 'top', 'width': '20vw', 'margin': '0vw 0vw 1vw 2vw'}),

    # second column
    html.Div(children=[

        ################################################################################################################
        # estimated parameters
        ################################################################################################################

        # title
        html.Label('Estimated Parameters', style={'color': '#2547A1', 'margin': '1vw 0vw 0.25vw 0vw'}),

        # table
        html.Div(id='model_estimation', children=[

            # placeholder
            html.Div(children=['Updating...'], style={'color': '#9d9d9d', 'font-size': '90%',
            'position': 'relative', 'top': '8.25vw', 'left': '15vw'})

        ], style={'height': '18.5vw', 'width': '37vw'}),

        ################################################################################################################
        # fitted values
        ################################################################################################################

        # title
        html.Label('Fitted Values', style={'color': '#2547A1', 'margin': '2vw 0vw 0.5vw 0vw'}),

        # plot
        html.Div(id='model_predictions', children=[

            # placeholder
            html.Div(children=['Updating...'], style={'color': '#9d9d9d', 'font-size': '90%',
            'position': 'relative', 'top': '8vw', 'left': '15vw'})

        ], style={'height': '20vw', 'width': '37vw'}),

    ], style={'display': 'inline-block', 'vertical-align': 'top', 'width': '37vw', 'margin': '0vw 0vw 0vw -2vw'}),

    # third column
    html.Div(children=[

        ################################################################################################################
        # autocorrelation function
        ################################################################################################################

        # title
        html.Label('Autocorrelation Function of Residuals', style={'color': '#2547A1', 'margin': '1vw 0vw 0.5vw 0vw'}),

        # plot
        html.Div(id='residuals_acf', children=[

            # placeholder
            html.Div(children=['Updating...'], style={'color': '#9d9d9d', 'font-size': '90%',
            'position': 'relative', 'top': '8vw', 'left': '15vw'})

        ], style={'height': '20vw', 'width': '37vw'}),

        ################################################################################################################
        # probability plot
        ################################################################################################################

        # title
        html.Label('Probability Plot of Residuals', style={'color': '#2547A1', 'margin': '0.25vw 0vw 0.5vw 0vw'}),

        # plot
        html.Div(id='residuals_probplot', children=[

            # placeholder
            html.Div(children=['Updating...'], style={'color': '#9d9d9d', 'font-size': '90%',
            'position': 'relative', 'top': '8vw', 'left': '15vw'})

        ], style={'height': '20vw', 'width': '37vw'}),

    ], style={'display': 'inline-block', 'vertical-align': 'top', 'width': '37vw', 'margin': '0vw 0vw 0vw 3vw'}),

])

@app.callback([Output('model_estimation', 'children'), Output('model_predictions', 'children'),
    Output('residuals_acf', 'children'), Output('residuals_probplot', 'children')], [Input('sample_data', 'value'),
    Input('uploaded_file', 'contents'), Input('uploaded_file', 'filename'), Input('transformation', 'value'),
    Input('differencing', 'value'), Input('auto_arima_order', 'value'), Input('p', 'value'), Input('d', 'value'),
    Input('q', 'value'), Input('auto_seasonal_order', 'value'), Input('P', 'value'), Input('D', 'value'),
    Input('Q', 'value'), Input('s', 'value'), Input('trend', 'value')])
def update_results(sample_data, uploaded_file, filename, transformation, differencing, auto_arima_order, p, d, q,
    auto_seasonal_order, P, D, Q, s, trend):

    ####################################################################################################################
    # load the data
    ####################################################################################################################

    if sample_data == ['true'] or uploaded_file is not None:

        if sample_data == ['true']:

            # load the sample data
            df = pd.read_csv('data/sample_data.csv', header=None)

        else:

            # load the provided data
            content_type, content_string = uploaded_file.split(',')
            decoded = base64.b64decode(content_string)

            if 'csv' in filename:

                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), header=None)

            elif 'xls' in filename:

                df = pd.read_excel(io.BytesIO(decoded), header=None)

            elif 'txt' in filename:

                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), delimiter=r'\s+', header=None)

        ################################################################################################################
        # extract the time series
        ################################################################################################################

        try:

            y = pd.Series(data=df.iloc[:, 1].astype(float).values, index=pd.DatetimeIndex(pd.to_datetime(df.iloc[:, 0],
            infer_datetime_format=True))).dropna()

        except:

            model_estimation = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
            'pre', 'color': '#8d8d8d'}), 'The input file should contain the dates in the first column and the values ',
            html.Br(), 'in the second column, without column headers. The accepted file formats ', html.Br(), 'are CSV,'
            ' XLS, XLSX, and TXT.'], style={'color': '#9d9d9d', 'font-size': '90%', 'text-align': 'justify',
            'padding': '8.25vw 3vw 0vw 3vw'})

            model_predictions = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
            'pre', 'color': '#8d8d8d'}), 'The input file should contain the dates in the first column and the values ',
            html.Br(), 'in the second column, without column headers. The accepted file formats ', html.Br(), 'are CSV,'
            ' XLS, XLSX, and TXT.'], style={'color': '#9d9d9d', 'font-size': '90%', 'text-align': 'justify',
            'padding': '8.25vw 3vw 0vw 3vw'})

            residuals_acf = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
            'pre', 'color': '#8d8d8d'}), 'The input file should contain the dates in the first column and the values ',
            html.Br(), 'in the second column, without column headers. The accepted file formats ', html.Br(), 'are CSV,'
            ' XLS, XLSX, and TXT.'], style={'color': '#9d9d9d', 'font-size': '90%', 'text-align': 'justify',
            'padding': '8.25vw 3vw 0vw 3vw'})

            residuals_probplot = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
            'pre', 'color': '#8d8d8d'}), 'The input file should contain the dates in the first column and the values ',
            html.Br(), 'in the second column, without column headers. The accepted file formats ', html.Br(), 'are CSV,'
            ' XLS, XLSX, and TXT.'], style={'color': '#9d9d9d', 'font-size': '90%', 'text-align': 'justify',
            'padding': '8.25vw 3vw 0vw 3vw'})

            return [model_estimation, model_predictions, residuals_acf, residuals_probplot]

        ################################################################################################################
        # transform the time series
        ################################################################################################################

        if transformation != 'none':

            if y.min() > 0:

                if transformation == 'logarithm':

                    y.iloc[:] = np.log(y)

                elif transformation == 'square-root':

                    y.iloc[:] = np.sqrt(y)

                elif transformation == 'box-cox':

                    y.iloc[:] = stats.boxcox(y.values)[0]

            else:

                model_estimation = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
                'pre', 'color': '#8d8d8d'}), 'The data must be positive for the selected transformation to be applied.'],
                style={'color': '#9d9d9d', 'font-size': '90%', 'position': 'relative', 'top': '8.25vw', 'left': '5vw'})

                model_predictions = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
                'pre', 'color': '#8d8d8d'}), 'The data must be positive for the selected transformation to be applied.'],
                style={'color': '#9d9d9d', 'font-size': '90%', 'position': 'relative', 'top': '8.25vw', 'left': '5vw'})

                residuals_acf = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
                'pre', 'color': '#8d8d8d'}), 'The data must be positive for the selected transformation to be applied.'],
                style={'color': '#9d9d9d', 'font-size': '90%', 'position': 'relative', 'top': '8.25vw', 'left': '5vw'})

                residuals_probplot = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
                'pre', 'color': '#8d8d8d'}), 'The data must be positive for the selected transformation to be applied.'],
                style={'color': '#9d9d9d', 'font-size': '90%', 'position': 'relative', 'top': '8.25vw', 'left': '5vw'})

                return [model_estimation, model_predictions, residuals_acf, residuals_probplot]

        ################################################################################################################
        # difference the time series
        ################################################################################################################

        # difference the time series
        if differencing != 0:

            y = y.diff(differencing).dropna()

            if len(y) == 0:

                model_estimation = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
                'pre', 'color': '#8d8d8d'}), 'Not enough data.'], style={'color': '#9d9d9d', 'font-size': '90%',
                'position': 'relative', 'top': '8.25vw', 'left': '15vw'})

                model_predictions = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
                'pre', 'color': '#8d8d8d'}), 'Not enough data.'], style={'color': '#9d9d9d', 'font-size': '90%',
                'position': 'relative', 'top': '8.25vw', 'left': '15vw'})

                residuals_acf = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
                'pre', 'color': '#8d8d8d'}), 'Not enough data.'], style={'color': '#9d9d9d', 'font-size': '90%',
                'position': 'relative', 'top': '8.25vw', 'left': '15vw'})

                residuals_probplot = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
                'pre', 'color': '#8d8d8d'}), 'Not enough data.'], style={'color': '#9d9d9d', 'font-size': '90%',
                'position': 'relative', 'top': '8.25vw', 'left': '15vw'})

                return [model_estimation, model_predictions, residuals_acf, residuals_probplot]

        ################################################################################################################
        # fit the SARIMA model
        ################################################################################################################
        try:

            # both the ARIMA order and the seasonal order are automatically identified;
            # the user has indicated that the time series is seasonal (s > 1)
            if auto_arima_order == ['true'] and auto_seasonal_order == ['true'] and s > 1:

                if trend == 'n':

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=False)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='n', order=opt.order, seasonal_order=opt.seasonal_order).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                elif trend == 'c':

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=True)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='c', order=opt.order, seasonal_order=opt.seasonal_order).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                elif trend == 't':

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=False, trend='t')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='t', order=opt.order, seasonal_order=opt.seasonal_order).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                else:

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=True, trend='ct')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='ct', order=opt.order, seasonal_order=opt.seasonal_order).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

            # the ARIMA order is automatically identified; the user has selected that the seasonal order
            # should also be automatically identified; however, the user has indicated that the time
            # series is non-seasonal (s = 1), and therefore the seasonality will be ignored
            elif auto_arima_order == ['true'] and auto_seasonal_order == ['true'] and s == 1:

                if trend == 'n':

                    # find the optimal order
                    opt = auto_arima(y=y, seasonal=False, with_intercept=False)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='n', order=opt.order, seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                elif trend == 'c':

                    # find the optimal order
                    opt = auto_arima(y=y, seasonal=False, with_intercept=True)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='c', order=opt.order, seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                elif trend == 't':

                    # find the optimal order
                    opt = auto_arima(y=y, seasonal=False, with_intercept=False, trend='t')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='t', order=opt.order, seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                else:

                    # find the optimal order
                    opt = auto_arima(y=y, seasonal=False, with_intercept=True, trend='ct')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='ct', order=opt.order, seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

            # the ARIMA order is automatically identified, while the seasonal order is provided
            # by the user; the user has indicated that the time series is seasonal (s > 1)
            # (this option is not recommended)
            elif auto_arima_order == ['true'] and auto_seasonal_order != ['true'] and s > 1:

                if trend == 'n':

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=False)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='n', order=opt.order, seasonal_order=(P, D, Q, s)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                elif trend == 'c':

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=True)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='c', order=opt.order, seasonal_order=(P, D, Q, s)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                elif trend == 't':

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=False, trend='t')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='t', order=opt.order, seasonal_order=(P, D, Q, s)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                else:

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=True, trend='ct')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='ct', order=opt.order, seasonal_order=(P, D, Q, s)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

            # the ARIMA order is automatically identified, while the seasonal order is provided
            # by the user; however, the user has indicated that the time series is non-seasonal
            # (s = 1), and therefore the seasonality will be ignored
            elif auto_arima_order == ['true'] and auto_seasonal_order != ['true'] and s == 1:

                if trend == 'n':

                    # find the optimal order
                    opt = auto_arima(y=y, seasonal=False, with_intercept=False)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='n', order=opt.order, seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                elif trend == 'c':

                    # find the optimal order
                    opt = auto_arima(y=y, seasonal=False, with_intercept=True)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='c', order=opt.order, seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                elif trend == 't':

                    # find the optimal order
                    opt = auto_arima(y=y, seasonal=False, with_intercept=False, trend='t')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='t', order=opt.order, seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

                else:

                    # find the optimal order
                    opt = auto_arima(y=y, seasonal=False, with_intercept=True, trend='ct')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='ct', order=opt.order, seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = opt.order[1]

            # the ARIMA order is provided by the user, while the seasonal order is automatically
            # identified; the user has indicated that the time series is seasonal (s > 1)
            # (this option is not recommended)
            elif auto_arima_order != ['true'] and auto_seasonal_order == ['true'] and s > 1:

                if trend == 'n':

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=False)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='n', order=(p, d, q), seasonal_order=opt.seasonal_order).fit(disp=0)

                    # extract the differencing order
                    diff = d

                elif trend == 'c':

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=True)

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='c', order=(p, d, q), seasonal_order=opt.seasonal_order).fit(disp=0)

                    # extract the differencing order
                    diff = d

                elif trend == 't':

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=False, trend='t')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='t', order=(p, d, q), seasonal_order=opt.seasonal_order).fit(disp=0)

                    # extract the differencing order
                    diff = d

                else:

                    # find the optimal order
                    opt = auto_arima(y=y, m=s, with_intercept=True, trend='ct')

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='ct', order=(p, d, q), seasonal_order=opt.seasonal_order).fit(disp=0)

                    # extract the differencing order
                    diff = d

            # the ARIMA order is provided by the user; the user has selected that the seasonal order
            # should be automatically identified; however, the user has indicated that the time series
            # is non-seasonal (s = 1), and therefore the seasonality will be ignored
            elif auto_arima_order != ['true'] and auto_seasonal_order == ['true'] and s == 1:

                if trend == 'n':

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='n', order=(p, d, q), seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = d

                elif trend == 'c':

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='c', order=(p, d, q), seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = d

                elif trend == 't':

                    # estimate the model
                    res = SARIMAX(endog=y, trend='t', order=(p, d, q), seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = d

                else:

                    # estimate the model
                    res = SARIMAX(endog=y, trend='ct', order=(p, d, q), seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = d

            # both the ARIMA order and the seasonal order are provided by the user;
            # the user has indicated that the time series is seasonal (s > 1)
            elif auto_arima_order != ['true'] and auto_seasonal_order != ['true'] and s > 1:

                if trend == 'n':

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='n', order=(p, d, q), seasonal_order=(P, D, Q, s)).fit(disp=0)

                    # extract the differencing order
                    diff = d

                elif trend == 'c':

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='c', order=(p, d, q), seasonal_order=(P, D, Q, s)).fit(disp=0)

                    # extract the differencing order
                    diff = d

                elif trend == 't':

                    # estimate the model
                    res = SARIMAX(endog=y, trend='t', order=(p, d, q), seasonal_order=(P, D, Q, s)).fit(disp=0)

                    # extract the differencing order
                    diff = d

                else:

                    # estimate the model
                    res = SARIMAX(endog=y, trend='ct', order=(p, d, q), seasonal_order=(P, D, Q, s)).fit(disp=0)

                    # extract the differencing order
                    diff = d

            # both the ARIMA order and the seasonal order are provided by the user;
            # however, the user has indicated that the time series is non-seasonal
            # (s = 1), and therefore the seasonality will be ignored
            else:

                if trend == 'n':

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='n', order=(p, d, q), seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = d

                elif trend == 'c':

                    # estimate the parameters
                    res = SARIMAX(endog=y, trend='c', order=(p, d, q), seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = d

                elif trend == 't':

                    # estimate the model
                    res = SARIMAX(endog=y, trend='t', order=(p, d, q), seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = d

                else:

                    # estimate the model
                    res = SARIMAX(endog=y, trend='ct', order=(p, d, q), seasonal_order=(0, 0, 0, 0)).fit(disp=0)

                    # extract the differencing order
                    diff = d

            # extract the model predictions
            predictions = res.get_prediction()

            # extract the model residuals
            residuals = res.resid

            # organize the model results
            results = pd.DataFrame(data={'Actual': y,
                                         'Predicted': predictions.predicted_mean,
                                         'Residual': residuals,
                                         '95% CI (Lower Bound)': predictions.conf_int(alpha=0.025).iloc[:, 0],
                                         '95% CI (Upper Bound)': predictions.conf_int(alpha=0.025).iloc[:, 1]},
                                   index=y.index)

            results = results.iloc[diff:, :]

            # extract the estimated parameters
            table_1 = pd.read_html(res.summary().tables[0].as_html(), header=None, index_col=None)[0]
            table_2 = pd.read_html(res.summary().tables[1].as_html(), header=None, index_col=None)[0]
            table_3 = pd.read_html(res.summary().tables[2].as_html(), header=None, index_col=None)[0]

            ############################################################################################################
            # generate the table with the estimated parameters
            ############################################################################################################

            # extract the model order
            model_order = table_1.iloc[1,1]

            # organize the first table
            table_1.iloc[1, 0] = 'No. Observations'
            table_1.iloc[1, 1] = len(y)
            table_1.iloc[2, 0] = 'BIC'
            table_1.iloc[2, 1] = format(table_1.iloc[3, 3], '.4f')
            table_1 = table_1.iloc[1:3, :]
            table_1.reset_index(inplace=True, drop=True)
            table_1.iloc[:, 3] = table_1.iloc[:, 3].astype(float).apply(lambda x: format(x, '.4f'))

            # organize the second table
            table_2 = table_2.iloc[:, [0, 1, 2, 4]].rename(columns={4: 3})
            table_2.iloc[1:, 1] = table_2.iloc[1:, 1].astype(float).apply(lambda x: format(x, '.4f'))
            table_2.iloc[1:, 2] = table_2.iloc[1:, 2].astype(float).apply(lambda x: format(x, '.4f'))
            table_2.iloc[1:, 3] = table_2.iloc[1:, 3].astype(float).apply(lambda x: format(x, '.4f'))

            # organize the third table
            table_3.iloc[:, 1] = table_3.iloc[:, 1].astype(float).apply(lambda x: format(x, '.4f'))
            table_3.iloc[:, 3] = table_3.iloc[:, 3].astype(float).apply(lambda x: format(x, '.4f'))

            # merge the three tables
            table = table_1.append(table_2).append(table_3)
            table.reset_index(inplace=True, drop=True)

            # generate the table with the estimated parameters
            model_estimation = html.Div(children=[

                html.Label(children=[model_order], style={'font-size': '90%', 'text-align': 'center',
                'margin': '0vw 0vw 0.25vw 0vw'}),

                dt.DataTable(data=table.to_dict(orient='records'),
                    columns=[{'id': str(x), 'name': str(x)} for x in list(table.columns)],
                    style_header={'display': 'none'}, style_table={'width': '36vw', 'min-width': '36vw',
                    'max-width': '36vw', 'overflow-x': 'hidden', 'overflow-y': 'hidden', 'border-top': 'double #9d9d9d',
                    'border-bottom': 'double #9d9d9d'},
                    style_cell={'font-family': 'Open Sans', 'font-size': '80%', 'text-align': 'left', 'width': '4vw',
                    'min-width': '4vw', 'max-width': '4vw', 'height': '1vw', 'min-height': '1vw', 'max-height': '1vw'},
                    style_data_conditional=[{'if': {'row_index': 2}, 'border-top': 'double #9d9d9d'},
                    {'if': {'row_index': table.shape[0] - 4}, 'border-top': 'double #9d9d9d'}])

            ], style={'height': '18.5vw', 'width': '37vw', 'overflow-x': 'auto', 'overflow-y': 'auto'})

            ############################################################################################################
            # generate the plot of the model predictions
            ############################################################################################################

            # generate the layout
            layout = go.Layout(font=dict(family='Open Sans', color='#737373'),
                               paper_bgcolor='white',
                               plot_bgcolor='white',
                               margin=dict(l=2, r=2, t=40, b=2, pad=0),
                               xaxis=dict(type='date', showgrid=False, zeroline=False, mirror=True, color='#737373',
                               linecolor='#d9d9d9'),
                               yaxis=dict(showgrid=False, zeroline=False, mirror=True, color='#737373',
                               linecolor='#d9d9d9'),
                               legend=dict(x=0, y=1.2, orientation='h', traceorder='reversed'))

            # generate the traces
            data = []

            data.append(go.Scatter(x=list(results.index),
                                   y=list(results['95% CI (Upper Bound)']),
                                   mode='lines',
                                   legendgroup='group',
                                   showlegend=False,
                                   line=dict(width=0.5, shape='spline', color='#4087BD'),
                                   hovertemplate='<b>95% Confidence Interval, Upper Bound</b><br>'
                                   'Date: %{x|%d %B %Y}<br>Value: %{y}<extra></extra>'))

            data.append(go.Scatter(x=list(results.index),
                                   y=list(results['95% CI (Lower Bound)']),
                                   mode='lines',
                                   legendgroup='group',
                                   name='95% Confidence Interval',
                                   line=dict(width=0.5, shape='spline', color='#4087BD'),
                                   fill='tonexty',
                                   fillcolor='rgba(64,135,189,0.1)',
                                   hovertemplate='<b>95% Confidence Interval, Lower Bound</b><br>'
                                   'Date: %{x|%d %B %Y}<br>Value: %{y}<extra></extra>'))

            data.append(go.Scatter(x=list(results.index),
                                   y=list(results['Predicted']),
                                   name='Predicted',
                                   mode='lines',
                                   line=dict(color='#4A9DC7', shape='spline', dash='dot', width=2),
                                   hovertemplate='<b>Predicted</b><br>Date: %{x|%d %B %Y}<br>'
                                   'Value: %{y}<extra></extra>'))

            data.append(go.Scatter(x=list(results.index),
                                   y=list(results['Actual']),
                                   name='Actual',
                                   mode='lines',
                                   line=dict(color='#484b4e', width=2),
                                   hovertemplate='<b>Actual</b><br>Date: %{x|%d %B %Y}<br>'
                                   'Value: %{y}<extra></extra>'))

            # generate the figure
            fig = go.Figure(data=data, layout=layout).to_dict()

            model_predictions = dcc.Graph(figure=fig, config={'responsive': True, 'autosizable': True},
            style={'width': '37vw', 'height': '20vw'})

            ############################################################################################################
            # generate the ACF plot of the residuals
            ############################################################################################################

            # calculate the autocorrelations
            nlags = results.shape[0] - 1
            rho = acf(x=results['Residual'], nlags=nlags, alpha=0.05, fft=False)

            # generate the layout
            layout = dict(font=dict(family='Open Sans', color='#737373'),
                          paper_bgcolor='white',
                          plot_bgcolor='white',
                          margin=dict(l=2, r=2, t=40, b=2, pad=0),
                          xaxis=dict(title='Lag', showgrid=False, zeroline=False, mirror=True, color='#737373',
                          linecolor='#d9d9d9'),
                          yaxis=dict(title='ACF', showgrid=False, zeroline=False, mirror=True, color='#737373',
                          linecolor='#d9d9d9', range=[-1.1, 1.5], tickmode='array', tickvals=[-1, -0.5, 0, 0.5, 1]),
                          legend=dict(x=0, y=1.2, orientation='h', traceorder='reversed'))

            # generate the traces
            data = []

            data.append(go.Scatter(x=list(range(nlags)),
                                   y=list(rho[1][:, 0] - rho[0]),
                                   fillcolor='rgba(213,200,97,0.2)',
                                   fill='tozeroy',
                                   mode='lines',
                                   showlegend=False,
                                   legendgroup='group',
                                   line=dict(color='#D5C861', shape='spline', width=0.5),
                                   hovertemplate='<b>95% Confidence Interval, Lower Bound</b><br>'
                                   'Lag: %{x}<br>Value: %{y: .2f}<extra></extra>'))

            data.append(go.Scatter(x=list(range(nlags)),
                                   y=list(rho[1][:, 1] - rho[0]),
                                   fillcolor='rgba(213,200,97,0.2)',
                                   fill='tozeroy',
                                   mode='lines',
                                   legendgroup='group',
                                   line=dict(color='#D5C861', shape='spline', width=0.5),
                                   name='95% Confidence Interval',
                                   hovertemplate='<b>95% Confidence Interval, Upper Bound</b><br>'
                                   'Lag: %{x}<br>Value: %{y: .2f}<extra></extra>'))

            data.append(go.Scatter(x=list(range(nlags)),
                                   y=[0] * nlags,
                                   mode='lines',
                                   showlegend=False,
                                   hoverinfo='none',
                                   line=dict(color='#9d9d9d', width=0.5)))

            data.append(go.Scatter(x=list(range(nlags)),
                                   y=list(rho[0]),
                                   mode='markers',
                                   name='ACF of Residuals',
                                   text=list(rho[0]),
                                   marker=dict(color='#484b4e', size=4),
                                   error_y=dict(type='data', width=0.5, symmetric=False, array=[0] * nlags,
                                   arrayminus=list(rho[0])),
                                   hovertemplate='<b>ACF</b><br>Lag: %{x}<br>Value: %{text: .2f}<extra></extra>'))

            # generate the figure
            fig = go.Figure(data=data, layout=layout)

            residuals_acf = dcc.Graph(figure=fig, config={'responsive': True, 'autosizable': True},
            style={'width': '37vw', 'height': '20vw'})

            ############################################################################################################
            # generate the probability plot of the residuals
            ############################################################################################################

            # calculate the quantiles
            probplot = stats.probplot(results['Residual'])

            # extract the corresponding dates
            dates = results.sort_values(by='Residual').index

            # generate the layout
            layout = go.Layout(font=dict(family='Open Sans', color='#737373'),
                               paper_bgcolor='white',
                               plot_bgcolor='white',
                               margin=dict(l=2, r=2, t=40, b=2, pad=0),
                               xaxis=dict(title='Theoretical Quantiles', showgrid=False, zeroline=False, mirror=True,
                               color='#737373', linecolor='#d9d9d9'),
                               yaxis=dict(title='Ordered Values', showgrid=False, zeroline=False, mirror=True,
                               color='#737373', linecolor='#d9d9d9'),
                               legend=dict(x=0, y=1.2, orientation='h'))

            # generate the traces
            data = []

            data.append(go.Scatter(x=probplot[0][0],
                                   y=probplot[0][1],
                                   text=dates,
                                   mode='markers',
                                   marker=dict(color='#484b4e', size=4),
                                   name='Residual',
                                   hovertemplate='<b>Residual</b><br>Date: %{text|%d %B %Y}<br>'
                                   'Value: %{y}<extra></extra>'))

            data.append(go.Scatter(x=probplot[0][0],
                                   y=probplot[1][1] + probplot[1][0] * probplot[0][0],
                                   mode='lines',
                                   line=dict(color='rgba(127, 25, 0, 0.75)', width=2),
                                   name='Best-fit line',
                                   hoverinfo='none'))

            # generate the figure
            fig = go.Figure(data=data, layout=layout)

            residuals_probplot = dcc.Graph(figure=fig, config={'responsive': True, 'autosizable': True},
            style={'width': '37vw', 'height': '20vw'})

        except:

            model_estimation = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
            'pre', 'color': '#8d8d8d'}), 'Convergence error.'], style={'color': '#9d9d9d', 'font-size': '90%',
            'position': 'relative', 'top': '8.25vw', 'left': '15vw'})

            model_predictions = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
            'pre', 'color': '#8d8d8d'}), 'Convergence error.'], style={'color': '#9d9d9d', 'font-size': '90%',
            'position': 'relative', 'top': '8.25vw', 'left': '15vw'})

            residuals_acf = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
            'pre', 'color': '#8d8d8d'}), 'Convergence error.'], style={'color': '#9d9d9d', 'font-size': '90%',
            'position': 'relative', 'top': '8.25vw', 'left': '15vw'})

            residuals_probplot = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space':
            'pre', 'color': '#8d8d8d'}), 'Convergence error.'], style={'color': '#9d9d9d', 'font-size': '90%',
            'position': 'relative', 'top': '8.25vw', 'left': '15vw'})

            return [model_estimation, model_predictions, residuals_acf, residuals_probplot]

    else:

        model_estimation = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space': 'pre',
        'color': '#8d8d8d'}), 'No data.'], style={'color': '#9d9d9d', 'font-size': '90%', 'position': 'relative',
        'top': '8.25vw', 'left': '15vw'})

        model_predictions = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space': 'pre',
        'color': '#8d8d8d'}), 'No data.'], style={'color': '#9d9d9d', 'font-size': '90%', 'position': 'relative',
        'top': '8.25vw', 'left': '15vw'})

        residuals_acf = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space': 'pre',
        'color': '#8d8d8d'}), 'No data.'], style={'color': '#9d9d9d', 'font-size': '90%', 'position': 'relative',
        'top': '8.25vw', 'left': '15vw'})

        residuals_probplot = html.Div(children=[html.Span('Error: ', style={'font-weight': 'bold', 'white-space': 'pre',
        'color': '#8d8d8d'}), 'No data.'], style={'color': '#9d9d9d', 'font-size': '90%', 'position': 'relative',
        'top': '8.25vw', 'left': '15vw'})

    return [model_estimation, model_predictions, residuals_acf, residuals_probplot]

if __name__ == '__main__':
    app.run_server(debug=False)