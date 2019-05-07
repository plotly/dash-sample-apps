import numpy as np
import datetime as dt
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

from dash.dependencies import Input, Output, State
from scipy.stats import rayleigh
from data.api import get_wind_data, get_wind_data_by_id


external_css = [
    "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
    "https://fonts.googleapis.com/css?family=Raleway:400,400i,700,700i",
    "https://fonts.googleapis.com/css?family=Product+Sans:400,400i,700,700i"
]

app = dash.Dash('streaming-wind-app', external_stylesheets=external_css)

app.layout = html.Div([
    html.Div([
        html.H2("Wind Speed Streaming"),
        html.Img(src="https://s3-us-west-1.amazonaws.com/plotly-tutorials/logo/new-branding/dash-logo-by-plotly-stripe-inverted.png"),
    ], className='banner'),
    html.Div([
        html.Div([
            html.H3("WIND SPEED (mph)")
        ], className='Title'),
        html.Div([
            dcc.Graph(id='wind-speed'),
        ], className='twelve columns wind-speed'),
        dcc.Interval(id='wind-speed-update', interval=1000, n_intervals=0),
    ], className='row wind-speed-row'),
    html.Div([
        html.Div([
            html.Div([
                html.H3("WIND SPEED HISTOGRAM")
            ], className='Title'),
            html.Div([
                dcc.Slider(
                    id='bin-slider',
                    min=1,
                    max=60,
                    step=1,
                    value=20,
                    updatemode='drag'
                ),
            ], className='histogram-slider'),
            html.P('# of Bins: Auto', id='bin-size', className='bin-size'),
            html.Div([
                dcc.Checklist(
                    id='bin-auto',
                    options=[
                        {'label': 'Auto', 'value': 'Auto'}
                    ],
                    values=['Auto']
                ),
            ], className='bin-auto'),
            dcc.Graph(id='wind-histogram'),
        ], className='seven columns wind-histogram'),
        html.Div([
            html.Div([
                html.H3("WIND DIRECTION")
            ], className='Title'),
            dcc.Graph(id='wind-direction'),
        ], className='five columns wind-polar')
    ], className='row wind-histo-polar')
], style={'padding': '0px 10px 15px 10px',
          'marginLeft': 'auto', 'marginRight': 'auto', "width": "900px",
          'boxShadow': '0px 0px 5px 5px rgba(204,204,204,0.4)'})


def get_current_time():
    """ Helper function to get the current time in seconds. """

    now = dt.datetime.now()
    total_time = (now.hour * 3600) + (now.minute * 60) + (now.second)
    return total_time


@app.callback(
    Output('wind-speed', 'figure'),
    [Input('wind-speed-update', 'n_intervals')]
)
def gen_wind_speed(interval):
    """
    Generate the wind speed graph.

    :params interval: update the graph based on an interval
    """

    total_time = get_current_time()
    df = get_wind_data(total_time-200, total_time)

    trace = go.Scatter(
        y=df['Speed'],
        line={'color': '#42C4F7'},
        hoverinfo='skip',
        error_y={
            'type': 'data',
            'array': df['SpeedError'],
            'thickness': 1.5,
            'width': 2,
            'color': '#B4E8FC'
        },
        mode='lines'
    )

    layout = go.Layout(
        height=450,
        xaxis={
            'range': [0, 200],
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'fixedrange': True,
            'tickvals': [0, 50, 100, 150, 200],
            'ticktext': ['200', '150', '100', '50', '0'],
            'title': 'Time Elapsed (sec)'
        },
        yaxis={
            'range': [min(0, min(df['Speed'])), max(45, max(df['Speed'])+max(df['SpeedError']))],
            'showline': False,
            'fixedrange': True,
            'zeroline': False,
            'nticks': max(6, round(df['Speed'].iloc[-1]/10))
        },
        margin={'t': 45, 'l': 50, 'r': 50}
    )

    return go.Figure(data=[trace], layout=layout)


@app.callback(
    Output('wind-direction', 'figure'),
    [Input('wind-speed-update', 'n_intervals')]
)
def gen_wind_direction(interval):
    """
    Generate the wind direction graph.

    :params interval: update the graph based on an interval
    """

    total_time = get_current_time()
    df = get_wind_data_by_id(total_time)
    val = df['Speed'].iloc[-1]
    direction = [0, (df['Direction'][0]-20), (df['Direction'][0]+20), 0]

    traces_scatterpolar = [
        {
            'r': [0, val, val, 0],
            'fillcolor': 'rgb(242, 196, 247)'
        },
        {
            'r': [0, val*0.65, val*0.65, 0],
            'fillcolor': 'rgb(242, 196, 247)'
        },
        {
            'r': [0, val*0.3, val*0.3, 0],
            'fillcolor': '#FAEBFC'
        }
    ]

    data = [
        go.Scatterpolar(
            r=traces['r'],
            theta=direction,
            mode='lines',
            fill='toself',
            fillcolor=traces['fillcolor'],
            line={
                'color': 'rgba(32, 32, 32, .6)',
                'width': 1
            }
        )
        for traces in traces_scatterpolar
    ]

    layout = go.Layout(
        autosize=True,
        width=275,
        margin={'t': 10, 'b': 10, 'r': 30, 'l': 40},
        polar={
            'bgcolor': '#F2F2F2',
            'radialaxis': {
                'range': [0, 45],
                'angle': 45,
                'dtick': 10
            },
            'angularaxis': {
                'showline': False,
                'tickcolor': 'white',
            }
        },
        showlegend=False,
    )

    return go.Figure(data=data, layout=layout)


@app.callback(
    Output('wind-histogram', 'figure'),
    [Input('wind-speed-update', 'n_intervals')],
    [
        State('wind-speed', 'figure'),
        State('bin-slider', 'value'),
        State('bin-auto', 'values')
    ]
)
def gen_wind_histogram(interval, wind_speed_figure, slider_value, auto_state):
    """
    Genererate wind histogram graph.

    :params interval: upadte the graph based on an interval
    :params wind_speed_figure: current wind speed graph
    :params slider_value: current slider value
    :params auto_state: current auto state
    """

    wind_val = []

    # Check to see whether wind-speed has been plotted yet
    if wind_speed_figure is not None:
        wind_val = wind_speed_figure['data'][0]['y']
    if 'Auto' in auto_state:
        bin_val = np.histogram(wind_val, bins=range(int(round(min(wind_val))),
                                                    int(round(max(wind_val)))))
    else:
        bin_val = np.histogram(wind_val, bins=slider_value)

    avg_val = float(sum(wind_val))/len(wind_val)
    median_val = np.median(wind_val)

    pdf_fitted = rayleigh.pdf(bin_val[1], loc=(avg_val)*0.55,
                              scale=(bin_val[1][-1] - bin_val[1][0])/3)

    y_val = pdf_fitted * max(bin_val[0]) * 20,
    y_val_max = max(y_val[0])
    bin_val_max = max(bin_val[0])

    trace = go.Bar(
        x=bin_val[1],
        y=bin_val[0],
        marker={
            'color': '#7F7F7F'
        },
        showlegend=False,
        hoverinfo='x+y'
    )

    traces_scatter = [
        {
            'line_dash': 'dash',
            'line_color': '#2E5266',
            'name': 'Average'
        },
        {
            'line_dash': 'dot',
            'line_color': '#BD9391',
            'name': 'Median'
        }
    ]

    scatter_data = [
        go.Scatter(
            x=[bin_val[int(len(bin_val)/2)]],
            y=[0],
            mode='lines',
            line={
                'dash': traces['line_dash'],
                'color': traces['line_color']
            },
            marker={
                'opacity': 0,
            },
            visible=True,
            name= traces['name']
        )
        for traces in traces_scatter
    ]
    
    trace3 = go.Scatter(
        mode='lines',
        line={
            'color': '#42C4F7'
        },
        y=y_val[0],
        x=bin_val[1][:len(bin_val[1])],
        name='Rayleigh Fit'
    )
    layout = go.Layout(
        xaxis={
            'title': 'Wind Speed (mph)',
            'showgrid': False,
            'showline': False,
            'fixedrange': True
        },
        yaxis={
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'title': 'Number of Samples',
            'fixedrange': True
        },
        margin={'t': 50, 'b': 20, 'r': 50},
        autosize=True,
        bargap=0.01,
        bargroupgap=0,
        hovermode='closest',
        legend={'x': 0.175, 'y': -0.2, 'orientation': 'h'},
        shapes=[
            {
                'xref': 'x',
                'yref': 'y',
                'y1': int(max(bin_val_max, y_val_max))+0.5,
                'y0': 0,
                'x0': avg_val,
                'x1': avg_val,
                'type': 'line',
                'line': {
                    'dash': 'dash',
                    'color': '#2E5266',
                    'width': 5
                }
            },
            {
                'xref': 'x',
                'yref': 'y',
                'y1': int(max(bin_val_max, y_val_max))+0.5,
                'y0': 0,
                'x0': median_val,
                'x1': median_val,
                'type': 'line',
                'line': {
                    'dash': 'dot',
                    'color': '#BD9391',
                    'width': 5
                }
            }
        ]
    )
    return go.Figure(data=[trace, scatter_data[0], scatter_data[1], trace3], layout=layout)


@app.callback(
    Output('bin-auto', 'values'),
    [Input('bin-slider', 'value')],
    [State('wind-speed', 'figure')]
)
def deselect_auto(slider_value, wind_speed_figure):
    """ Toggle the auto checkbox. """

    if (wind_speed_figure is not None and
            len(wind_speed_figure['data'][0]['y']) > 5):
        return ['']
    else:
        return ['Auto']


@app.callback(
    Output('bin-size', 'children'),
    [Input('bin-auto', 'values')],
    [State('bin-slider', 'value')]
)
def show_num_bins(autoValue, slider_value):
    """ Display the number of bins. """

    if 'Auto' in autoValue:
        return '# of Bins: Auto'
    else:
        return '# of Bins: ' + str(int(slider_value))


if __name__ == '__main__':
    app.run_server()
