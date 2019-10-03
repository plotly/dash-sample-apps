# -*- coding: utf-8 -*-
import flask
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input, State
import pandas as pd
from datetime import datetime
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from io import BytesIO
import base64
from os import path, getcwd
from PIL import Image
from wordcloud import WordCloud, ImageColorGenerator, STOPWORDS



external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
filename = "customer_complaints_narrative_sample.csv"

"""
#  Helpful functions
"""

def load_data(filename, n):
    print("loading data...")
    df = pd.read_csv(filename, header=0, skiprows=lambda i: i % n != 0)
    print("loading data...done")
    return df

def sample_data(dataframe, float_percent):
    return dataframe.sample(frac=float_percent, random_state=1)

def crunch_data(dataframe):
    companyCounts = dataframe['Company'].value_counts()
    values = companyCounts.keys().tolist()
    counts = companyCounts.tolist()
    companyCounts_sample = companyCounts[:sample_size]
    values_sample = companyCounts_sample.keys().tolist()
    counts_sample = companyCounts_sample.tolist()

def calculate_full(dataframe):
    companyCounts = dataframe['Company'].value_counts()
    values = companyCounts.keys().tolist()
    counts = companyCounts.tolist()
    return values, counts

def calculate_sample(dataframe, sample_size, time_values):
    print("got:", str(time_values))
    dataframe['Date received'] =  pd.to_datetime(dataframe['Date received'], format='%m/%d/%Y')
    if time_values is not None:
        dataframe = dataframe[(dataframe['Date received'] >= datetime.strptime(str(time_values[0]), '%Y')) & (dataframe['Date received'] <= datetime.strptime(str(time_values[1])+"/12/31", '%Y/%m/%d'))]
    companyCounts = dataframe['Company'].value_counts()
    companyCounts_sample = companyCounts[:sample_size]
    values_sample = companyCounts_sample.keys().tolist()
    counts_sample = companyCounts_sample.tolist()
    all_complaints = dataframe['Consumer complaint narrative'].dropna().values
    return values_sample, counts_sample, all_complaints

def calculate_dates(dataframe):
    date_data = dataframe['Date received']
    date_data = date_data.apply(lambda x : x[-4:])
    max_date = date_data.max() # who doesn't like a good magic number to grab the year from a string?
    min_date = date_data.min() # - " -
    print(min_date, max_date)
    return max_date, min_date

def calculate_n(inp):
    return int(100 / inp)

def make_marks(mini, maxi):
    mini = int(mini)
    maxi = int(maxi)
    ret = {}
    i = 0
    while mini+i <= maxi:
        ret[mini+i] = str(mini+i)
        i+=1
    return ret

def make_options(values):
    ret = []
    for value in values:
        ret.append({'label': value, 'value': value})
    return ret

def make_n_marks():
    # TODO: Johan. there must be another way, but for now this will do.
    ret = {}
    i = 0
    while i <= 100:
        if i % 10 == 0:
            ret[i]=str(i)+"%"
        i+=1
    return ret

def fig_to_uri(in_fig, close_all=True, **save_args):
    # type: (plt.Figure) -> str
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)

def create_wordcloud(df):
    complaints_text  = list(df["Consumer complaint narrative"].dropna().values)

    ## join all documents in corpus
    text = " ".join(list(complaints_text))
    print("Complaints received")
    print(len(complaints_text))

    d = getcwd()
    mask = np.array(Image.open(path.join(d, "thumbs-down.png")))

    STOPWORDS.add('XXXX')
    STOPWORDS.add('XX')
    STOPWORDS.add('xx')
    STOPWORDS.add('xxxx')
    # TODO exclude name of all banks here
    STOPWORDS.add('wells')
    STOPWORDS.add('fargo')

    wc = WordCloud(background_color="white", stopwords=STOPWORDS, max_words=1000, mask=mask,
                max_font_size=90, random_state=42, contour_width=1, contour_color='#119DFF')
    wc.generate(text)

    # create wordcloud shape from image
    fig = plt.figure(figsize=[8,8])
    ax = plt.imshow(wc.recolor(), interpolation="bilinear")
    plt.axis("off")
    out_url = fig_to_uri(fig, bbox_inches='tight')
    return out_url

"""
#  Page layout and contents
"""

global_df = load_data(filename, 1) 

left_column = [html.Label('Please select percentage of dataset (higher is more accurate but also slower)'),
               dcc.Slider(id="n-selection-slider",
                          min=1,
                          max=100,
                          step=1,
                          marks=make_n_marks(),
                          value=10),
               html.Label("_"),
               html.Label('Please select a bank, using either the dropdown or the plot below.'),
               dcc.Dropdown(id="bank_drop"),
               html.Div(dcc.Graph(id='bank-sample')),
               html.Label("Please select time frame"),
               html.Div(dcc.RangeSlider(id='time-window-slider')),
              ]

right_column = [html.Label('RIGHT COLUMN TOP'),
                html.Div(id='output-container-range-slider'),
                html.Div(id='debug-1'),
                html.Div([html.Img(id = 'bank-wordcloud', src = '')], id='wordcloud_div'),
                html.Div(id='debug-2'),
                html.Div(id='debug-3'),
               ]

server = flask.Flask(__name__)
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, server=server)
app.layout = html.Div(className='container', children=[
                html.Div(className='', children=[
                    html.Div(className='five columns', children=left_column),
                    html.Div(className='seven columns', children=right_column),
                    ])
                ]
            )

"""
#  Callbacks
"""

@app.callback(
    [Output('time-window-slider', 'marks'),
     Output('time-window-slider', 'min'),
     Output('time-window-slider', 'max'),
     Output('time-window-slider', 'step'),
     Output('time-window-slider', 'value')],
     [Input('n-selection-slider', 'value')])
def populate_time_slider(value):
    # print("repopulating time-window-slider")
    max_date, min_date = calculate_dates(global_df)
    print(max_date, min_date)
    #print(make_marks(min_date, max_date), int(min_date), max(min_date), 1, [int(min_date), int(max_date)])
    # {2012: '2012', 2013: '2013', 2014: '2014', 2015: '2015', 2016: '2016', 2017: '2017'} 2012 2 1 [2012, 2017]
    return make_marks(min_date, max_date), int(min_date), int(max_date), 1, [int(min_date), int(max_date)]

@app.callback(Output('bank_drop', 'options'), [Input('time-window-slider', 'value')])
def populate_bank_dropdown(time_values):
    print("bank_drop: TODO USE THE TIME VALUES TO LIMIT THE DATASET")
    values, counts = calculate_full(global_df)
    # print("repopulating dropdown")
    return make_options(values)

@app.callback(Output('bank-sample', 'figure'), [Input('n-selection-slider', 'value'), Input('time-window-slider', 'value')])
def set_n(n_value, time_values):
    n = float(n_value/100)
    print("redrawing bank-sample...")
    # print("n set to: ", n)
    # print("time_values: ", time_values)
    # print("bank-sample: TODO USE THE TIME VALUES TO LIMIT THE DATASET")
    sample_size = 20
    local_df = sample_data(global_df, n)
    values_sample, counts_sample, texts = calculate_sample(local_df, sample_size, time_values)
    print("redrawing bank-sample...done")
    return {'data': [{'x': values_sample, 'y': counts_sample, 'type': 'bar', 'name': ''}],'layout': {'title': 'Top 20 offenders'}}

@app.callback(Output('output-container-range-slider', 'children'), [Input('time-window-slider', 'value')])
def update_output(value):
    # print(value)
    return 'You have selected time frame: "{}"'.format(value)

@app.callback(
[Output('debug-1', 'children'),
Output('bank-wordcloud', 'src')],
[Input("bank-sample", "clickData"),
Input('bank_drop', 'value'),
Input('time-window-slider', 'value'),
Input('n-selection-slider', 'value')])
def update_wordcloud(value_click, value_drop, time_values, n_selection):
    if value_drop:
        selected_bank = value_drop
    elif value_click: 
        selected_bank = value_click['points'][0]['x']
    else:
        return ["",""]
    print("redrawing bank-wordcloud...")
    n = float(n_selection/100)
    print("got time window:", str(time_values))
    print("got n_selection:", str(n_selection), str(n))
    # df = global_df
    # sample the dataset according to the slider
    
    local_df = sample_data(global_df, n)
    if time_values is not None:
        local_df['Date received'] =  pd.to_datetime(local_df['Date received'], format='%m/%d/%Y')
        local_df = local_df[(local_df['Date received'] >= datetime.strptime(str(time_values[0]), '%Y')) & (local_df['Date received'] <= datetime.strptime(str(time_values[1])+"/12/31", '%Y/%m/%d'))]
    local_df = local_df[local_df["Company"]==selected_bank]
    wordcloud_url = create_wordcloud(local_df)
    print("redrawing bank-wordcloud...done")
    return([update_debug(selected_bank, 'graph'), wordcloud_url])


@app.callback(Output('bank_drop', 'value'), [Input("bank-sample", "clickData")])
def update_bank_click(value):
    if value is not None:
        selected_bank = value['points'][0]['x']
        return(selected_bank)

def update_debug(input_value, source):
    if input_value is not None:
        return "picked value %s from %s" % (input_value, source)


if __name__ == '__main__':
    app.run_server(debug=True)