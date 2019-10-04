# -*- coding: utf-8 -*-
import flask
import dash
import dash_table
import matplotlib.colors as mcolors
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from dash.dependencies import Output, Input, State
from datetime import datetime
import pandas as pd
import numpy as np
from wordcloud import WordCloud, STOPWORDS
import re
from ldacomplaints import lda_analysis

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
    
    return values_sample, counts_sample 

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

def add_stopwords(selected_bank):
    STOPWORDS.add('XXXX')
    STOPWORDS.add('XX')
    STOPWORDS.add('xx')
    STOPWORDS.add('xxxx')
    selected_bank_words = re.findall(r"[\w']+", selected_bank)
    for word in selected_bank_words: 
        #STOPWORDS.add(word)
        STOPWORDS.add(word.lower())

    print("ADD STOPWORDS")
    print(selected_bank_words)
    return STOPWORDS

def populate_lda_scatter(tsne_lda, lda_model, topic_num, df_dominant_topic):
    topic_top3words = [(i, topic) for i, topics in lda_model.show_topics(formatted=False) 
                                 for j, (topic, wt) in enumerate(topics) if j < 3]

    df_top3words_stacked = pd.DataFrame(topic_top3words, columns=['topic_id', 'words'])
    df_top3words = df_top3words_stacked.groupby('topic_id').agg(', \n'.join)
    df_top3words.reset_index(level=0,inplace=True)
  
    tsne_df = pd.DataFrame({'tsne_x': tsne_lda[:,0], 'tsne_y':  tsne_lda[:,1], 'topic_num': topic_num, 'doc_num': df_dominant_topic["Document_No"]})
    mycolors = np.array([color for name, color in mcolors.TABLEAU_COLORS.items()])

    # Plot and embed in ipython notebook!
    # for each topic create separate trace
    traces=[]
    for topic_id in df_top3words['topic_id']:
        #print('Topic: {} \nWords: {}'.format(idx, topic))
        tsne_df_f = tsne_df[tsne_df.topic_num==topic_id]
        cluster_name = ', '.join(df_top3words[df_top3words['topic_id']==topic_id]['words'].to_list())
        trace = go.Scatter(name = cluster_name, x=tsne_df_f['tsne_x'], y=tsne_df_f['tsne_y'], mode = 'markers', hovertext= tsne_df_f['doc_num'],
                marker=dict(
                    size=6,
                    color = mycolors[tsne_df_f['topic_num']], #set color equal to a variable
                    colorscale='Viridis',
                    showscale=False
            )
        )
        traces.append(trace)

    layout = go.Layout({'title': 'Topic analysis using LDA'})

    return {'data': traces,'layout': layout}

def plotly_wordcloud(df):
    complaints_text  = list(df["Consumer complaint narrative"].dropna().values)

    ## join all documents in corpus
    text = " ".join(list(complaints_text))

    wc = WordCloud(stopwords = set(STOPWORDS),
                   max_words = 100,
                   max_font_size = 90)
    wc.generate(text)
    
    word_list=[]
    freq_list=[]
    fontsize_list=[]
    position_list=[]
    orientation_list=[]
    color_list=[]

    for (word, freq), fontsize, position, orientation, color in wc.layout_:
        word_list.append(word)
        freq_list.append(freq)
        fontsize_list.append(fontsize)
        position_list.append(position)
        orientation_list.append(orientation)
        color_list.append(color)
        
    # get the positions
    x=[]
    y=[]
    for i in position_list:
        x.append(i[0])
        y.append(i[1]) 
            
    # get the relative occurence frequencies
    new_freq_list = []
    for i in freq_list:
        new_freq_list.append(i*80)
    
    trace = go.Scatter(x=x, 
                       y=y, 
                       textfont = dict(size=new_freq_list,
                                       color=color_list),
                       hoverinfo='text',
                       textposition='top center',
                       hovertext=['{0} - {1}'.format(w, f) for w, f in zip(word_list, freq_list)],
                       mode='text',  
                       text=word_list
                      )
    
    layout = go.Layout({'xaxis': {'showgrid': False, 'showticklabels': False, 'zeroline': False, 'automargin': True, 'range': [-100, 250]},
                        'yaxis': {'showgrid': False, 'showticklabels': False, 'zeroline': False, 'automargin': True, 'range': [-100, 450]},
                        'margin': dict(t=100, b=100, l=100, r=100, pad=4),
                        'title': 'Most frequent words in a complaint'})
    
    return {'data': [trace],'layout': layout}

"""
#  Page layout and contents
"""

global_df = load_data(filename, 1) 

left_column = [ html.H3(children='Select bank & dataset size'),
                html.Label('Select percentage of dataset (higher is more accurate but also slower)',
                        style={'marginTop': 30}),
               dcc.Slider(id="n-selection-slider",
                          min=1,
                          max=100,
                          step=1,
                          marks=make_n_marks(),
                          value=10),
               html.Label('Select a bank, using either the dropdown or the plot below.', 
                          style={'marginTop': 50}),
               dcc.Dropdown(id='bank-drop', clearable=False, style={'marginBottom': 50}),
               html.Label("Select time frame"),
               html.Div(dcc.RangeSlider(id='time-window-slider'), style={'marginBottom': 50})
              ]

right_column = [html.Label('RIGHT COLUMN TOP', style={'display': 'none'}),
                html.Div(id='output-container-range-slider', style={'display': 'none'}),
                html.Div(id='debug-1', style={'display': 'none'}),
                html.Div(dcc.Graph(id='bank-sample')),
                html.Div(id='debug-2', style={'display': 'none'}),
                html.Div(id='debug-3', style={'display': 'none'})
               ]

server = flask.Flask(__name__)
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, server=server)
app.layout = html.Div(className='container', children=[
                html.H1(children='Banks customers complaints'),
                html.Div(className='', children=[
                    html.Div(className='five columns', children=left_column, style={'border': '1px solid', 'background': '#f9fafe',  'padding': '20px'}),
                    html.Div(className='seven columns offset-by-one-half.column', children=right_column),
                    html.Div(className='twelve columns', children=dcc.Loading(id="loading-wordcloud", children=[dcc.Graph(id='bank-wordcloud')], type="default")),
                    html.Div(className='twelve columns', 
                        children=dcc.Loading(id="loading-lda", children=[dcc.Graph(id='tsne-lda')], type="default")),
                    html.Div(className='twelve columns', 
                        children=dcc.Loading(id="loading-table", children=[dash_table.DataTable(
                                                                    id='lda-table',
                                                                    style_cell_conditional=[
                                                                        {
                                                                            'if': {'column_id': 'Text'},
                                                                            'textAlign': 'left',
                                                                        }
                                                                    ],
                                                                    style_cell={'padding': '5px'},
                                                                    style_header={
                                                                        'backgroundColor': 'white',
                                                                         'fontWeight': 'bold'
                                                                    },
                                                                    style_data={
                                                                        'whiteSpace': 'normal',
                                                                        'height': 'auto'
                                                                    },
                                                                    filter_action="native",
                                                                    page_action='native',
                                                                    page_current= 0,
                                                                    page_size= 5,
                                                                    columns=[],
                                                                    data=[])], type="default"))
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

@app.callback(Output('bank-drop', 'options'), [Input('time-window-slider', 'value')])
def populate_bank_dropdown(time_values):
    print("bank-drop: TODO USE THE TIME VALUES TO LIMIT THE DATASET")
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
    values_sample, counts_sample = calculate_sample(local_df, sample_size, time_values)
    print("redrawing bank-sample...done")
    return {'data': [{'x': values_sample, 'y': counts_sample, 'type': 'bar', 'name': ''}],'layout': {'title': 'Top 20 offenders'}}

@app.callback(Output('output-container-range-slider', 'children'), [Input('time-window-slider', 'value')])
def update_output(value):
    # print(value)
    return 'You have selected time frame: "{}"'.format(value)

@app.callback([Output('lda-table', 'data'),
Output('lda-table', 'columns'),
Output('tsne-lda', 'figure')],
[Input("bank-sample", "clickData"),
Input('bank-drop', 'value'),
Input('time-window-slider', 'value'),
Input('n-selection-slider', 'value')])
def update_lda_table(bank_click, value_drop, time_values, n_selection):
    if value_drop:
        selected_bank = value_drop
        print(selected_bank)
    elif bank_click: 
        selected_bank = bank_click['points'][0]['x']
    else:
        return [[],[], {}]

    print("redrawing lda table...")
    n = float(n_selection/100)
    print("got time window:", str(time_values))
    print("got n_selection:", str(n_selection), str(n))


    # sample the dataset according to the slider
    local_df = sample_data(global_df, n)
    if time_values is not None:
        local_df['Date received'] =  pd.to_datetime(local_df['Date received'], format='%m/%d/%Y')
        local_df = local_df[(local_df['Date received'] >= datetime.strptime(str(time_values[0]), '%Y')) & (local_df['Date received'] <= datetime.strptime(str(time_values[1])+"/12/31", '%Y/%m/%d'))]
    local_df = local_df[local_df["Company"]==selected_bank]
    
    add_stopwords(selected_bank)

    complaints_text  = list(local_df["Consumer complaint narrative"].dropna().values)
    tsne_lda, lda_model, topic_num,  df_dominant_topic = lda_analysis(complaints_text, list(STOPWORDS))
    
    lda_scatter_figure = populate_lda_scatter(tsne_lda, lda_model, topic_num, df_dominant_topic)
    
    columns = [{'name': i, 'id': i} for i in df_dominant_topic.columns]
    data = df_dominant_topic.to_dict('records')
    
    return (data, columns, lda_scatter_figure)

@app.callback(
[Output('debug-1', 'children'),
Output('bank-wordcloud', 'figure')],
[Input("bank-sample", "clickData"),
Input('bank-drop', 'value'),
Input('time-window-slider', 'value'),
Input('n-selection-slider', 'value')])
def update_wordcloud(value_click, value_drop, time_values, n_selection):
    if value_drop:
        selected_bank = value_drop
    elif value_click: 
        selected_bank = value_click['points'][0]['x']
    else:
        return ["",{}]
    print("redrawing bank-wordcloud...")
    n = float(n_selection/100)
    print("got time window:", str(time_values))
    print("got n_selection:", str(n_selection), str(n))

    # sample the dataset according to the slider
    local_df = sample_data(global_df, n)
    if time_values is not None:
        local_df['Date received'] =  pd.to_datetime(local_df['Date received'], format='%m/%d/%Y')
        local_df = local_df[(local_df['Date received'] >= datetime.strptime(str(time_values[0]), '%Y')) & (local_df['Date received'] <= datetime.strptime(str(time_values[1])+"/12/31", '%Y/%m/%d'))]
    local_df = local_df[local_df["Company"]==selected_bank]
    
    add_stopwords(selected_bank)
    wordcloud = plotly_wordcloud(local_df)
    
    print("redrawing bank-wordcloud...done")
    return([update_debug(selected_bank, 'graph'), wordcloud])


@app.callback([Output('lda-table', 'filter_query'), Output('lda-table', 'style_data_conditional')], [Input("tsne-lda", "clickData")])
def filter_table_on_scatter_click(tsne_click):
    if tsne_click is not None:
        selected_complaint = tsne_click['points'][0]['hovertext']
        filter_query = '{Document_No} eq ' + str(selected_complaint)
        print(filter_query)
        styled_data = [
                        {
                        'if': {'filter_query': filter_query},
                        "backgroundColor": "#3D9970",
                        'color': 'white'
                        }
                    ]
        return(filter_query, styled_data)
    else:
        return ['', []]

@app.callback(Output('bank-drop', 'value'), [Input("bank-sample", "clickData")])
def update_bank_click(value):
    if value is not None:
        selected_bank = value['points'][0]['x']
        return(selected_bank)
    else:
        return('CITIBANK, N.A.')

def update_debug(input_value, source):
    if input_value is not None:
        return "picked value %s from %s" % (input_value, source)


if __name__ == '__main__':
    app.run_server(debug=False)
