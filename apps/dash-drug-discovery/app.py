import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_core_components as dcc
import pandas as pd
import flask
from flask_cors import CORS
import os

app = dash.Dash('drug-discovery')
server = app.server

df = pd.read_csv('small_molecule_drugbank.csv').drop(['Unnamed: 0'],axis=1)

if 'DYNO' in os.environ:
    app.scripts.append_script({
        'external_url': 'https://cdn.rawgit.com/chriddyp/ca0d8f02a1659981a0ea7f013a378bbd/raw/e79f3f789517deec58f41251f7dbb6bee72c44ab/plotly_ga.js'
    })

def add_markers( figure_data, molecules, plot_type = 'scatter3d' ):
    indices = []
    drug_data = figure_data[0]
    for m in molecules:
        hover_text = drug_data['text']
        for i in range(len(hover_text)):
            if m == hover_text[i]:
                indices.append(i)

    if plot_type == 'histogram2d':
        plot_type = 'scatter'

    traces = []
    for point_number in indices:
        trace = dict(
            x = [ drug_data['x'][point_number] ],
            y = [ drug_data['y'][point_number] ],
            marker = dict(
                color = 'red',
                size = 16,
                opacity = 0.6,
                symbol = 'cross'
            ),
            type = plot_type
        )

        if plot_type == 'scatter3d':
            trace['z'] = [ drug_data['z'][point_number] ]

        traces.append(trace)

    return traces

BACKGROUND = 'rgb(230, 230, 230)'

COLORSCALE = [ [0, "rgb(244,236,21)"], [0.3, "rgb(249,210,41)"], [0.4, "rgb(134,191,118)"],
                [0.5, "rgb(37,180,167)"], [0.65, "rgb(17,123,215)"], [1, "rgb(54,50,153)"] ]

def scatter_plot_3d(
        x = df['PKA'],
        y = df['LOGP'],
        z = df['SOL'],
        size = df['MW'],
        color = df['MW'],
        xlabel = 'LogP',
        ylabel = 'pkA',
        zlabel = 'Solubility (mg/ml)',
        plot_type = 'scatter3d',
        markers = [] ):

    def axis_template_3d( title, type='linear' ):
        return dict(
            showbackground = True,
            backgroundcolor = BACKGROUND,
            gridcolor = 'rgb(255, 255, 255)',
            title = title,
            type = type,
            zerolinecolor = 'rgb(255, 255, 255)'
        )

    def axis_template_2d(title):
        return dict(
            xgap = 10, ygap = 10,
            backgroundcolor = BACKGROUND,
            gridcolor = 'rgb(255, 255, 255)',
            title = title,
            zerolinecolor = 'rgb(255, 255, 255)',
            color = '#444'
        )

    def blackout_axis( axis ):
        axis['showgrid'] = False
        axis['zeroline'] = False
        axis['color']  = 'white'
        return axis

    data = [ dict(
        x = x,
        y = y,
        z = z,
        mode = 'markers',
        marker = dict(
                colorscale = COLORSCALE,
                colorbar = dict( title = "Molecular<br>Weight" ),
                line = dict( color = '#444' ),
                reversescale = True,
                sizeref = 45,
                sizemode = 'diameter',
                opacity = 0.7,
                size = size,
                color = color,
            ),
        text = df['NAME'],
        type = plot_type,
    ) ]

    layout = dict(
        font = dict( family = 'Raleway' ),
        hovermode = 'closest',
        margin = dict( r=20, t=0, l=0, b=0 ),
        showlegend = False,
        scene = dict(
            xaxis = axis_template_3d( xlabel ),
            yaxis = axis_template_3d( ylabel ),
            zaxis = axis_template_3d( zlabel, 'log' ),
            camera = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=0.08, y=2.2, z=0.08)
            )
        )
    )

    if plot_type in ['histogram2d', 'scatter']:
        layout['xaxis'] = axis_template_2d(xlabel)
        layout['yaxis'] = axis_template_2d(ylabel)
        layout['plot_bgcolor'] = BACKGROUND
        layout['paper_bgcolor'] = BACKGROUND
        del layout['scene']
        del data[0]['z']

    if plot_type == 'histogram2d':
        # Scatter plot overlay on 2d Histogram
        data[0]['type'] = 'scatter'
        data.append( dict(
                x = x,
                y = y,
                type = 'histogram2d',
                colorscale = 'Greys',
                showscale = False
            ) )
        layout['plot_bgcolor'] = 'black'
        layout['paper_bgcolor'] = 'black'
        layout['xaxis'] = blackout_axis(layout['xaxis'])
        layout['yaxis'] = blackout_axis(layout['yaxis'])
        layout['font']['color'] = 'white'

    if len(markers) > 0:
        data = data + add_markers( data, markers, plot_type = plot_type )

    return dict( data=data, layout=layout )


BACKGROUND = 'rgb(230, 230, 230)'

COLORSCALE = [ [0, "rgb(244,236,21)"], [0.3, "rgb(249,210,41)"], [0.4, "rgb(134,191,118)"],
                [0.5, "rgb(37,180,167)"], [0.65, "rgb(17,123,215)"], [1, "rgb(54,50,153)"] ]

def scatter_plot_3d(
        x = df['PKA'],
        y = df['LOGP'],
        z = df['SOL'],
        size = df['MW'],
        color = df['MW'],
        xlabel = 'LogP',
        ylabel = 'pkA',
        zlabel = 'Solubility (mg/ml)',
        plot_type = 'scatter3d',
        markers = [] ):

    def axis_template_3d( title, type='linear' ):
        return dict(
            showbackground = True,
            backgroundcolor = BACKGROUND,
            gridcolor = 'rgb(255, 255, 255)',
            title = title,
            type = type,
            zerolinecolor = 'rgb(255, 255, 255)'
        )

    def axis_template_2d(title):
        return dict(
            xgap = 10, ygap = 10,
            backgroundcolor = BACKGROUND,
            gridcolor = 'rgb(255, 255, 255)',
            title = title,
            zerolinecolor = 'rgb(255, 255, 255)',
            color = '#444'
        )

    def blackout_axis( axis ):
        axis['showgrid'] = False
        axis['zeroline'] = False
        axis['color']  = 'white'
        return axis

    data = [ dict(
        x = x,
        y = y,
        z = z,
        mode = 'markers',
        marker = dict(
                colorscale = COLORSCALE,
                colorbar = dict( title = "Molecular<br>Weight" ),
                line = dict( color = '#444' ),
                reversescale = True,
                sizeref = 45,
                sizemode = 'diameter',
                opacity = 0.7,
                size = size,
                color = color,
            ),
        text = df['NAME'],
        type = plot_type,
    ) ]

    layout = dict(
        font = dict( family = 'Raleway' ),
        hovermode = 'closest',
        margin = dict( r=20, t=0, l=0, b=0 ),
        showlegend = False,
        scene = dict(
            xaxis = axis_template_3d( xlabel ),
            yaxis = axis_template_3d( ylabel ),
            zaxis = axis_template_3d( zlabel, 'log' ),
            camera = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=0.08, y=2.2, z=0.08)
            )
        )
    )

    if plot_type in ['histogram2d', 'scatter']:
        layout['xaxis'] = axis_template_2d(xlabel)
        layout['yaxis'] = axis_template_2d(ylabel)
        layout['plot_bgcolor'] = BACKGROUND
        layout['paper_bgcolor'] = BACKGROUND
        del layout['scene']
        del data[0]['z']

    if plot_type == 'histogram2d':
        # Scatter plot overlay on 2d Histogram
        data[0]['type'] = 'scatter'
        data.append( dict(
                x = x,
                y = y,
                type = 'histogram2d',
                colorscale = 'Greys',
                showscale = False
            ) )
        layout['plot_bgcolor'] = 'black'
        layout['paper_bgcolor'] = 'black'
        layout['xaxis'] = blackout_axis(layout['xaxis'])
        layout['yaxis'] = blackout_axis(layout['yaxis'])
        layout['font']['color'] = 'white'

    if len(markers) > 0:
        data = data + add_markers( data, markers, plot_type = plot_type )

    return dict( data=data, layout=layout )

def make_dash_table( selection ):
    ''' Return a dash defintion of an HTML table for a Pandas dataframe '''
    df_subset = df.loc[df['NAME'].isin(selection)]
    table = []
    for index, row in df_subset.iterrows():
        html_row = []
        for i in range(len(row)):
            if i == 0 or i == 6:
                html_row.append( html.Td([ row[i] ]) )
            elif i == 1:
                html_row.append( html.Td([ html.A( href=row[i], children='Datasheet' )]))
            elif i == 2:
                html_row.append( html.Td([ html.Img( src=row[i] )]))
        table.append( html.Tr( html_row ) )
    return table

FIGURE = scatter_plot_3d()
STARTING_DRUG = 'Levobupivacaine'
DRUG_DESCRIPTION = df.loc[df['NAME'] == STARTING_DRUG]['DESC'].iloc[0]
DRUG_IMG = df.loc[df['NAME'] == STARTING_DRUG]['IMG_URL'].iloc[0]

app.layout = html.Div([
    # Row 1: Header and Intro text

    html.Div([
        html.Img(src="https://s3-us-west-1.amazonaws.com/plotly-tutorials/logo/new-branding/dash-logo-by-plotly-stripe.png",
                style={
                    'height': '100px',
                    'float': 'right',
                    'position': 'relative',
                    'bottom': '40px',
                    'left': '50px'
                },
                ),
        html.H2('Dash',
                style={
                    'position': 'relative',
                    'top': '0px',
                    'left': '10px',
                    'font-family': 'Dosis',
                    'display': 'inline',
                    'font-size': '6.0rem',
                    'color': '#4D637F'
                }),
        html.H2('for',
                style={
                    'position': 'relative',
                    'top': '0px',
                    'left': '20px',
                    'font-family': 'Dosis',
                    'display': 'inline',
                    'font-size': '4.0rem',
                    'color': '#4D637F'
                }),
        html.H2('Drug Discovery',
                style={
                    'position': 'relative',
                    'top': '0px',
                    'left': '27px',
                    'font-family': 'Dosis',
                    'display': 'inline',
                    'font-size': '6.0rem',
                    'color': '#4D637F'
                }),
    ], className='row twelve columns', style={'position': 'relative', 'right': '15px'}),

    html.Div([
        html.Div([
            html.Div([
                html.P('HOVER over a drug in the graph to the right to see its structure to the left.'),
                html.P('SELECT a drug in the dropdown to add it to the drug candidates at the bottom.')
            ], style={'margin-left': '10px'}),
            dcc.Dropdown(id='chem_dropdown',
                        multi=True,
                        value=[ STARTING_DRUG ],
                        options=[{'label': i, 'value': i} for i in df['NAME'].tolist()]),
            ], className='twelve columns' )

    ], className='row' ),

    # Row 2: Hover Panel and Graph

    html.Div([
        html.Div([

            html.Img(id='chem_img', src=DRUG_IMG ),

            html.Br(),

            html.A(STARTING_DRUG,
                  id='chem_name',
                  href="https://www.drugbank.ca/drugs/DB01002",
                  target="_blank"),

            html.P(DRUG_DESCRIPTION,
                  id='chem_desc',
                  style=dict( maxHeight='400px', fontSize='12px' )),

        ], className='three columns', style=dict(height='300px') ),

        html.Div([

            dcc.RadioItems(
                id = 'charts_radio',
                options=[
                    dict( label='3D Scatter', value='scatter3d' ),
                    dict( label='2D Scatter', value='scatter' ),
                    dict( label='2D Histogram', value='histogram2d' ),
                ],
                labelStyle = dict(display='inline'),
                value='scatter3d'
            ),

            dcc.Graph(id='clickable-graph',
                      style=dict(width='700px'),
                      hoverData=dict( points=[dict(pointNumber=0)] ),
                      figure=FIGURE ),

        ], className='nine columns', style=dict(textAlign='center')),


    ], className='row' ),

    html.Div([
        html.Table( make_dash_table( [STARTING_DRUG] ), id='table-element' )
    ])

], className='container')


@app.callback(
    Output('clickable-graph', 'figure'),
    [Input('chem_dropdown', 'value'),
    Input('charts_radio', 'value')])
def highlight_molecule(chem_dropdown_values, plot_type):
    return scatter_plot_3d( markers = chem_dropdown_values, plot_type = plot_type )


@app.callback(
    Output('table-element', 'children'),
    [Input('chem_dropdown', 'value')])
def update_table(chem_dropdown_value):
    table = make_dash_table( chem_dropdown_value )
    return table


def dfRowFromHover( hoverData ):
    ''' Returns row for hover point as a Pandas Series '''
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                point_number = firstPoint['pointNumber']
                molecule_name = str(FIGURE['data'][0]['text'][point_number]).strip()
                return df.loc[df['NAME'] == molecule_name]
    return pd.Series()


@app.callback(
    Output('chem_name', 'children'),
    [Input('clickable-graph', 'hoverData')])
def return_molecule_name(hoverData):
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                point_number = firstPoint['pointNumber']
                molecule_name = str(FIGURE['data'][0]['text'][point_number]).strip()
                return molecule_name


@app.callback(
    dash.dependencies.Output('chem_name', 'href'),
    [dash.dependencies.Input('clickable-graph', 'hoverData')])
def return_href(hoverData):
    row = dfRowFromHover(hoverData)
    if row.empty:
        return
    datasheet_link = row['PAGE'].iloc[0]
    return datasheet_link


@app.callback(
    Output('chem_img', 'src'),
    [Input('clickable-graph', 'hoverData')])
def display_image(hoverData):
    row = dfRowFromHover(hoverData)
    if row.empty:
        return
    img_src = row['IMG_URL'].iloc[0]
    return img_src


@app.callback(
    Output('chem_desc', 'children'),
    [Input('clickable-graph', 'hoverData')])
def display_molecule(hoverData):
    row = dfRowFromHover(hoverData)
    if row.empty:
        return
    description = row['DESC'].iloc[0]
    return description


external_css = ["https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
                "//fonts.googleapis.com/css?family=Raleway:400,300,600",
                "//fonts.googleapis.com/css?family=Dosis:Medium",
                "https://cdn.rawgit.com/plotly/dash-app-stylesheets/0e463810ed36927caf20372b6411690692f94819/dash-drug-discovery-demo-stylesheet.css"]


for css in external_css:
    app.css.append_css({"external_url": css})


if __name__ == '__main__':
    app.run_server()
