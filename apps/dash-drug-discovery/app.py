import os
import dash
import pandas as pd
import pathlib
import dash_html_components as html
import dash_core_components as dcc

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

app = dash.Dash(__name__)
server = app.server

DATA_PATH = pathlib.Path(__file__).parent.joinpath("data").resolve()
ASSETS_PATH = pathlib.Path(__file__, "/assets")

# read from datasheet
df = pd.read_csv(DATA_PATH.joinpath("small_molecule_drugbank.csv")).drop(['Unnamed: 0'],axis=1)

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
                colorbar = {
                    "title": "Molecular<br>Weight",
                    "yanchor": "left",
                },
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

def make_dash_table(selection):
    """ Return a dash defintion of an HTML table for a Pandas dataframe. """

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

    html.Div([
        html.Img(src=str(ASSETS_PATH.joinpath("dash-logo-stripe.png").resolve()), className="app__banner__img")
    ], className="app__banner"),

    html.Div([
        html.Div([
            html.Div([
                html.H3("dash for drug discovery", className="uppercase"),
                html.Span("Hover ", className="uppercase bold"),
                html.Span("over a drug in the graph to see its structure."), 
                html.Br(),
                html.Span("Select ", className="uppercase bold"),
                html.Span("a drug in the dropdown to add it to the drug candidates at the bottom."), 
            ]),
        ], className="app__header"),

        html.Div([
            dcc.Dropdown(
                id='chem_dropdown',
                multi=True,
                value=[ STARTING_DRUG ],
                options=[{'label': i, 'value': i} for i in df['NAME'].tolist()]),
        ], className="app__dropdown"),

        html.Div([

            html.Div([
                dcc.RadioItems(
                    id = 'charts_radio',
                    options=[
                        {
                            "label": "3D Scatter",
                            "value": "scatter3d"
                        },
                        {
                            "label": "2D Scatter",
                            "value": "scatter"
                        },
                        {
                            "label": "2D Histogram",
                            "value": "histogram2d"
                        }
                    ],
                    labelClassName = "radio__labels",
                    value='scatter3d',
                    className="radio__group"
                ),
                dcc.Graph(
                    id='clickable-graph',
                    hoverData={"points": [{"pointNumber": 0}]},
                    figure=FIGURE 
                ),
            ], className="two-thirds column app__content__graph"),


            html.Div([
                html.Img(id='chem_img', src=DRUG_IMG ),
                html.Br(),
                html.A(STARTING_DRUG,
                    id='chem_name',
                    href="https://www.drugbank.ca/drugs/DB01002",
                    target="_blank"),
                html.P(DRUG_DESCRIPTION,
                    id='chem_desc'),
            ], className="one-third column app__content__info"),


        ], className="container card app__content")

    ], className="app__container"),

    html.Div([
        html.Table( make_dash_table([STARTING_DRUG]), id='table-element' )
    ])

], className="")


def df_row_from_hover(hoverData):
    """ Returns row for hover point as a Pandas Series. """

    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                point_number = firstPoint['pointNumber']
                molecule_name = str(FIGURE['data'][0]['text'][point_number]).strip()
                return df.loc[df['NAME'] == molecule_name]
    return pd.Series()


@app.callback(
    Output('clickable-graph', 'figure'),
    [   
        Input('chem_dropdown', 'value'),
        Input('charts_radio', 'value')
    ]
)
def highlight_molecule(chem_dropdown_values, plot_type):
    """
    Selected chemical dropdown values handler.

    :params chem_dropdown_values: selected dropdown values
    :params plot_type: selected plot graph
    """

    return scatter_plot_3d(markers = chem_dropdown_values, plot_type = plot_type)


@app.callback(
    Output('table-element', 'children'),
    [Input('chem_dropdown', 'value')])
def update_table(chem_dropdown_value):
    """
    Update the table rows.

    :params chem_dropdown_values: selected dropdown values
    """

    return make_dash_table(chem_dropdown_value)


@app.callback(
    [
        Output('chem_name', 'children'),
        Output('chem_name', 'href'),
        Output('chem_img', 'src'),
        Output('chem_desc', 'children')
    ],
    [Input('clickable-graph', 'hoverData')]
)
def chem_info_on_hover(hoverData):
    """ 
    Display chemical information on graph hover.
    Update the image, link, description.

    :params hoverData: data on graph hover
    """

    if hoverData is None:
        raise PreventUpdate

    try:
        row = df_row_from_hover(hoverData)
        if row.empty:
            raise Exception
        return row["NAME"].iloc[0], row['PAGE'].iloc[0], row['IMG_URL'].iloc[0], row['DESC'].iloc[0]

    except Exception as error:
        raise PreventUpdate


if __name__ == '__main__':
    app.run_server(debug=True)
