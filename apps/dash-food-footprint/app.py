import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_daq as daq
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px

path = 'https://raw.githubusercontent.com/InesRoque3/GroupV_project2/main/data/'

emissions = pd.read_csv(path + "emissions_with_origin.csv")
productions = pd.read_csv(path + "productions.csv")
water = pd.read_csv(path + "water_use.csv")
global_emissions = pd.read_csv(path + "Global_Emissions.csv")

top10 = emissions.sort_values("Total_emissions")[-10:]
top10_vegetal = emissions[emissions.Origin=='Vegetal'].sort_values("Total_emissions")[-10:]
top8_animal = emissions[emissions.Origin=='Animal'].sort_values("Total_emissions")



radio_ani_veg = dbc.RadioItems(
        id='ani_veg', 
        className='radio',
        options=[dict(label='Animal', value=0), dict(label='Vegetal', value=1), dict(label='Total', value=2)],
        value=2, 
        inline=True
    )

dict_ = {'Apples':'Apples', 'Bananas':'Bananas', 'Barley':'Barley', 'Beet Sugar':'Sugar beet', 'Berries & Grapes':'Berries & Grapes', 'Brassicas':'Brassicas', 
        'Cane Sugar':'Sugar cane', 'Cassava':'Cassava', 'Citrus Fruit':'Citrus', 'Coffee':'Coffee beans', 'Groundnuts':'Groundnuts','Maize':'Maize', 'Nuts':'Nuts', 
        'Oatmeal':'Oats', 'Olive Oil':'Olives', 'Onions & Leeks':'Onions & Leeks','Palm Oil':'Oil palm fruit', 'Peas':'Peas', 'Potatoes':'Potatoes', 'Rapeseed Oil':'Rapeseed',
        'Rice':'Rice', 'Root Vegetables':'Roots and tubers', 'Soymilk':'Soybeans', 'Sunflower Oil':'Sunflower seed', 'Tofu':'Soybeans','Tomatoes':'Tomatoes', 
        'Wheat & Rye':'Wheat & Rye', 'Dark Chocolate':'Cocoa, beans', 'Milk': 'Milk', 'Eggs': 'Eggs','Poultry Meat': 'Poultry Meat', 'Pig Meat': 'Pig Meat', 
        'Seafood (farmed)': 'Seafood (farmed)', 'Cheese': 'Cheese', 'Lamb & Mutton': 'Lamb & Mutton', 'Beef (beef herd)': 'Beef (beef herd)'}

options_veg = [dict(label=key, value=dict_[key]) for key in top10_vegetal['Food product'].tolist()[::-1] if key in dict_.keys()]
options_an = [dict(label=val, value=val) for val in top8_animal["Food product"].tolist()[::-1]]
options_total = [dict(label=key, value=dict_[key]) for key in top10['Food product'].tolist()[::-1] if key in dict_.keys()]

bar_colors = ['#ebb36a','#6dbf9c']
bar_options = [top8_animal, top10_vegetal, top10]

drop_map = dcc.Dropdown(
        id = 'drop_map',
        clearable=False,
        searchable=False, 
        style= {'margin': '4px', 'box-shadow': '0px 0px #ebb36a', 'border-color': '#ebb36a'}        
    )

drop_continent = dcc.Dropdown(
        id = 'drop_continent',
        clearable=False, 
        searchable=False, 
        options=[{'label': 'World', 'value': 'world'},
                {'label': 'Europe', 'value': 'europe'},
                {'label': 'Asia', 'value': 'asia'},
                {'label': 'Africa', 'value': 'africa'},
                {'label': 'North america', 'value': 'north america'},
                {'label': 'South america', 'value': 'south america'}],
        value='world', 
        style= {'margin': '4px', 'box-shadow': '0px 0px #ebb36a', 'border-color': '#ebb36a'}
    )

slider_map = daq.Slider(
        id = 'slider_map',
        handleLabel={"showCurrentValue": True,"label": "Year"},
        marks = {str(i):str(i) for i in [1990,1995,2000,2005,2010,2015]},
        min = 1990,
        size=450, 
        color='#4B9072'
    )

fig_water = px.sunburst(water, path=['Origin', 'Category', 'Product'], values='Water Used', color='Category', 
                        color_discrete_sequence = px.colors.sequential.haline_r).update_traces(hovertemplate = '%{label}<br>' + 'Water Used: %{value} L')

fig_water = fig_water.update_layout({'margin' : dict(t=0, l=0, r=0, b=10),
                        'paper_bgcolor': '#F9F9F8',
                        'font_color':'#363535'
                    })

fig_gemissions = px.sunburst(global_emissions, path = ['Emissions', 'Group','Subgroup'], values = 'Percentage of food emissions', 
                    color = 'Group', color_discrete_sequence = px.colors.sequential.Peach_r).update_traces(hovertemplate = '%{label}<br>' + 'Global Emissions: %{value}%', textinfo = "label + percent entry") 

fig_gemissions = fig_gemissions.update_layout({'margin' : dict(t=0, l=0, r=0, b=10),
                        'paper_bgcolor': '#F9F9F8',
                        'font_color':'#363535'})


#------------------------------------------------------ APP ------------------------------------------------------ 

app = dash.Dash(__name__)

server = app.server

app.layout = html.Div([

    html.Div([
        html.H1(children='FOOD FOOTPRINT'),
        html.Label('We are interested in investigating the food products that have the biggest impact on environment. Here you can understand which are the products whose productions emit more greenhouse gases and associate this with each supply chain step, their worldwide productions, and the water use.', 
                    style={'color':'rgb(33 36 35)'}), 
        html.Img(src=app.get_asset_url('supply_chain.png'), style={'position': 'relative', 'width': '180%', 'left': '-83px', 'top': '-20px'}),
    ], className='side_bar'),

    html.Div([
        html.Div([
            html.Div([
                html.Label("Choose the Product's Origin:"), 
                html.Br(),
                html.Br(),
                radio_ani_veg
            ], className='box', style={'margin': '10px', 'padding-top':'15px', 'padding-bottom':'15px'}),

            html.Div([
                html.Div([

                    html.Div([   
                        html.Label(id='title_bar'),           
                        dcc.Graph(id='bar_fig'), 
                        html.Div([              
                            html.P(id='comment')
                        ], className='box_comment'),
                    ], className='box', style={'padding-bottom':'15px'}),

                    html.Div([
                        html.Img(src=app.get_asset_url('Food.png'), style={'width': '100%', 'position':'relative', 'opacity':'80%'}),
                    ]),

                ], style={'width': '40%'}),


                html.Div([

                    html.Div([
                    html.Label(id='choose_product', style= {'margin': '10px'}),
                    drop_map,
                    ], className='box'),

                    html.Div([
                        html.Div([
                            html.Label('Emissions measured as kg of CO2 per kg of product', style={'font-size': 'medium'}),
                            html.Br(),
                            html.Br(),
                            html.Div([
                                html.Div([
                                    html.H4('Land use', style={'font-weight':'normal'}),
                                    html.H3(id='land_use')
                                ],className='box_emissions'),

                                html.Div([
                                    html.H4('Animal Feed', style={'font-weight':'normal'}),
                                    html.H3(id='animal_feed')
                                ],className='box_emissions'),
                            
                                html.Div([
                                    html.H4('Farm', style={'font-weight':'normal'}),
                                    html.H3(id='farm')
                                ],className='box_emissions'),

                                html.Div([
                                    html.H4('Processing', style={'font-weight':'normal'}),
                                    html.H3(id='processing')
                                ],className='box_emissions'),
                            
                                html.Div([
                                    html.H4('Transport', style={'font-weight':'normal'}),
                                    html.H3(id='transport')
                                ],className='box_emissions'),

                                html.Div([
                                    html.H4('Packaging', style={'font-weight':'normal'}),
                                    html.H3(id='packging')
                                ],className='box_emissions'),
                            
                                html.Div([
                                    html.H4('Retail', style={'font-weight':'normal'}),
                                    html.H3(id='retail')
                                ],className='box_emissions'),
                            ], style={'display': 'flex'}),

                        ], className='box', style={'heigth':'10%'}),

                        html.Div([ 
                            html.Div([
                                
                                html.Div([
                                    html.Br(),
                                    html.Label(id='title_map', style={'font-size':'medium'}), 
                                    html.Br(),
                                    html.Label('These quantities refer to the raw material used to produce the product selected above', style={'font-size':'9px'}),
                                ], style={'width': '70%'}),
                                html.Div([

                                ], style={'width': '5%'}),
                                html.Div([
                                    drop_continent, 
                                    html.Br(),
                                    html.Br(), 
                                ], style={'width': '25%'}),
                            ], className='row'),
                            
                            dcc.Graph(id='map', style={'position':'relative', 'top':'-50px'}), 

                            html.Div([
                                slider_map
                            ], style={'margin-left': '15%', 'position':'relative', 'top':'-38px'}),
                            
                        ], className='box', style={'padding-bottom': '0px'}), 
                    ]),
                ], style={'width': '60%'}),           
            ], className='row'),

            html.Div([
                html.Div([
                    html.Label("3. Global greenhouse gas emissions from food production, in percentage", style={'font-size': 'medium'}),
                    html.Br(),
                    html.Label('Click on it to know more!', style={'font-size':'9px'}),
                    html.Br(), 
                    html.Br(), 
                    dcc.Graph(figure=fig_gemissions)
                ], className='box', style={'width': '40%'}), 
                html.Div([
                    html.Label("4. Freshwater withdrawals per kg of product, in Liters", style={'font-size': 'medium'}),
                    html.Br(),
                    html.Label('Click on it to know more!', style={'font-size':'9px'}),
                    html.Br(), 
                    html.Br(), 
                    dcc.Graph(figure=fig_water)
                ], className='box', style={'width': '63%'}), 
            ], className='row'),

            html.Div([
                html.Div([
                    html.P(['GroupV', html.Br(),'Ana Carrelha (20200631), Inês Melo (20200624), Inês Roque (20200644), Ricardo Nunes(20200611)'], style={'font-size':'12px'}),
                ], style={'width':'60%'}), 
                html.Div([
                    html.P(['Sources ', html.Br(), html.A('Our World in Data', href='https://ourworldindata.org/', target='_blank'), ', ', html.A('Food and Agriculture Organization of the United Nations', href='http://www.fao.org/faostat/en/#data', target='_blank')], style={'font-size':'12px'})
                ], style={'width':'37%'}),
            ], className = 'footer', style={'display':'flex'}),
        ], className='main'),
    ]),
])


#------------------------------------------------------ Callbacks ------------------------------------------------------

@app.callback(
    [
        Output('title_bar', 'children'),
        Output('bar_fig', 'figure'),
        Output('comment', 'children'),
        Output('drop_map', 'options'),
        Output('drop_map', 'value'),
        Output('choose_product', 'children')
    ],
    [
        Input('ani_veg', 'value')
    ], 
)
def bar_chart(top10_select):

    ################## Top10 Plot ##################
    title = '1. Greenhouse emissions (kg CO2 per kg of product)'
    df = bar_options[top10_select]

    if top10_select==2:
        bar_fig = dict(type='bar',
            x=df.Total_emissions,
            y=df["Food product"],
            orientation='h',
            marker_color=['#ebb36a' if x=='Animal' else '#6dbf9c' for x in df.Origin])
    else:
        bar_fig = dict(type='bar',
            x=df.Total_emissions,
            y=df["Food product"],
            orientation='h',
            marker_color=bar_colors[top10_select])

    ################## Dropdown Bar ##################
    if top10_select==0:
        options_return = options_an
        product_chosen = "2. Choose an animal product:" 
        comment = ["Look at the beef production emissions! Each kilogram of beef produces almost 60 kg of CO2.", html.Br(), html.Br()]
    elif top10_select==1:
        options_return = options_veg
        product_chosen = "2. Choose a vegetal product:" 
        comment = ["Did you know that dark chocolate and coffee are the vegetal-based products that emit more gases?", html.Br(), html.Br()]
    else:
        options_return = options_total
        product_chosen = "2. Choose an animal or vegetal product:" 
        comment = "Check the difference between animal and vegetal-based products! Beef (top1 animal-based emitter) produces around 3 times more emissions than dark chocolate (top1 plant-based emitter)."

    return title, \
            go.Figure(data=bar_fig, layout=dict(height = 300, font_color = '#363535', paper_bgcolor='rgba(0,0,0,0)',plot_bgcolor='rgba(0,0,0,0)', margin=dict(l=20, r=20, t=30, b=20), margin_pad=10)), \
            comment, \
            options_return, \
            options_return[0]['value'], \
            product_chosen



@app.callback(
    [ 
        Output('slider_map', 'max'),
        Output('slider_map', 'value'),
    ],
    [
        Input('drop_map', 'value')
    ]
)

def update_slider(product):
    year = productions[productions['Item']==product]['Year'].max()
    return year, year



@app.callback(
    [
        Output('land_use', 'children'),
        Output('animal_feed', 'children'),
        Output('farm', 'children'),
        Output('processing', 'children'),
        Output('transport', 'children'),
        Output('packging', 'children'),
        Output('retail', 'children'),
        Output('title_map', 'children'),
        Output('map', 'figure')
    ],
    [
        Input('drop_map', 'value'),
        Input('slider_map', 'value'), 
        Input('drop_continent', 'value')
    ],
    [State("drop_map","options")]
)

def update_map(drop_map_value, year, continent, opt):

    ################## Emissions datset ##################
    
    the_label = [x['label'] for x in opt if x['value'] == drop_map_value]

    data_emissions = emissions[emissions["Food product"]==the_label[0]]
    land_use_str = str(np.round(data_emissions["Land use change"].values[0],2))
    animal_feed_str = str(np.round(data_emissions["Animal Feed"].values[0],2))
    farm_str = str(np.round(data_emissions["Farm"].values[0],2))
    processing_str = str(np.round(data_emissions["Processing"].values[0],2))
    transport_str = str(np.round(data_emissions["Transport"].values[0],2))
    packging_str = str(np.round(data_emissions["Packging"].values[0],2))
    retail_str = str(np.round(data_emissions["Retail"].values[0],2))

    ################## Choroplet Plot ##################
    
    prod1 = productions[(productions['Item']== drop_map_value) & (productions['Year']== year)]

    title ='Production quantities of {}, by country'.format(prod1['Item'].unique()[0])  #font_color = '#363535',
    data_slider = []
    data_each_yr = dict(type='choropleth',
                        locations = prod1['Area'],
                        locationmode='country names',
                        autocolorscale = False,
                        z=np.log(prod1['Value'].astype(float)),
                        zmin=0, 
                        zmax = np.log(productions[productions['Item']== drop_map_value]['Value'].max()),
                        colorscale = ["#ffe2bd", "#006837"],   
                        marker_line_color= 'rgba(0,0,0,0)',
                        colorbar= {'title':'Tonnes (log)'},#Tonnes in logscale
                        colorbar_lenmode='fraction',
                        colorbar_len=0.8,
                        colorbar_x=1,
                        colorbar_xanchor='left',
                        colorbar_y=0.5,
                        name='')
    data_slider.append(data_each_yr)
 
    layout = dict(geo=dict(scope=continent,
                            projection={'type': 'natural earth'},
                            bgcolor= 'rgba(0,0,0,0)'),
                    margin=dict(l=0,
                                r=0,
                                b=0,
                                t=30,
                                pad=0),
                    paper_bgcolor='rgba(0,0,0,0)',
                    plot_bgcolor='rgba(0,0,0,0)')
 
    fig_choropleth = go.Figure(data=data_slider, layout=layout)
    fig_choropleth.update_geos(showcoastlines=False, showsubunits=False,showframe=False)

    return land_use_str, \
        animal_feed_str, \
        farm_str, \
        processing_str, \
        transport_str, \
        packging_str, \
        retail_str, \
        title, \
        fig_choropleth


if __name__ == '__main__':
    app.run_server(debug=True)
