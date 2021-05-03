import json
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import urllib.request as urllib
#from urllib.request import urlopen
import pandas as pd
import plotly.graph_objs as go
import statsmodels.api as sm
import plotly.figure_factory as ff
import numpy as np

mapbox_access_token = "pk.eyJ1Ijoic3RlZmZlbmhpbGwiLCJhIjoiY2ttc3p6ODlrMG1ybzJwcG10d3hoaDZndCJ9.YE2gGNJiw6deBuFgHRHPjg"
path = 'https://raw.githubusercontent.com/FranzMichaelFrank/health_eu/main/'
df = pd.read_csv(path + 'food_supply.csv', dtype={"id": str})
df_scatter = pd.read_csv(path + 'scatter_data.csv', dtype={"id": str})
clusters = pd.read_csv(path + "final_clusters.csv")
box_cluster = pd.read_csv(path + "box_cluster.csv")

url = path + "european-union-countries.geojson"
response = urllib.urlopen(url)
european_union = json.loads(response.read())

def cz(s):
    if s == "Czech Republic":
        return "Czechia"
    else:
        return s
clusters = clusters.rename(columns={"cluster_food":"Food clusters", "cluster_health":"Health clusters", "cluster_total":"Food & Health clusters"})
clusters["Food clusters"] = clusters["Food clusters"] + 1
clusters["Food clusters"] = clusters["Food clusters"].apply(str)
clusters["Health clusters"] = clusters["Health clusters"] + 1
clusters["Health clusters"] = clusters["Health clusters"].apply(str)
clusters["Food & Health clusters"] = clusters["Food & Health clusters"] + 1
clusters["Food & Health clusters"] = clusters["Food & Health clusters"].apply(str)

clusters["Country"] = clusters["Country"].apply(cz)
df["Country"] = df["Country"].apply(cz)

layout = dict(
    autosize=True,
    #automargin=True,
    margin=dict(l=30, r=30, b=20, t=40),
    hovermode="closest",
    plot_bgcolor="#F9F9F9",
    paper_bgcolor="#F9F9F9",
    legend=dict(font=dict(size=10), orientation="h"),
    title="Satellite Overview",
    mapbox=dict(
        accesstoken=mapbox_access_token,
        style="light",
        center=dict(lon=-78.05, lat=42.54),
        zoom=7,
    ),
)

# design for mapbox
bgcolor = "#f3f3f1"  # mapbox light map land color
row_heights = [150, 500, 300]
template = {"layout": {"paper_bgcolor": bgcolor, "plot_bgcolor": bgcolor}}

def blank_fig(height):
    """
    Build blank figure with the requested height
    """
    return {
        "data": [],
        "layout": {
            "height": height,
            "template": template,
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
        },
    }

columns = df.iloc[:,2:22].columns
health_cols = ["Obesity", "Diabetes Prevalence", "Cardiovascular Death Rate", "Life Expectancy", "Health Expenditure" ]

# Create controls
behaviour_options = [
    dict(label=country, value=country)
    for country in columns]



food_options_ = ["Alcoholic Beverages", "Animal fats", "Cereals - Excluding Beer", "Eggs", 
                "Fish, Seafood", "Fruits - Excluding Wine", "Meat", "Milk - Excluding Butter",
                "Offals", "Oilcrops", "Pulses", "Spices", "Starchy Roots","Stimulants", 
                "Sugar & Sweeteners", "Treenuts", "Vegetable Oils",
                "Vegetables" ]

food_options = [
    dict(label=country, value=country)
    for country in food_options_]




dropdown_behaviour = dcc.Dropdown(
        id='candidate_radio',
        options=behaviour_options,
        value=columns[0]
    )

ots = ["Alcoholic Beverages", "Animal fats", "Cereals - Excluding Beer", "Eggs", 
                "Fish, Seafood", "Fruits - Excluding Wine", "Meat", "Milk - Excluding Butter",
                "Offals", "Oilcrops", "Pulses", "Spices", "Starchy Roots","Stimulants", 
                "Sugar & Sweeteners", "Treenuts", "Vegetable Oils",
                "Vegetables","Obesity", "Diabetes Prevalence", "Cardiovascular Death Rate", "Life Expectancy", "Health Expenditure" ]
ots_ = [
    dict(label=country, value=country)
    for country in ots]

ots_behaviour = dcc.Dropdown(
        id='box_dd',
        options=ots_,
        value="Alcoholic Beverages"
    )

radio_food_behaviour = dcc.RadioItems(
                            id='nutrition_types',
                            options=food_options,
                            value="Alcoholic Beverages",
                            labelStyle={'display': 'block', "text-align": "justify"}
                            
                        )

                        
corr_options = [
    dict(label=country, value=country)
    for country in ["Obesity", "Diabetes Prevalence", "Cardiovascular Death Rate", "Life Expectancy", "Health Expenditure"]]
cor_behav = dcc.Dropdown(
                            id='cor_behave',
                            options=corr_options,
                            value="Obesity"#,
                            #labelStyle={'display': 'block', "text-align": "justify"}
                        )

clust_opt = ["Food clusters", "Health clusters", "Food & Health clusters"]
clust_options = [
    dict(label=country, value=country)
    for country in clust_opt]
radio_clust_behaviour = dcc.RadioItems(
                            id='clust_types',
                            options=clust_options,
                            value="Food clusters",
                            labelStyle={"display":"inline-flex", "text-align":"center"},
                            style={"padding-left": "10%"}
                        )

app = dash.Dash(__name__)

###
df_new = pd.read_csv('https://plotly.github.io/datasets/country_indicators.csv')
available_indicators = df_new['Indicator Name'].unique()
###


## FF ##
foods = ["Alcoholic Beverages", "Animal fats", "Cereals - Excluding Beer", "Eggs", 
                "Fish, Seafood", "Fruits - Excluding Wine", "Meat", "Milk - Excluding Butter",
                "Offals", "Oilcrops", "Pulses", "Spices", "Starchy Roots","Stimulants", 
                "Sugar & Sweeteners", "Treenuts", "Vegetable Oils",
                "Vegetables" ]

foods_r = ['Vegetables', 'Vegetable Oils', 'Treenuts', 'Sugar & Sweeteners', 'Stimulants', 'Starchy Roots', 'Spices', 'Pulses', 'Oilcrops', 'Offals', 'Milk - Excluding Butter', 'Meat', 'Fruits - Excluding Wine', 'Fish, Seafood', 'Eggs', 'Cereals - Excluding Beer', 'Animal fats', 'Alcoholic Beverages']
foods_health = ['Vegetables', 'Vegetable Oils', 'Treenuts', 'Sugar & Sweeteners', 'Stimulants', 'Starchy Roots', 'Spices', 'Pulses', 'Oilcrops', 'Offals', 'Milk - Excluding Butter', 'Meat', 'Fruits - Excluding Wine', 'Fish, Seafood', 'Eggs', 'Cereals - Excluding Beer', 'Animal fats', 'Alcoholic Beverages', "Obesity", "Diabetes Prevalence", "Cardiovascular Death Rate", "Life Expectancy", "Health Expenditure" ]

df_corr_r = df_scatter[foods_health]
df_corr_round = df_corr_r.corr()[["Obesity"]].T[foods_r].T.round(2)
#, "Diabetes Prevalence", "Cardiovascular Death Rate", "Life Expectancy", "Health Expenditure" 
fig_cor = ff.create_annotated_heatmap(
            z=df_corr_round.to_numpy(),
            x=df_corr_round.columns.tolist(),
            y=df_corr_round.index.tolist(),
            zmax=1, zmin=-1,
            showscale=True,
            hoverongaps=True
            )
fig_cor.update_layout(legend=dict(
    orientation="h",
    yanchor="bottom",
    y=1.02,
    xanchor="right",
    x=1
))
#fig_cor.update_layout(yaxis_tickangle=-45)
fig_cor.update_layout(xaxis_tickangle=0)
fig_cor.update_layout(title_text='',height=600)



health = df_scatter[["Country","Obesity", "Diabetes Prevalence", "Cardiovascular Death Rate", "Life Expectancy", "Health Expenditure" ]]



fig_bar = go.Figure()
fig_bar.add_trace(go.Bar(
    x=health["Country"],
    y=health["Obesity"],
    name="Obesity",
    marker_color='#0d0887'
))
fig_bar.add_trace(go.Bar(
    x=health["Country"],
    y=health["Diabetes Prevalence"],
    name="Diabetes Prevalence",
    marker_color='#7201a8'
))

fig_bar.add_trace(go.Bar(
    x=health["Country"],
    y=health["Cardiovascular Death Rate"],
    name="Cardiovascular Death Rate",
    marker_color='#bd3786'
))

fig_bar.add_trace(go.Bar(
    x=health["Country"],
    y=health["Life Expectancy"],
    name="Life Expectancy",
    marker_color='#ed7953'
))

fig_bar.add_trace(go.Bar(
    x=health["Country"],
    y=health["Health Expenditure"],
    name="Health Expenditure",
    marker_color='#fdca26'
))
# Here we modify the tickangle of the xaxis, resulting in rotated labels.
fig_bar.update_layout(barmode='group', xaxis_tickangle=-45)
fig_bar.update_layout(plot_bgcolor='white')
fig_bar.update_yaxes(showline=True, linewidth=2, linecolor='black', gridcolor='grey')
fig_bar.update_xaxes(showline=True, linewidth=2, linecolor='black')



## FF ##

# Create app layout
app.layout = html.Div(
    [
        dcc.Store(id="aggregate_data"),
        # empty Div to trigger javascript file for graph resizing
        html.Div(id="output-clientside"),
        html.Div(
            [
                html.Div(
                    [
                        html.Img(
                            src=app.get_asset_url("Nova_IMS.png"),
                            id="plotly-image",
                            style={
                                "height": "60px",
                                "width": "auto",
                                "margin-bottom": "25px",
                            },
                        )
                    ],
                    className="one-third column",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.H4(
                                    "Food consumption characteristics of the European Union",
                                    style={"font-weight": "bold"},
                                ),
                                html.H5(
                                    "Analysis of the relationship between nutritional patterns and \n the health status within the countries", style={"margin-top": "0px"}
                                ),
                            ]
                        )
                    ],
                    className="three column",
                    id="title",
                ),
                html.Div(
                    # create empty div for align cente
                    className="one-third column",
                ),
            ],
            id="header",
            className="row flex-display",
            style={"margin-bottom": "25px"},
        ),
        html.Div(
            [
                html.Div(
                    [

                        html.H6("Consumption by food type", style={"margin-top": "0","font-weight": "bold","text-align": "center"}),

                        html.P(
                            "The cultures and customs of the 27 EU countries differ widely. The same applies to their eating and drinking habits." 
                            " The map on the right explores the food supply in kilograms per capita per year.",
                            className="control_label",style={"text-align": "justify"}
                        ),
                        html.P(),
                        html.P("Select a food category", className="control_label",style={"text-align": "center","font-weight":"bold"}),
                        radio_food_behaviour,

                    ],
                    className="pretty_container four columns",
                    id="cross-filter-options",
                    style={"text-align": "justify"},

                ),

                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.P(id="well_text"), html.P("Maximum",style={"text-align": "center","font-weight":"bold"}),
                                        html.P(id="max_name",style={"text-align": "center"}),
                                        html.P(id="max_value",style={"text-align": "center"}),
                                    ],
                                        className="mini_container",
                                        id="wells",
                                ),
                                html.Div(
                                    [html.P(id="gasText"), html.P("Minimum",style={"text-align": "center","font-weight":"bold"}),
                                    html.P(id="min_name",style={"text-align": "center"}),
                                     html.P(id="min_value",style={"text-align": "center"}),
                                    ],
                                    className="mini_container",
                                    id="gas"
                                ),
                                html.Div(
                                    [html.P(id="oilText"), html.P("Mean",style={"text-align": "center","font-weight":"bold"}),
                                     html.P(id="mean",style={"text-align": "center"}),
                                     html.P("Standard deviation",style={"text-align": "center","font-weight":"bold"}), 
                                     html.P(id="st_dev",style={"text-align": "center"})],
                                    #,
                                    className="mini_container",
                                    id="oil",
                                ),
                            ],
                            id="info-container",
                            className="row container-display",
                        ),
                        html.Div(
                            [dcc.Graph(id="choropleth")],
                            #id="countGraphContainer",
                            className="pretty_container",
                        ),
                    ],
                    id="right-column",
                    className="eight columns",
                ),
            ],
            className="row flex-display",
        ),

        html.Div(
            [
                html.H6("General health information about the countries", style={"margin-top": "0","font-weight": "bold","text-align": "center"}),
                html.P("Similarly to nutrition, the health status also varies from country to country. The bar chart below shows the differences between the countries in terms of the following variables: prevalence of obesity in the adult population in % (Obesity), prevalence of diabetes in the adult population in % (Diabetes Prevalence), cardiovascular death rate per 100,000 population (Cardiovascular Death Rate), average life expectancy in years (Life Expectancy) and the expenditure of the government on the country's health system in % of the respective GDP (Health Expenditure).", className="control_label",style={"text-align": "justify"}),
                
                html.Div([dcc.Graph(id="bar_chart", figure=fig_bar)],className="pretty_container twelve columns"),
            ],
            className="row pretty_container",
        ),

        html.Div(
            [   
                html.Div(
                    #cor_behav,
                    
                    [
                    html.H6("Exploring correlations", style={"margin-top": "0","font-weight": "bold","text-align": "center"}),
                    html.P("In the heatmap below, the correlations between the 5 health variables and the 18 food variables can be explored.", className="control_label",style={"text-align": "justify"}),
                    #html.P("""<br>"""),
                    html.P("Select a health variable", style={"font-weight": "bold", "text-align": "center"}),
                    cor_behav,

                    html.Div([dcc.Graph(id="cor_ma")],className="pretty_container twelve columns"),
                    ],#,cor_behav,
                    className="pretty_container four columns",
                    

                ),

                
                html.Div(
                    [   html.H6("Analysing the correlations between food consumption and health", style={"margin-top": "0","font-weight": "bold","text-align": "center"}),
                        html.P("Below, the correlations can be analysed in more detail. It is important to note that correlation in this case does not necessarily mean causation. For example, the wealth level of a country might often influence the variables, which is why it is represented through the size of the dots as GDP per capita. Furthermore, outliers should also be watched out for.", className="control_label",style={"text-align": "justify"}),

                        html.Div([
                        html.P("Select a food category", className="control_label", style={"font-weight": "bold", "text-align": "center"}),
                        dcc.Dropdown(
                            id='xaxis-column',
                            options=[{'label': i, 'value': i} for i in food_options_],
                            value="Alcoholic Beverages"#,className="pretty_container four columns",
                        ),
                        dcc.RadioItems(
                            id='xaxis-type',
                            options=[{'label': i, 'value': i} for i in ['Box', 'Violin']],
                            value='Box',
                            labelStyle={'display': 'inline-block'},#,className="pretty_container four columns",
                            style={"padding-left": "34%"}
                        ),
                        ],className="pretty_container sixish columns",),

                        html.Div([
                        html.P("Select a health variable", className="control_label", style={"font-weight": "bold", "text-align": "center"}),
                        dcc.Dropdown(
                            id='yaxis-column',
                            options=[{'label': i, 'value': i} for i in health_cols],
                            value=health_cols[0]
                        ),
                        dcc.RadioItems(
                            id='yaxis-type',
                            options=[{'label': i, 'value': i} for i in ['Box', 'Violin']],
                            value='Box',
                            labelStyle={'display': 'inline-block'},
                            style={"padding-left": "34%"}
                        ),
                        ],className="pretty_container sixish columns",),

                        html.Div([
                        dcc.Graph(id="indicator-graphic"),
                        ],className="pretty_container almost columns",),
                    ],
                    
                    className="pretty_container eight columns",
                ),

                
            ],
            className="row flex-display",
        ),

        

        html.Div(
            [
                
                    html.H6("K-means clustering", style={"margin-top": "0","font-weight": "bold","text-align": "center"}),
                    html.P("Finally, k-means clustering is carried out. The criterion can be selected on the left side, whereby either the 18 food variables, the 5 health variables or all of them in combination may be chosen for the clustering. On the right side, the resulting clusters can then be compared with respect to a chosen variable.",className="control_label",style={"text-align": "justify"}),
                    html.Div([ 
                    html.P("Select a clustering criterion", className="control_label",style={"text-align": "center","font-weight":"bold"}),
                    radio_clust_behaviour,
                    
                    dcc.Graph(id="cluster_map")
                    ],className="pretty_container sixish columns",), 

                    html.Div([ 
                    
                    html.P("Select a variable for cluster comparison", className="control_label",style={"text-align": "center","font-weight":"bold"}),
                    ots_behaviour,
                    dcc.Graph(id="boxes"),
                    ],className="pretty_container sixish columns",), 

                    
                
            ],
            className="row pretty_container",
        ),

        html.Div(
            [
                html.H6("Authors", style={"margin-top": "0","font-weight": "bold","text-align": "center"}),
                
                html.P("Maximilian Maukner (m20200645@novaims.unl.pt)  -  Ehsan Meisami Fard (m20201050@novaims.unl.pt)  -  Franz Michael Frank (m20200618@novaims.unl.pt)  -  Steffen Hillmann (m20200589@novaims.unl.pt)", style={"text-align": "center", "font-size":"10pt"}),
                
            ],
            className="row pretty_container",
        ),
        html.Div(
            [
                html.H6("Sources", style={"margin-top": "0","font-weight": "bold","text-align": "center"}),
                dcc.Markdown(
                    """\
                         - Eurostat: https://ec.europa.eu/eurostat/databrowser/view/HLTH_SHA11_HF__custom_227597/bookmark/table?lang=en&bookmarkId=1530a1e6-767e-4661-9e15-0ed2f7fae0d5
                         - Food and Agriculture Organization of the United Nations: http://www.fao.org/faostat/en/#data/FBS
                         - Opendatasoft: https://data.opendatasoft.com/explore/dataset/european-union-countries@public/export/
                         - Our World in Data: https://covid.ourworldindata.org/data/owid-covid-data.csv?v=2021-03-11
                        """
                ,style={"font-size":"10pt"}),
                
            ],
            className="row pretty_container",
        ),



        
    ],
    id="mainContainer",
    style={"display": "flex", "flex-direction": "column"},
)

colors = ['#0d0887', '#46039f', '#7201a8', '#9c179e', '#bd3786', '#d8576b', '#ed7953', '#fb9f3a', '#fdca26', '#f0f921']
colors2 = ['#fdca26', '#ed7953', '#bd3786', '#7201a8','#0d0887']

####
@app.callback(
    Output("choropleth", "figure"),
    [Input("nutrition_types", "value")])
def display_choropleth(candi):
    fig = px.choropleth_mapbox(
        df, geojson=european_union, color=candi,
        locations="iso_a3", featureidkey="properties.gu_a3",hover_name = "Country" ,opacity=0.7,  #hover_data = [],
        center={"lat": 56.5, "lon": 11},
        zoom=2.5)
    fig.update_layout(
        margin={"r":0,"t":0,"l":0,"b":0},
        mapbox_accesstoken=mapbox_access_token)

    return fig

@app.callback(
    Output("boxes", "figure"),
    [Input("clust_types", "value"), Input("box_dd", "value")])
def display_boxes(clust, vari):
    clusts = ["Food", "Health", "Food & Health"]

    
    box_cluster = pd.read_csv(path + "box_cluster.csv").rename(columns={"Food":"Food clusters", "Health":"Health clusters", "Food & Health":"Food & Health clusters"})
    box_cluster["Food clusters"] = box_cluster["Food clusters"] + 1
    box_cluster["Health clusters"] = box_cluster["Health clusters"] + 1
    box_cluster["Food & Health clusters"] = box_cluster["Food & Health clusters"] + 1
    box_cluster = box_cluster.sort_values(clust)
    fig = px.box(box_cluster, x=clust, y=vari, color=clust, color_discrete_sequence=['#fdca26', '#ed7953', '#bd3786', '#7201a8', '#0d0887'])
    fig.update_layout(plot_bgcolor='white')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_layout(legend=dict(yanchor='bottom', y=0))
    
    return fig

@app.callback(
    Output("cluster_map", "figure"),
    [Input("clust_types", "value")])
def display_cluster_map(type_clust):
    fig = px.choropleth_mapbox(
    clusters.sort_values(type_clust),
    locations="id",
    geojson=european_union,featureidkey="properties.gu_a3",
    color=type_clust,
    hover_name="Country",
    hover_data=[type_clust],
    #hover_data=["Life Expectancy"],
    #title="clusters by food",
    mapbox_style="carto-positron", 
    center={"lat": 56.5, "lon": 11},
    zoom=2.5,
    opacity=0.7,
    color_discrete_sequence = colors2
)   
    
    fig.update_layout(legend=dict(yanchor='bottom', y=0),
        margin={"r":0,"t":0,"l":0,"b":0},
        mapbox_accesstoken=mapbox_access_token)

    return fig

## FF ##
@app.callback(
    Output("cor_ma", "figure"),
    [Input("cor_behave", "value")])
def display_cor_ma(var):
    foods = ["Alcoholic Beverages", "Animal fats", "Cereals - Excluding Beer", "Eggs", 
                "Fish, Seafood", "Fruits - Excluding Wine", "Meat", "Milk - Excluding Butter",
                "Offals", "Oilcrops", "Pulses", "Spices", "Starchy Roots","Stimulants", 
                "Sugar & Sweeteners", "Treenuts", "Vegetable Oils",
                "Vegetables" ]

    foods_r = ['Vegetables', 'Vegetable Oils', 'Treenuts', 'Sugar & Sweeteners', 'Stimulants', 'Starchy Roots', 'Spices', 'Pulses', 'Oilcrops', 'Offals', 'Milk - Excluding Butter', 'Meat', 'Fruits - Excluding Wine', 'Fish, Seafood', 'Eggs', 'Cereals - Excluding Beer', 'Animal fats', 'Alcoholic Beverages']
    foods_health = ['Vegetables', 'Vegetable Oils', 'Treenuts', 'Sugar & Sweeteners', 'Stimulants', 'Starchy Roots', 'Spices', 'Pulses', 'Oilcrops', 'Offals', 'Milk - Excluding Butter', 'Meat', 'Fruits - Excluding Wine', 'Fish, Seafood', 'Eggs', 'Cereals - Excluding Beer', 'Animal fats', 'Alcoholic Beverages', "Obesity", "Diabetes Prevalence", "Cardiovascular Death Rate", "Life Expectancy", "Health Expenditure" ]

    df_corr_r = df_scatter[foods_health]
    df_corr_round = df_corr_r.corr()[[var]].T[foods_r].T.round(2)
    #, "Diabetes Prevalence", "Cardiovascular Death Rate", "Life Expectancy", "Health Expenditure" 
    fig_cor = ff.create_annotated_heatmap(
            z=df_corr_round.to_numpy(),
            x=df_corr_round.columns.tolist(),
            y=df_corr_round.index.tolist(),
            zmax=1, zmin=-1,
            showscale=True,
            hoverongaps=True,
            ygap=3,
            )
    fig_cor.update_layout(yaxis=dict(showgrid=False), xaxis=dict(showgrid=False),
        
        legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1
    ))
    #fig_cor.update_layout(yaxis_tickangle=-45)
    fig_cor.update_layout(xaxis_tickangle=0)
    fig_cor.update_layout(title_text='',height=550)#

    return fig_cor
## FF ##

@app.callback(
    [
        Output("max_name", "children"),
        Output("max_value", "children"),
        Output("min_name", "children"),
        Output("min_value", "children"),
        Output("mean", "children"),
        Output("st_dev", "children"),
    ],
    [
        Input("nutrition_types", "value"),
    ]
)
def indicator(auswahl):
    max_id = df[auswahl].idxmax()
    min_id = df[auswahl].idxmin()

    max_value = df[auswahl].max()
    max_value = str(max_value)

    max_name = df.loc[max_id, 'Country']
    min_value = df[auswahl].min()
    min_value = str(min_value)

    min_name = df.loc[min_id, 'Country']
    mean = df[auswahl].mean()
    st_dev = df[auswahl].std()
    st_dev = round(st_dev, 2)
    st_dev = str(st_dev)
    mean = round(mean, 2)
    mean = str(mean)


    return "Country: " + max_name, max_value + " kg per capita per year",  \
           "Country: " + min_name, min_value + " kg per capita per year", \
           mean + " kg per capita per year", \
            st_dev + " kg per capita per year"


@app.callback(
    Output('indicator-graphic', 'figure'),
    Input('xaxis-column', 'value'),
    Input('yaxis-column', 'value'),
    Input('xaxis-type', 'value'),
    Input('yaxis-type', 'value'))
def update_graph(xaxis_column_name, yaxis_column_name, xaxis_type, yaxis_type):

    #col_name = str(yaxis_column_name) + " (above Average)"
    col_name = " "
    df_scatter[col_name] = (df_scatter[yaxis_column_name] > df_scatter[yaxis_column_name].mean())

    def aa(inp):
        if inp == True:
            return yaxis_column_name + " above average"
        else:
            return yaxis_column_name + " below average"

    df_scatter[col_name] = df_scatter[col_name].apply(aa)

    if yaxis_type == "Box":
        type_y = "box"
    else:
        type_y = "violin"
    

    if xaxis_type == "Box":
        type_x = "box"
    else:
        type_x = "violin"

    fig = px.scatter(df_scatter, x=xaxis_column_name, y =yaxis_column_name, size="GDP per Capita", color=col_name, hover_name="Country", 
                    log_x=False,marginal_x = type_x,marginal_y = type_y, template="simple_white", color_discrete_sequence=["#0d0887", "#9c179e"])



    # linear regression
    regline = sm.OLS(df_scatter[yaxis_column_name],sm.add_constant(df_scatter[xaxis_column_name])).fit().fittedvalues

    # add linear regression line for whole sample
    fig.add_traces(go.Scatter(x=df_scatter[xaxis_column_name], y=regline,
                            mode = 'lines',
                            marker_color='#fb9f3a',
                            name='OLS Trendline')
                            )

    fig.update_layout(legend=dict(orientation="h",xanchor='center',x=0.5,yanchor='top',y=-0.2))

    fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')

    #fig.update_xaxes(title=xaxis_column_name,
     #                type='linear' if xaxis_type == 'Linear' else 'log')

    #fig.update_yaxes(title=yaxis_column_name,
     #                type='linear' if yaxis_type == 'Linear' else 'log')

    return fig


server = app.server

if __name__ == '__main__':
    app.run_server(debug=True)
