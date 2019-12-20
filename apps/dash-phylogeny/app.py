# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

from dash.dependencies import Input, Output
from utils import *

app = dash.Dash(__name__)
server = app.server

virus_name = "measles"
species = ["Avian", "Ebola", "Lassa", "Measles", "Mumps", "Zika"]
tree_fig = {}
mapbox_access_token = "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNqdnBvNDMyaTAxYzkzeW5ubWdpZ2VjbmMifQ.TXcBE-xg9BFdV2ocecc_7g"


tree_file, metadata_file, metadata_file_stat = create_paths_file(
    virus_name, level1="", level2="", level3=""
)

# To know the minimum and the maximum values of date for slicer
df_stat_metadata = pd.read_csv(metadata_file_stat)
min_date, max_date = min_max_date(df_stat_metadata)

# create the dictionary of slider
marks_data = slicer(min_date, max_date)
min_max_date_value = [min_date, max_date]

fig = create_tree(virus_name, tree_file, metadata_file, "Country")
tree_fig[tree_file] = fig

fig_map_bubble = create_map_bubble_year(
    virus_name, metadata_file_stat, 2, min_date, max_date
)

fig_curve_line = create_curve_line(df_stat_metadata, virus_name, min_date, max_date)


######################################### MAIN APP #########################################
app.layout = html.Div(
    [
        # Banner display
        html.Div(
            className="header-title",
            children=[
                html.H2(
                    id="title",
                    children="Phylogeny trees and global spread of 6 viruses",
                ),
                html.Div(
                    id="learn_more",
                    children=[
                        html.Img(className="logo", src=app.get_asset_url("logo.png")),
                        html.Button(className="btn", children=["Learn More"]),
                    ],
                ),
            ],
        ),
        html.Div(
            id="grid",
            children=[
                html.Div(
                    id="controls",
                    className="row div-row div-card",
                    children=[
                        html.Div(
                            id="dataset-picker",
                            children=[
                                html.Div(
                                    className="six columns",
                                    children=[
                                        html.H6(children="Dataset"),
                                        dcc.Dropdown(
                                            id="d_virus-name",
                                            options=[
                                                {
                                                    "label": species[i],
                                                    "value": species[i],
                                                }
                                                for i in range(len(species))
                                            ],
                                            value="Measles",
                                        ),
                                        html.Div(id="output-container"),
                                    ],
                                ),
                                # Strain dropdown picker
                                html.Div(
                                    className="row",
                                    children=[
                                        html.Div(
                                            className="four columns",
                                            children=[
                                                html.Div(
                                                    children=[
                                                        html.Div(
                                                            id="controls-container_mumps",
                                                            children=[
                                                                dcc.Dropdown(
                                                                    id="d_mumps",
                                                                    options=[
                                                                        {
                                                                            "label": i,
                                                                            "value": i,
                                                                        }
                                                                        for i in [
                                                                            "global",
                                                                            "na",
                                                                        ]
                                                                    ],
                                                                    value="global",
                                                                )
                                                            ],
                                                            style={"display": "none"},
                                                        ),
                                                        html.Div(
                                                            id="controls-container_dengue",
                                                            children=[
                                                                dcc.Dropdown(
                                                                    id="d_dengue",
                                                                    options=[
                                                                        {
                                                                            "label": i,
                                                                            "value": i,
                                                                        }
                                                                        for i in [
                                                                            "all",
                                                                            "denv1",
                                                                            "denv2",
                                                                            "denv3",
                                                                            "denv4",
                                                                        ]
                                                                    ],
                                                                    value="all",
                                                                )
                                                            ],
                                                            style={"display": "none"},
                                                        ),
                                                        html.Div(
                                                            id="controls-container_lassa",
                                                            children=[
                                                                dcc.Dropdown(
                                                                    id="d_lassa",
                                                                    options=[
                                                                        {
                                                                            "label": i,
                                                                            "value": i,
                                                                        }
                                                                        for i in [
                                                                            "s",
                                                                            "l",
                                                                        ]
                                                                    ],
                                                                    value="s",
                                                                )
                                                            ],
                                                            style={"display": "none"},
                                                        ),
                                                        html.Div(
                                                            id="controls-container_avian",
                                                            children=[
                                                                dcc.Dropdown(
                                                                    id="d_avian_opt1",
                                                                    options=[
                                                                        {
                                                                            "label": i,
                                                                            "value": i,
                                                                        }
                                                                        for i in [
                                                                            "h7n9"
                                                                        ]
                                                                    ],
                                                                    value="h7n9",
                                                                ),
                                                                dcc.Dropdown(
                                                                    id="d_avian_opt2",
                                                                    options=[
                                                                        {
                                                                            "label": i,
                                                                            "value": i,
                                                                        }
                                                                        for i in [
                                                                            "ha",
                                                                            "mp",
                                                                            "na",
                                                                            "ns",
                                                                            "np",
                                                                            "pa",
                                                                            "pb2",
                                                                            "pb1",
                                                                        ]
                                                                    ],
                                                                    value="ha",
                                                                ),
                                                            ],
                                                            style={"display": "none"},
                                                        ),
                                                        html.Div(
                                                            id="controls-container_flu",
                                                            children=[
                                                                dcc.Dropdown(
                                                                    id="d_flu_opt1",
                                                                    options=[
                                                                        {
                                                                            "label": i,
                                                                            "value": i,
                                                                        }
                                                                        for i in [
                                                                            "h3n2",
                                                                            "h1n1pdm",
                                                                            "vic",
                                                                            "yam",
                                                                        ]
                                                                    ],
                                                                    value="h3n2",
                                                                ),
                                                                dcc.Dropdown(
                                                                    id="d_flu_opt2",
                                                                    options=[
                                                                        {
                                                                            "label": i,
                                                                            "value": i,
                                                                        }
                                                                        for i in [
                                                                            "ha",
                                                                            "na",
                                                                        ]
                                                                    ],
                                                                    value="ha",
                                                                ),
                                                                dcc.Dropdown(
                                                                    id="d_flu_opt3",
                                                                    options=[
                                                                        {
                                                                            "label": i,
                                                                            "value": i,
                                                                        }
                                                                        for i in [
                                                                            "2y",
                                                                            "3y",
                                                                            "6y",
                                                                            "12y",
                                                                        ]
                                                                    ],
                                                                    value="3y",
                                                                ),
                                                            ],
                                                            style={"display": "none"},
                                                        ),
                                                    ]
                                                )
                                            ],
                                        )
                                    ],
                                ),
                            ],
                        ),
                        html.Div(
                            className="six columns",
                            children=[
                                html.H6(children="Data Range"),
                                html.Div(
                                    id="id-slicer",
                                    children=[
                                        dcc.RangeSlider(
                                            id="id-year",
                                            min=min_date,
                                            max=max_date,
                                            step=1,
                                            marks=marks_data,
                                            value=min_max_date_value,
                                            vertical=True,
                                        )
                                    ],
                                ),
                            ],
                        ),
                    ],
                ),
                dcc.Graph(
                    id="curve-line-graph", className="div-card", figure=fig_curve_line
                ),
                dcc.Graph(id="phylogeny-graph", className="div-card", figure=fig),
                dcc.Graph(id="map-graph", className="div-card", figure=fig_map_bubble),
                dcc.Graph(id="histo-graph", className="div-card"),
            ],
        ),
    ]
)


######################################### UPDATING FIGURES #########################################
@app.callback(Output("output-container", "children"), [Input("d_virus-name", "value")])
def _update_legend_gene(virus_name):
    return 'You have selected "{}" virus'.format(virus_name)


@app.callback(
    Output("controls-container_mumps", "style"), [Input("d_virus-name", "value")]
)
def _update_mumps_option(virus_name):
    if virus_name == "Mumps":
        return {"display": "block"}
    else:
        return {"display": "none"}


@app.callback(
    Output("controls-container_dengue", "style"), [Input("d_virus-name", "value")]
)
def _update_dengue_option(virus_name):
    if virus_name == "Dengue":
        return {"display": "block"}
    else:
        return {"display": "none"}


@app.callback(
    Output("controls-container_lassa", "style"), [Input("d_virus-name", "value")]
)
def _update_lassa_option(virus_name):
    if virus_name == "Lassa":
        return {"display": "block"}
    else:
        return {"display": "none"}


@app.callback(
    Output("controls-container_avian", "style"), [Input("d_virus-name", "value")]
)
def _update_avian_option(virus_name):
    if virus_name == "Avian":
        return {"display": "block"}
    else:
        return {"display": "none"}


@app.callback(
    Output("controls-container_flu", "style"), [Input("d_virus-name", "value")]
)
def _update_flu_option(virus_name):
    if virus_name == "Flu":
        return {"display": "block"}
    else:
        return {"display": "none"}


# Running the server
if __name__ == "__main__":
    app.run_server(debug=True)
