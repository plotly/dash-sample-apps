# -*- coding: utf-8 -*-
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import flask
import plotly.plotly as py
from plotly import graph_objs as go
import math
from app import app, server, sf_manager
from apps import opportunities, cases, leads

app.layout = html.Div(
    [
        # header
        html.Div([

            html.Span("CRM App using Salesforce API", className='app-title'),
            
            html.Div(
                html.Img(src='https://s3-us-west-1.amazonaws.com/plotly-tutorials/logo/new-branding/dash-logo-by-plotly-stripe-inverted.png',height="100%")
                ,style={"float":"right","height":"100%"})
            ],
            className="row header"
            ),

        # tabs
        html.Div([

            dcc.Tabs(
                id="tabs",
                style={"height":"20","verticalAlign":"middle"},
                children=[
                    dcc.Tab(label="Opportunities", value="opportunities_tab"),
                    dcc.Tab(label="Leads", value="leads_tab"),
                    dcc.Tab(id="cases_tab",label="Cases", value="cases_tab"),
                ],
                value="leads_tab",
            )

            ],
            className="row tabs_div"
            ),
       
                
        # divs that save dataframe for each tab
        html.Div(
                sf_manager.get_opportunities().to_json(orient="split"),  # opportunities df
                id="opportunities_df",
                style={"display": "none"},
            ),
        html.Div(sf_manager.get_leads().to_json(orient="split"), id="leads_df", style={"display": "none"}), # leads df
        html.Div(sf_manager.get_cases().to_json(orient="split"), id="cases_df", style={"display": "none"}), # cases df



        # Tab content
        html.Div(id="tab_content", className="row", style={"margin": "2% 3%"}),
        
        html.Link(href="https://use.fontawesome.com/releases/v5.2.0/css/all.css",rel="stylesheet"),
        html.Link(href="https://cdn.rawgit.com/plotly/dash-app-stylesheets/2d266c578d2a6e8850ebce48fdb52759b2aef506/stylesheet-oil-and-gas.css",rel="stylesheet"),
        html.Link(href="https://fonts.googleapis.com/css?family=Dosis", rel="stylesheet"),
        html.Link(href="https://fonts.googleapis.com/css?family=Open+Sans", rel="stylesheet"),
        html.Link(href="https://fonts.googleapis.com/css?family=Ubuntu", rel="stylesheet"),
        html.Link(href="https://cdn.rawgit.com/amadoukane96/8a8cfdac5d2cecad866952c52a70a50e/raw/cd5a9bf0b30856f4fc7e3812162c74bfc0ebe011/dash_crm.css", rel="stylesheet")
    ],
    className="row",
    style={"margin": "0%"},
)


@app.callback(Output("tab_content", "children"), [Input("tabs", "value")])
def render_content(tab):
    if tab == "opportunities_tab":
        return opportunities.layout
    elif tab == "cases_tab":
        return cases.layout
    elif tab == "leads_tab":
        return leads.layout
    else:
        return opportunities.layout

if __name__ == "__main__":
    app.run_server(debug=True)
