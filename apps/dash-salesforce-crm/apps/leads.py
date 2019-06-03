# -*- coding: utf-8 -*-
import json
import math

import pandas as pd
import flask
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import plotly.plotly as py
from plotly import graph_objs as go

from app import app, indicator, millify, df_to_table, sf_manager

states = ["AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DC", "DE", "FL", "GA", 
          "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", 
          "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", 
          "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", 
          "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"]



# returns choropleth map figure based on status filter
def choropleth_map(status, df):
    if status == "open":
        df = df[
            (df["Status"] == "Open - Not Contacted")
            | (df["Status"] == "Working - Contacted")
        ]

    elif status == "converted":
        df = df[df["Status"] == "Closed - Converted"]

    elif status == "lost":
        df = df[df["Status"] == "Closed - Not Converted"]

    df = df.groupby("State").count()
    
    scl = [[0.0, "rgb(38, 78, 134)"], [1.0, "#0091D5"]] # colors scale

    data = [
        dict(
            type="choropleth",
            colorscale=scl,
            locations=df.index,
            z=df["Id"],
            locationmode="USA-states",
            marker=dict(line=dict(color="rgb(255,255,255)", width=2)),
        )
    ]

    layout = dict(
        geo=dict(
            scope="usa",
            projection=dict(type="albers usa"),
            lakecolor="rgb(255, 255, 255)",
        ),
        margin=dict(l=10, r=10, t=0, b=0),
    )
    return dict(data=data, layout=layout)


# returns pie chart that shows lead source repartition
def lead_source(status, df):
    if status == "open":
        df = df[
            (df["Status"] == "Open - Not Contacted")
            | (df["Status"] == "Working - Contacted")
        ]

    elif status == "converted":
        df = df[df["Status"] == "Closed - Converted"]

    elif status == "lost":
        df = df[df["Status"] == "Closed - Not Converted"]

    nb_leads = len(df.index)
    types = df["LeadSource"].unique().tolist()
    values = []

    # compute % for each leadsource type
    for case_type in types:
        nb_type = df[df["LeadSource"] == case_type].shape[0]
        values.append(nb_type / nb_leads * 100)

    trace = go.Pie(
        labels=types,
        values=values,
        marker={"colors": ["#264e86", "#0074e4", "#74dbef", "#eff0f4"]},
    )

    layout = dict(margin=dict(l=15, r=10, t=0, b=65), legend=dict(orientation="h"))
    return dict(data=[trace], layout=layout)



def converted_leads_count(period, df):
    df["CreatedDate"] = pd.to_datetime(df["CreatedDate"], format="%Y-%m-%d")
    df = df[df["Status"] == "Closed - Converted"]

    df = (
        df.groupby([pd.Grouper(key="CreatedDate", freq=period)])
        .count()
        .reset_index()
        .sort_values("CreatedDate")
    )

    trace = go.Scatter(
        x=df["CreatedDate"],
        y=df["Id"],
        name="converted leads",
        fill="tozeroy",
        fillcolor="#e6f2ff",
    )

    data = [trace]

    layout = go.Layout(
        xaxis=dict(showgrid=False),
        margin=dict(l=33, r=25, b=37, t=5, pad=4),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return {"data": data, "layout": layout}


def modal():
    return html.Div(
        html.Div(
            [
                html.Div(
                    [   

                        # modal header
                        html.Div(
                            [
                                html.Span(
                                    "New Lead",
                                    style={
                                        "color": "#506784",
                                        "fontWeight": "bold",
                                        "fontSize": "20",
                                    },
                                ),
                                html.Span(
                                    "×",
                                    id="leads_modal_close",
                                    n_clicks=0,
                                    style={
                                        "float": "right",
                                        "cursor": "pointer",
                                        "marginTop": "0",
                                        "marginBottom": "17",
                                    },
                                ),
                            ],
                            className="row",
                            style={"borderBottom": "1px solid #C8D4E3"},
                        ),

                        # modal form
                        html.Div(
                            [
                                html.P(
                                    [
                                        "Company Name",
                                        
                                    ],
                                    style={
                                        "float": "left",
                                        "marginTop": "4",
                                        "marginBottom": "2",
                                    },
                                    className="row",
                                ),
                                dcc.Input(
                                    id="new_lead_company",
                                    # placeholder="Enter company name",
                                    type="text",
                                    value="",
                                    style={"width": "100%"},
                                ),
                                html.P(
                                    "Company State",
                                    style={
                                        "textAlign": "left",
                                        "marginBottom": "2",
                                        "marginTop": "4",
                                    },
                                ),
                                dcc.Dropdown(
                                    id="new_lead_state",
                                    options=[
                                        {"label": state, "value": state}
                                        for state in states
                                    ],
                                    value="NY",
                                ),
                                html.P(
                                    "Status",
                                    style={
                                        "textAlign": "left",
                                        "marginBottom": "2",
                                        "marginTop": "4",
                                    },
                                ),
                                dcc.Dropdown(
                                    id="new_lead_status",
                                    options=[
                                        {
                                            "label": "Open - Not Contacted",
                                            "value": "Open - Not Contacted",
                                        },
                                        {
                                            "label": "Working - Contacted",
                                            "value": "Working - Contacted",
                                        },
                                        {
                                            "label": "Closed - Converted",
                                            "value": "Closed - Converted",
                                        },
                                        {
                                            "label": "Closed - Not Converted",
                                            "value": "Closed - Not Converted",
                                        },
                                    ],
                                    value="Open - Not Contacted",
                                ),
                                html.P(
                                    "Source",
                                    style={
                                        "textAlign": "left",
                                        "marginBottom": "2",
                                        "marginTop": "4",
                                    },
                                ),
                                dcc.Dropdown(
                                    id="new_lead_source",
                                    options=[
                                        {"label": "Web", "value": "Web"},
                                        {
                                            "label": "Phone Inquiry",
                                            "value": "Phone Inquiry",
                                        },
                                        {
                                            "label": "Partner Referral",
                                            "value": "Partner Referral",
                                        },
                                        {
                                            "label": "Purchased List",
                                            "value": "Purchased List",
                                        },
                                        {"label": "Other", "value": "Other"},
                                    ],
                                    value="Web",
                                ),
                            ],
                            className="row",
                            style={"padding": "2% 8%"},
                        ),

                        # submit button
                        html.Span(
                            "Submit",
                            id="submit_new_lead",
                            n_clicks=0,
                            className="button button--primary add"
                        ),
                    ],
                    className="modal-content",
                    style={"textAlign": "center"},
                )
            ],
            className="modal",
        ),
        id="leads_modal",
        style={"display": "none"},
    )


layout = [

    # top controls
    html.Div(
        [
            html.Div(
                dcc.Dropdown(
                    id="converted_leads_dropdown",
                    options=[
                        {"label": "By day", "value": "D"},
                        {"label": "By week", "value": "W-MON"},
                        {"label": "By month", "value": "M"},
                    ],
                    value="D",
                    clearable=False,
                ),
                className="two columns",
            ),
            html.Div(
                dcc.Dropdown(
                    id="lead_source_dropdown",
                    options=[
                        {"label": "All status", "value": "all"},
                        {"label": "Open leads", "value": "open"},
                        {"label": "Converted leads", "value": "converted"},
                        {"label": "Lost leads", "value": "lost"},
                    ],
                    value="all",
                    clearable=False,
                ),
                className="two columns",
            ),

            # add button
            html.Div(
                html.Span(
                    "Add new",
                    id="new_lead",
                    n_clicks=0,
                    className="button button--primary",
                    style={
                        "height": "34",
                        "background": "#119DFF",
                        "border": "1px solid #119DFF",
                        "color": "white",
                    },
                ),
                className="two columns",
                style={"float": "right"},
            ),
        ],
        className="row",
        style={"marginBottom": "10"},
    ),

    # indicators row div
    html.Div(
        [
            indicator(
                "#00cc96", "Converted Leads", "left_leads_indicator"
            ),
            indicator(
                "#119DFF", "Open Leads", "middle_leads_indicator"
            ),
            indicator(
                "#EF553B",
                "Conversion Rates",
                "right_leads_indicator",
            ),
        ],
        className="row",
    ),

    # charts row div
    html.Div(
        [
            html.Div(
                [
                    html.P("Leads count per state" ),
                    dcc.Graph(
                        id="map",
                        style={"height": "90%", "width": "98%"},
                        config=dict(displayModeBar=False),
                    ),
                ],
                className="four columns chart_div"
            ),

            html.Div(
                [
                    html.P("Leads by source"),
                    dcc.Graph(
                        id="lead_source",
                        style={"height": "90%", "width": "98%"},
                        config=dict(displayModeBar=False),
                    ),
                ],
                className="four columns chart_div"
            ),

            html.Div(
                [
                    html.P("Converted Leads count"),
                    dcc.Graph(
                        id="converted_leads",
                        style={"height": "90%", "width": "98%"},
                        config=dict(displayModeBar=False),
                    ),
                ],
                className="four columns chart_div"
            ),
        ],
        className="row",
        style={"marginTop": "5"},
    ),

    # table div
    html.Div(
        id="leads_table",
        className="row",
        style={
            "maxHeight": "350px",
            "overflowY": "scroll",
            "padding": "8",
            "marginTop": "5",
            "backgroundColor":"white",
            "border": "1px solid #C8D4E3",
            "borderRadius": "3px"
        },
    ),


    modal(),
]


# updates left indicator based on df updates
@app.callback(
    Output("left_leads_indicator", "children"), [Input("leads_df", "children")]
)
def left_leads_indicator_callback(df):
    df = pd.read_json(df, orient="split")
    converted_leads = len(df[df["Status"] == "Closed - Converted"].index)
    return converted_leads


# updates middle indicator based on df updates
@app.callback(
    Output("middle_leads_indicator", "children"), [Input("leads_df", "children")]
)
def middle_leads_indicator_callback(df):
    df = pd.read_json(df, orient="split")
    open_leads = len(
        df[
            (df["Status"] == "Open - Not Contacted")
            | (df["Status"] == "Working - Contacted")
        ].index
    )
    return open_leads


# updates right indicator based on df updates
@app.callback(
    Output("right_leads_indicator", "children"), [Input("leads_df", "children")]
)
def right_leads_indicator_callback(df):
    df = pd.read_json(df, orient="split")
    converted_leads = len(df[df["Status"] == "Closed - Converted"].index)
    lost_leads = len(df[df["Status"] == "Closed - Not Converted"].index)
    conversion_rates = converted_leads / (converted_leads + lost_leads) * 100
    conversion_rates = "%.2f" % conversion_rates + "%"
    return conversion_rates


# update pie chart figure based on dropdown's value and df updates
@app.callback(
    Output("lead_source", "figure"),
    [Input("lead_source_dropdown", "value"), Input("leads_df", "children")],
)
def lead_source_callback(status, df):
    df = pd.read_json(df, orient="split")
    return lead_source(status, df)


# update heat map figure based on dropdown's value and df updates
@app.callback(
    Output("map", "figure"),
    [Input("lead_source_dropdown", "value"), Input("leads_df", "children")],
)
def map_callback(status, df):
    df = pd.read_json(df, orient="split")
    return choropleth_map(status, df)


# update table based on dropdown's value and df updates
@app.callback(
    Output("leads_table", "children"),
    [Input("lead_source_dropdown", "value"), Input("leads_df", "children")],
)
def leads_table_callback(status, df):
    df = pd.read_json(df, orient="split")
    if status == "open":
        df = df[
            (df["Status"] == "Open - Not Contacted")
            | (df["Status"] == "Working - Contacted")
        ]
    elif status == "converted":
        df = df[df["Status"] == "Closed - Converted"]
    elif status == "lost":
        df = df[df["Status"] == "Closed - Not Converted"]
    df = df[["CreatedDate", "Status", "Company", "State", "LeadSource"]]
    return df_to_table(df)


# update pie chart figure based on dropdown's value and df updates
@app.callback(
    Output("converted_leads", "figure"),
    [Input("converted_leads_dropdown", "value"), Input("leads_df", "children")],
)
def converted_leads_callback(period, df):
    df = pd.read_json(df, orient="split")
    return converted_leads_count(period, df)


# hide/show modal
@app.callback(Output("leads_modal", "style"), [Input("new_lead", "n_clicks")])
def display_leads_modal_callback(n):
    if n > 0:
        return {"display": "block"}
    return {"display": "none"}


# reset to 0 add button n_clicks property 
@app.callback(
    Output("new_lead", "n_clicks"),
    [Input("leads_modal_close", "n_clicks"), Input("submit_new_lead", "n_clicks")],
)
def close_modal_callback(n, n2):
    return 0


# add new lead to salesforce and stores new df in hidden div
@app.callback(
    Output("leads_df", "children"),
    [Input("submit_new_lead", "n_clicks")],
    [
        State("new_lead_status", "value"),
        State("new_lead_state", "value"),
        State("new_lead_company", "value"),
        State("new_lead_source", "value"),
        State("leads_df", "children"),
    ],
)
def add_lead_callback(n_clicks, status, state, company, source, current_df):
    if n_clicks > 0:
        if company == "":
            company = "Not named yet"
        query = {
            "LastName": company,
            "Company": company,
            "Status": status,
            "State": state,
            "LeadSource": source,
        }
        sf_manager.add_lead(query)
        df = sf_manager.get_leads()
        return df.to_json(orient="split")

    return current_df
