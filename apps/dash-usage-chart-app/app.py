# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.figure_factory as ff
import numpy as np
import pathlib

from dash.dependencies import Input, Output, State
from datetime import datetime

external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

df = pd.read_csv("data/data.csv")


df["datetime"] = pd.to_datetime(df["datetime"])
df["endtime"] = df["datetime"] + pd.to_timedelta(df["time_on_page"], unit="s")
df.sort_values(by=["datetime"], inplace=True)
df.reset_index(inplace=True)

active_sessions = []
ganttData = []
layers = []
# Loop through each row
for i, r in df.iterrows():
    got_one = False
    layer = 0
    # Loop through active sessions
    for j, end in enumerate(active_sessions):
        # if the session to add happens after the active session
        if r.datetime >= end:
            active_sessions[j] = r.datetime + pd.to_timedelta(r.time_on_page, unit="s")
            got_one = True
            break
        # there is an overlap so make a new layer
        else:
            layer = layer + 1
    # If there is no overlapping time or sessions aren't switched
    if not got_one:
        active_sessions.append(r.datetime + pd.to_timedelta(r.time_on_page, unit="s"))

    # Add to Gantt Data
    ganttData.append(
        dict(
            Task="{}".format(layer + 1),
            Start=str(r.datetime),
            Finish=str(r.endtime),
        )
    )
    layers.append(layer+1)

ganttFig = ff.create_gantt(
                ganttData,
                group_tasks=True,
                showgrid_x=True,
                showgrid_y=True,
                title="Number of Concurrent Users over Time",
                width=1000
            )


app.layout = html.Div(
    children=[
        html.Div(className="header", children=[
            html.Img(
                src=app.get_asset_url('dash-logo-stripe.png'),
                style={'display': 'inline-block',
                      'maxHeight':'65px'}
            ),
            html.H1(children="Usage Chart App", style={'weight':'bold'}),
            html.P("""The Dash Usage Chart App visualizes 
                    number of concurrent users using a certain app over time. 
                    The goal of the app is to help you determine the number of
                    resources (ex. Python workers) you need to service your app.
                    
                    The Gantt Chart allows you to view the number of users using an app
                    at one time. The table below the chart details all the user start and 
                    end time of """)
        ]),
        html.Div(className="row", children=[
                dcc.Graph(
                    id="chart-usage",
                    figure=ganttFig,
                ),
                dash_table.DataTable(
                    id='table',
                    columns=[{"name": "Start", "id": "datetime"},
                            {"name": "End", "id": "endtime"},
                            {"name": "Time on Page (in Seconds)", "id": "time_on_page"}],
                    data=df.to_dict("rows"),
                    style_table={'maxHeight':'500px','overflowY':'scroll'},
                )
            ]
        ),
        
    ]
)


if __name__ == "__main__":
    app.run_server(debug=True)
