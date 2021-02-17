#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Aug31 2020

@author: doug hagey
"""
# This script takes game events and creates a team by team KPI report in the form of a radar graph

import pandas as pd
import plotly.graph_objects as go

# import dash
# import dash_core_components as dcc
# import dash_html_components as html

# Disable non-applicable warnings
pd.options.mode.chained_assignment = None  # default='warn'

# 'Normalize' values (not in a traditional sense, but in terms of percentages during the match)
def normalize_events(df):
    # Rename and trim down the df to columns we want
    df = df.set_index("Team")
    df = df[
        [
            "SET PIECE",
            "PASS",
            "BALL LOST",
            "RECOVERY",
            "CHALLENGE",
            "SHOT",
            "INTERCEPTION",
            "CROSS",
            "DEEP BALL",
            "FREE KICK",
        ]
    ]
    df.columns = [
        "Set Plays",
        "Passes",
        "Balls Lost",
        "Recoveries",
        "Challenges",
        "Shots",
        "Interceptions",
        "Crosses",
        "Long Balls",
        "Free Kicks",
    ]

    # 'Normalize' values from different scales by scaling them as a percentage of max. Use try/except to
    # avoid division by zero
    try:
        df_normalized = (df) / ((df.max() + df.min()))
    except:
        df_normalized = (df + 0.1) / ((df.max() + df.min()))

    # Replace nan's with zeros
    standard_df_normalized = df_normalized.fillna(0)
    return df_normalized


# Main function
def team_radar_builder(filename, team_id):
    # Read in the events data file
    data_file = "data/" + filename
    events_df = pd.read_csv(data_file)
    events_df = events_df[["Team", "Type", "Subtype"]]

    # Count up the events for each team and create a new pivoted dataframe to hold the results
    type_counted_df = (
        events_df.groupby(["Team", "Type"]).size().to_frame("count").reset_index()
    )
    subtype_counted_df = (
        events_df.groupby(["Team", "Subtype"]).size().to_frame("count").reset_index()
    )
    df1_pivoted = pd.pivot_table(
        type_counted_df, values="count", index="Team", columns=["Type"], aggfunc=sum
    ).reset_index()
    df2_pivoted = pd.pivot_table(
        subtype_counted_df,
        values="count",
        index="Team",
        columns=["Subtype"],
        aggfunc=sum,
    ).reset_index()
    df_pivoted = pd.merge(df1_pivoted, df2_pivoted, on="Team")
    df_pivoted = df_pivoted.fillna(0)
    normalized_df = normalize_events(df_pivoted)

    theta_values = list(normalized_df.columns)

    # Create initial figure
    fig = go.Figure()
    colormap = {"Home": "lightskyblue", "Away": "mintcream"}
    normalized_df["Team"] = normalized_df.index
    team_row_normalized = normalized_df.loc[normalized_df["Team"] == team_id]
    team_row_normalized.drop("Team", axis=1, inplace=True)
    team_row = team_row_normalized.reset_index(drop=True)
    r_values = team_row.iloc[0].values.flatten().tolist()

    # Providing a way to use normalized values for the graph but actual values for hover
    team_row_value = df_pivoted.loc[df_pivoted["Team"] == team_id]
    team_row_value.drop("Team", axis=1, inplace=True)
    team_row_value = team_row_value[
        [
            "SET PIECE",
            "PASS",
            "BALL LOST",
            "RECOVERY",
            "CHALLENGE",
            "SHOT",
            "INTERCEPTION",
            "CROSS",
            "DEEP BALL",
            "FREE KICK",
        ]
    ]
    team_row_value = team_row_value.reset_index(drop=True)
    pki_values = team_row_value.iloc[0].values.flatten().tolist()

    fig.add_trace(
        go.Scatterpolar(
            r=r_values,
            theta=theta_values,
            fill="toself",
            name=team_id,
            opacity=0.25,
            fillcolor=colormap[team_id],
            hovertext=pki_values,
            hovertemplate="%{theta}:" + " %{hovertext}<br>",
        )
    )

    # Making some small markers that users can hover over to get more info (else they might not know where to hover)
    fig.update_traces(
        mode="lines+markers",
        line_color=colormap[team_id],
        marker=dict(color="white", symbol="circle", size=4),
    )

    fig.update_layout(
        title="Match KPIs",
        polar=dict(
            bgcolor="#2A2A2A",
            radialaxis=dict(visible=True, range=[0, 1], showticklabels=False,),
        ),
        showlegend=False,
        autosize=False,
    )

    fig.update_layout(
        template="plotly_dark",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
    )

    fig.update_layout(modebar=dict(bgcolor="rgba(0, 0, 0, 0)"))
    fig.update_layout(margin=dict(l=55, r=55, b=30, t=45))

    return fig


# This part is here in case anybody wants to execute it standalone without Dash, especially for
# troubleshooting purposes
"""if __name__ == "__main__":
    filename= input('Filename: ')
    team = input('Team Id: ')
    fig = team_radar_builder(filename, team)

    app = dash.Dash()
    app.layout = html.Div([
        dcc.Graph(figure=fig)
    ])

    app.run_server(debug=False, use_reloader=False)"""
