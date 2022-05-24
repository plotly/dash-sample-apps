from dash import html, dcc, Input, Output, State, callback
import dash_bootstrap_components as dbc
import pandas as pd

from constants import salesforce_manager
from utils.helper_functions import millify, top_open_opportunities, top_lost_opportunities
from utils.components import opportunities_modal, opportunities_controls, opportunities_data_cards, opportunities_graphs
from utils.graphs import converted_opportunities, opportunities_heat_map_fig

layout = html.Div([
    dbc.Row(opportunities_controls),
    dbc.Row(opportunities_data_cards),
    dbc.Row(opportunities_graphs),
    opportunities_modal(),
])


@callback(
    Output("left_opportunities_indicator", "children"),
    Output("middle_opportunities_indicator", "children"),
    Output("right_opportunities_indicator", "children"),
    Output("opportunities_heatmap", "figure"),
    Output("converted_count", "figure"),
    Output("top_open_opportunities", "data"), 
    Output("top_open_opportunities", "columns"), 
    Output("top_lost_opportunities", "data"), 
    Output("top_lost_opportunities", "columns"), 
    Input("opportunities_df", "data"),
    Input("heatmap_dropdown", "value"),
    Input("converted_opportunities_dropdown", "value"),
    Input("source_dropdown", "value"),
)
def update_graphs(df, stage, period, source):
    df = pd.read_json(df, orient="split")

    ## Data Cards
    left_opportunities_indicator = millify(str(df[df["IsWon"] == 1]["Amount"].sum()))
    middle_opportunities_indicator = millify(str(df[(df["IsClosed"] == 0)]["Amount"].sum()))
    right_opportunities_indicator = millify(str(df[(df["IsWon"] == 0) & (df["IsClosed"] == 1)]["Amount"].sum()))
    
    ## Figures
    fig_converted_opportunities = converted_opportunities(period, source, df)
    fig_heatmap = opportunities_heat_map_fig(df, stage)

    ## Table
    table_data_top_ten, table_cols_top_ten = top_open_opportunities(df)
    table_data_top_lost, table_cols_top_lost =  top_lost_opportunities(df)

    return left_opportunities_indicator, middle_opportunities_indicator, right_opportunities_indicator, \
            fig_heatmap, fig_converted_opportunities, \
            table_data_top_ten, table_cols_top_ten, table_data_top_lost, table_cols_top_lost

@callback(
    Output("opportunities_modal", "is_open"),
    Input("new_opportunity", "n_clicks"), 
    State("opportunities_modal", "is_open"),
)
def toggle_opportunities_modal(n1, is_open):
    if n1:
        return not is_open
    return is_open


@callback(
    Output("opportunities_df", "data"),
    Input("submit_new_opportunity", "n_clicks"),
    State("new_opportunity_name", "value"),
    State("new_opportunity_stage", "value"),
    State("new_opportunity_amount", "value"),
    State("new_opportunity_probability", "value"),
    State("new_opportunity_date", "date"),
    State("new_opportunity_type", "value"),
    State("new_opportunity_source", "value"),
    State("opportunities_df", "data"),
    prevent_initial_call=True,
)
def add_new_opportunity(
    n_clicks, name, stage, amount, probability, date, o_type, source, current_df
):
    if n_clicks > 0:
        if name == "":
            name = "Not named yet"
        query = {
            "Name": name,
            "StageName": stage,
            "Amount": amount,
            "Probability": probability,
            "CloseDate": date,
            "Type": o_type,
            "LeadSource": source,
        }
        salesforce_manager.add_opportunity(query)
        df = salesforce_manager.get_opportunities()
        return df.to_json(orient="split")

    return current_df