from dash import html, dcc, Input, Output, State, callback
import dash_bootstrap_components as dbc
import pandas as pd

from constants import salesforce_manager
from utils.components import cases_modal, cases_controls, cases_data_cards, cases_graphs
from utils.graphs import cases_by_account, cases_by_period, cases_pie_chart


layout = html.Div([
    dbc.Row(cases_controls),
    dbc.Row(cases_data_cards),
    dbc.Row(cases_graphs),
    cases_modal(salesforce_manager),
])

@callback(
    Output("left_cases_indicator", "children"), 
    Output("middle_cases_indicator", "children"),
    Output("right_cases_indicator", "children"),
    Output("cases_by_account", "figure"),
    Output("cases_reasons", "figure"),
    Output("cases_types", "figure"),
    Output("cases_by_period", "figure"),
    Input("cases_df", "data"),
    Input("cases_period_dropdown", "value"),
    Input("origin_dropdown", "value"),
    Input("priority_dropdown", "value"),
)
def update_graphs(df, period, origin, priority):
    df = pd.read_json(df, orient="split")

    ## Data cards
    left_cases_indicator = len(df[(df["Priority"] == "Low") & (df["Status"] == "New")]["Priority"].index)
    middle_cases_indicator = len(df[(df["Priority"] == "Medium") & (df["Status"] == "New")]["Priority"].index)
    right_cases_indicator = len(df[(df["Priority"] == "High") & (df["Status"] == "New")]["Priority"].index)

    ## Figures
    fig_by_account = cases_by_account(df, salesforce_manager)
    fig_pie_chart = cases_pie_chart(df, "Reason", priority, origin)
    fig_pie_chart_h = cases_pie_chart(df, "Type", priority, origin, h_orientation=True)
    fig_by_period = cases_by_period(df, period, priority)

    return left_cases_indicator, middle_cases_indicator, right_cases_indicator, \
            fig_by_account, fig_pie_chart, fig_pie_chart_h, fig_by_period

@callback(
    Output("cases_modal", "is_open"),
    Input("new_case", "n_clicks"), 
    State("cases_modal", "is_open"),
)
def toggle_cases_modal(n1, is_open):
    if n1:
        return not is_open
    return is_open

@callback(
    Output("cases_df", "data"),
    Input("submit_new_case", "n_clicks"),
    State("new_case_account", "value"),
    State("new_case_origin", "value"),
    State("new_case_reason", "value"),
    State("new_case_subject", "value"),
    State("new_case_contact", "value"),
    State("new_case_type", "value"),
    State("new_case_status", "value"),
    State("new_case_description", "value"),
    State("new_case_priority", "value"),
    State("cases_df", "data"),
    prevent_initial_call=True,
)
def add_new_case(n_clicks, account_id, origin, reason, subject, contact_id, case_type, status, description, priority, current_df):
    if n_clicks > 0:
        query = {
            "AccountId": account_id,
            "Origin": origin,
            "Reason": reason,
            "Subject": subject,
            "ContactId": contact_id,
            "Type": case_type,
            "Status": status,
            "Description": description,
            "Priority": priority,
        }

        salesforce_manager.add_case(query)
        df = salesforce_manager.get_cases()
        return df.to_json(orient="split")

    return current_df
