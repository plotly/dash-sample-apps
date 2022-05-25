from dash import html, Input, Output, State, callback
import dash_bootstrap_components as dbc
import pandas as pd

from constants import salesforce_manager
from utils.components import leads_modal, leads_controls, leads_data_cards, leads_graphs
import utils.figures as figs

layout = html.Div([
    dbc.Row(leads_controls),
    dbc.Row(leads_data_cards),
    dbc.Row(leads_graphs),
    leads_modal(),
])

@callback(
    Output("left_leads_indicator", "children"), 
    Output("middle_leads_indicator", "children"),
    Output("right_leads_indicator", "children"),
    Output("lead_source", "figure"),
    Output("leads_map", "figure"),
    Output("converted_leads", "figure"),
    Output("leads_table", "data"),
    Output("leads_table", "columns"),
    Input("leads_df", "data"),
    Input("lead_source_dropdown", "value"),
    Input("converted_leads_dropdown", "value"),
)
def update_graphs(df, status, period):
    df = pd.read_json(df, orient="split")

    ## Data Cards
    left_leads_indicator = len(df[df["Status"] == "Closed - Converted"].index)
    middle_leads_indicator = len( df[ (df["Status"] == "Open - Not Contacted") | (df["Status"] == "Working - Contacted") ].index)

    lost_leads = len(df[df["Status"] == "Closed - Not Converted"].index)
    right_leads_indicator = left_leads_indicator / (left_leads_indicator + lost_leads) * 100

    ## Figures
    fig_lead_source = figs.lead_source(status, df)
    fig_leads_choropleth_map = figs.leads_choropleth_map(status, df)
    fig_converted_leads_count = figs.converted_leads_count(period, df)

    ## Table
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
    table_data = df.to_dict('records')
    table_cols = [{"name": i, "id": i} for i in df.columns]

    return left_leads_indicator, middle_leads_indicator, f"{right_leads_indicator:.2f}%", \
            fig_lead_source, fig_leads_choropleth_map, fig_converted_leads_count, \
            table_data, table_cols

@callback(
    Output("leads_modal", "is_open"),
    Input("new_lead", "n_clicks"), 
    State("leads_modal", "is_open"),
)
def toggle_cases_modal(n1, is_open):
    if n1:
        return not is_open
    return is_open


@callback(
    Output("leads_df", "data"),
    Input("submit_new_lead", "n_clicks"),
    State("new_lead_status", "value"),
    State("new_lead_state", "value"),
    State("new_lead_company", "value"),
    State("new_lead_source", "value"),
    State("leads_df", "data"),
    prevent_initial_call=True,
)
def add_new_lead(n_clicks, status, state, company, source, current_df):
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
        salesforce_manager.add_lead(query)
        df = salesforce_manager.get_leads()
        return df.to_json(orient="split")

    return current_df
