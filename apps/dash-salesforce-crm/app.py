from dash import Dash, html, dcc, Input, Output
import dash_bootstrap_components as dbc
from panels import opportunities, cases, leads

from constants import salesforce_manager
from utils.components import Header

app = Dash(
    __name__, 
    external_stylesheets=[dbc.themes.CYBORG],
    title="CRM Salesforce"
)
server = app.server

app.layout = dbc.Container([
        Header(app),
        dbc.Tabs([
            dbc.Tab(opportunities.layout, label="Opportunities"),
            dbc.Tab(leads.layout, label="Leads"),
            dbc.Tab(cases.layout, label="Cases"),
        ]),
        
        dcc.Store(id="opportunities_df", data=salesforce_manager.get_opportunities().to_json(orient="split")),
        dcc.Store(id="leads_df", data=salesforce_manager.get_leads().to_json(orient="split")),
        dcc.Store(id="cases_df", data=salesforce_manager.get_cases().to_json(orient="split")),
    ],
    fluid=True,
)

if __name__ == "__main__":
    app.run_server(debug=True)
