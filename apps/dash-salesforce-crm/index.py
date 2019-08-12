import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from app import sf_manager, app
from panels import opportunities, cases, leads


server = app.server

app.layout = html.Div(
    [
        html.Div(
            className="row header",
            children=[
                html.Button(id="menu", children=dcc.Markdown("&#8801")),
                html.Span(
                    className="app-title",
                    children=[
                        dcc.Markdown("**CRM App**"),
                        html.Span(
                            id="subtitle",
                            children=dcc.Markdown("&nbsp using Salesforce API"),
                            style={"font-size": "1.8rem", "margin-top": "15px"},
                        ),
                    ],
                ),
                html.Img(src=app.get_asset_url("logo.png")),
                html.A(
                    id="learn_more",
                    children=html.Button("Learn More"),
                    href="https://plot.ly/dash/",
                ),
            ],
        ),
        html.Div(
            id="tabs",
            className="row tabs",
            children=[
                dcc.Link("Opportunities", href="/"),
                dcc.Link("Leads", href="/"),
                dcc.Link("Cases", href="/"),
            ],
        ),
        html.Div(
            id="mobile_tabs",
            className="row tabs",
            style={"display": "none"},
            children=[
                dcc.Link("Opportunities", href="/"),
                dcc.Link("Leads", href="/"),
                dcc.Link("Cases", href="/"),
            ],
        ),
        dcc.Store(  # opportunities df
            id="opportunities_df",
            data=sf_manager.get_opportunities().to_json(orient="split"),
        ),
        dcc.Store(  # leads df
            id="leads_df", data=sf_manager.get_leads().to_json(orient="split")
        ),
        dcc.Store(
            id="cases_df", data=sf_manager.get_cases().to_json(orient="split")
        ),  # cases df
        dcc.Location(id="url", refresh=False),
        html.Div(id="tab_content"),
        html.Link(
            href="https://use.fontawesome.com/releases/v5.2.0/css/all.css",
            rel="stylesheet",
        ),
        html.Link(
            href="https://fonts.googleapis.com/css?family=Dosis", rel="stylesheet"
        ),
        html.Link(
            href="https://fonts.googleapis.com/css?family=Open+Sans", rel="stylesheet"
        ),
        html.Link(
            href="https://fonts.googleapis.com/css?family=Ubuntu", rel="stylesheet"
        ),
    ],
    className="row",
    style={"margin": "0%"},
)

# Update the index


@app.callback(
    [
        Output("tab_content", "children"),
        Output("tabs", "children"),
        Output("mobile_tabs", "children"),
    ],
    [Input("url", "pathname")],
)
def display_page(pathname):
    tabs = [
        dcc.Link("Opportunities", href="/dash-salesforce-crm/opportunities"),
        dcc.Link("Leads", href="/dash-salesforce-crm/leads"),
        dcc.Link("Cases", href="/dash-salesforce-crm/cases"),
    ]
    if pathname == "/dash-salesforce-crm/opportunities":
        tabs[0] = dcc.Link(
            dcc.Markdown("**&#9632 Opportunities**"),
            href="/dash-salesforce-crm/opportunities",
        )
        return opportunities.layout, tabs, tabs
    elif pathname == "/dash-salesforce-crm/cases":
        tabs[2] = dcc.Link(
            dcc.Markdown("**&#9632 Cases**"), href="/dash-salesforce-crm/cases"
        )
        return cases.layout, tabs, tabs
    tabs[1] = dcc.Link(
        dcc.Markdown("**&#9632 Leads**"), href="/dash-salesforce-crm/leads"
    )
    return leads.layout, tabs, tabs


@app.callback(
    Output("mobile_tabs", "style"),
    [Input("menu", "n_clicks")],
    [State("mobile_tabs", "style")],
)
def show_menu(n_clicks, tabs_style):
    if n_clicks:
        if tabs_style["display"] == "none":
            tabs_style["display"] = "flex"
        else:
            tabs_style["display"] = "none"
    return tabs_style


if __name__ == "__main__":
    app.run_server(debug=True)
