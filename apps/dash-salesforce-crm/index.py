import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import flask
import plotly.plotly as py
import pathlib

from app import app, server, sf_manager
from apps import opportunities, cases, leads

app.layout = html.Div(
    [
        html.Div(
            className="row header",
            children=[
                html.Span(
                    dcc.Markdown("**CRM App** using Salesforce API"),
                    className="app-title",
                ),
                html.Img(src=app.get_asset_url("logo.png")),
                html.A(
                    children=html.Button("Learn More"), href="https://plot.ly/dash/"
                ),
            ],
        ),
        html.Div(
            id="tabs",
            className="row tabs_div",
            children=[
                dcc.Link(dcc.Markdown("Opportunities"), href="/opportunities"),
                dcc.Link(dcc.Markdown("Leads"), href="/leads"),
                dcc.Link(dcc.Markdown("Cases"), href="/cases"),
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
        dash.dependencies.Output("tab_content", "children"),
        dash.dependencies.Output("tabs", "children"),
    ],
    [dash.dependencies.Input("url", "pathname")],
)
def display_page(pathname):
    print("display_page({})".format(pathname))
    tabs = [
        dcc.Link(dcc.Markdown("Opportunities"), href="/opportunities"),
        dcc.Link(dcc.Markdown("Leads"), href="/leads"),
        dcc.Link(dcc.Markdown("Cases"), href="/cases"),
    ]
    if pathname == "/opportunities":
        tabs[0] = dcc.Link(
            dcc.Markdown("**&#9632 Opportunities**"), href="/opportunities"
        )
        return opportunities.layout, tabs
    elif pathname == "/cases":
        tabs[2] = dcc.Link(dcc.Markdown("**&#9632 Cases**"), href="/leads")
        return cases.layout, tabs
    tabs[1] = dcc.Link(dcc.Markdown("**&#9632 Leads**"), href="/leads")
    return leads.layout, tabs


# @app.callback(Output("tab_content", "children"), [Input("tabs", "value")])
# def render_content(tab):
#     if tab == "opportunities_tab":
#         return opportunities.layout
#     elif tab == "cases_tab":
#         return cases.layout
#     elif tab == "leads_tab":
#         return leads.layout
#     else:
#         return opportunities.layout


if __name__ == "__main__":
    app.run_server(debug=True)
