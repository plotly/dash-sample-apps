import dash
import dash_core_components as dcc
import dash_html_components as html
from app import sf_manager, app
from panels import opportunities, cases, leads

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
                            children=dcc.Markdown(
                                "&nbsp using Salesforce API"),
                            style={"font-size": "1.8rem", "margin-top": "15px"}
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
                dcc.Link(dcc.Markdown("Opportunities"), href="/opportunities"),
                dcc.Link(dcc.Markdown("Leads"), href="/leads"),
                dcc.Link(dcc.Markdown("Cases"), href="/cases"),
            ],
        ),
        html.Div(
            id="mobile_tabs",
            className="row tabs",
            style={"display": "none"},
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
        dash.dependencies.Output("mobile_tabs", "children"),
    ],
    [dash.dependencies.Input("url", "pathname")],
)
def display_page(pathname):
    # print("display_page({})".format(pathname))
    tabs = [
        dcc.Link(dcc.Markdown("Opportunities"), href="/opportunities"),
        dcc.Link(dcc.Markdown("Leads"), href="/leads"),
        dcc.Link(dcc.Markdown("Cases"), href="/cases"),
    ]
    if pathname == "/opportunities":
        tabs[0] = dcc.Link(
            dcc.Markdown("**&#9632 Opportunities**"), href="/opportunities"
        )
        return opportunities.layout, tabs, tabs
    elif pathname == "/cases":
        tabs[2] = dcc.Link(dcc.Markdown("**&#9632 Cases**"), href="/leads")
        return cases.layout, tabs, tabs
    tabs[1] = dcc.Link(dcc.Markdown("**&#9632 Leads**"), href="/leads")
    return leads.layout, tabs, tabs


@app.callback(
    dash.dependencies.Output("mobile_tabs", "style"),
    [dash.dependencies.Input("menu", "n_clicks")],
    [dash.dependencies.State("mobile_tabs", "style")],
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
