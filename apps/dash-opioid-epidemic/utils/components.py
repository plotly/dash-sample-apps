from dash import html, dcc
import dash_bootstrap_components as dbc

from constants import YEARS, mapbox_access_token, mapbox_style


def header(app, header_color, header, subheader=None, header_background_color="transparent"):
    left_headers = html.Div(
        [
            html.Div(header, className="header-title"),
            html.Div(subheader, className="subheader-title"),
        ],
        style={"color": header_color}
    )

    logo = html.Img(src=app.get_asset_url("images/plotly-logo-dark-theme.png"))
    logo_link = html.A(logo, href="https://plotly.com/get-demo/", target="_blank")
    demo_link = html.A(
        "LEARN MORE",
        href="https://plotly.com/dash/",
        target="_blank",
        className="demo-button",
    )
    right_logos = html.Div([demo_link, logo_link], className="header-logos")

    return html.Div([left_headers, right_logos], className="header", style={"background-color": header_background_color})


def choropleth_card(county_choropleth_id):
    return dbc.Row([
            dbc.Card(dbc.CardBody([
                html.P("Drag the slider to change the year:"),
                dcc.Slider(
                    id="years-slider",
                    min=min(YEARS),
                    max=max(YEARS),
                    value=min(YEARS),
                    step=1,
                    marks={
                        str(year): {
                            "label": str(year),
                            "style": {"color": "#7fafdf"},
                        }
                        for year in YEARS
                    },
                ),
            ])),
            dbc.Card(dbc.CardBody([
                html.P(f"Heatmap of age adjusted mortality rates from poisonings in year {min(YEARS)}", id="heatmap-title"),
                dcc.Graph(
                    id=county_choropleth_id,
                    figure=dict(
                        layout=dict(
                            mapbox=dict(
                                layers=[],
                                accesstoken=mapbox_access_token,
                                style=mapbox_style,
                                center=dict(lat=38.72490, lon=-95.61446),
                                pitch=0,
                                zoom=3.5,
                            ),
                            autosize=True,
                        ),
                    ),
                ),
            ])),
        ],
    )


def slider_graph_card(selected_data_id):
    return dbc.Card(dbc.CardBody([
            html.P(children="Select chart:"),
            dcc.Dropdown(
                options=[
                    {
                        "label": "Histogram of total number of deaths (single year)",
                        "value": "show_absolute_deaths_single_year",
                    },
                    {
                        "label": "Histogram of total number of deaths (1999-2016)",
                        "value": "absolute_deaths_all_time",
                    },
                    {
                        "label": "Age-adjusted death rate (single year)",
                        "value": "show_death_rate_single_year",
                    },
                    # { # Takes too long to compute for many counties
                    #     "label": "Trends in age-adjusted death rate (1999-2016)",
                    #     "value": "death_rate_all_time",
                    # },
                ],
                value="show_death_rate_single_year",
                id="chart-dropdown",
            ),
            dcc.Loading(
                dcc.Graph(
                    id=selected_data_id,
                    figure=dict(
                        data=[dict(x=0, y=0)],
                        layout=dict(
                            paper_bgcolor="#F4F4F8",
                            plot_bgcolor="#F4F4F8",
                            autofill=True,
                            margin=dict(t=75, r=50, b=100, l=50),
                        ),
                    ),
                ),
                type="dot"
            )
        ],
    ))
