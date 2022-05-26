from dash import html, dcc

from constants import (
    YEARS,
    mapbox_access_token,
    mapbox_style,
)


def Header(app, header, subheader=None):
    left_headers = html.Div(
        [
            html.Div(header, className="header-title"),
            html.Div(subheader, className="subheader-title"),
        ]
    )

    logo = html.Img(src=app.get_asset_url("images/plotly-logo.png"))
    link = html.A(logo, href="https://plotly.com/dash/", target="_blank")
    demo_link = html.A(
        "LEARN MORE",
        href="https://plotly.com/dash/",
        target="_blank",
        className="demo-button",
    )
    right_logos = html.Div([demo_link, link], className="header-logos")

    return html.Div([left_headers, right_logos], className="header")


def choropleth_card(county_choropleth_id):
    return html.Div(
        id="left-column",
        children=[
            html.Div(
                id="slider-container",
                children=[
                    html.P(
                        id="slider-text",
                        children="Drag the slider to change the year:",
                    ),
                    dcc.Slider(
                        id="years-slider",
                        min=min(YEARS),
                        max=max(YEARS),
                        value=min(YEARS),
                        marks={
                            str(year): {
                                "label": str(year),
                                "style": {"color": "#7fafdf"},
                            }
                            for year in YEARS
                        },
                    ),
                ],
            ),
            html.Div(
                id="heatmap-container",
                children=[
                    html.P(
                        "Heatmap of age adjusted mortality rates \
                            from poisonings in year {0}".format(
                            min(YEARS)
                        ),
                        id="heatmap-title",
                    ),
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
                ],
            ),
        ],
    )


def slider_graph_card(selected_data_id):
    return html.Div(
        id="graph-container",
        children=[
            html.P(id="chart-selector", children="Select chart:"),
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
                    {
                        "label": "Trends in age-adjusted death rate (1999-2016)",
                        "value": "death_rate_all_time",
                    },
                ],
                value="show_death_rate_single_year",
                id="chart-dropdown",
            ),
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
        ],
    )
