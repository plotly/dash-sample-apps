from dash import html, dcc
from constants import state_list, cost_metric, state_map, data_dict, init_region
from utils.figures import generate_procedure_plot


def header(
    app, header_color, header, subheader=None, header_background_color="transparent"
):
    left_headers = html.Div(
        [
            html.Div(header, className="header-title"),
            html.Div(subheader, className="subheader-title"),
        ],
        style={"color": header_color},
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

    return html.Div(
        [left_headers, right_logos],
        className="header",
        style={"background-color": header_background_color},
    )


def build_upper_left_panel():
    return html.Div(
        id="upper-left",
        className="six columns",
        children=[
            html.P(
                className="section-title",
                children="Choose hospital on the map or procedures from the list below to see costs",
            ),
            html.Div(
                className="control-row-1",
                children=[
                    html.Div(
                        id="state-select-outer",
                        children=[
                            html.Label("Select a State"),
                            dcc.Dropdown(
                                id="state-select",
                                options=[{"label": i, "value": i} for i in state_list],
                                value=state_list[1],
                            ),
                        ],
                    ),
                    html.Div(
                        id="select-metric-outer",
                        children=[
                            html.Label("Choose a Cost Metric"),
                            dcc.Dropdown(
                                id="metric-select",
                                options=[{"label": i, "value": i} for i in cost_metric],
                                value=cost_metric[0],
                            ),
                        ],
                    ),
                ],
            ),
            html.Div(
                id="region-select-outer",
                className="control-row-2",
                children=[
                    html.Label("Pick a Region"),
                    html.Div(
                        id="checklist-container",
                        children=dcc.Checklist(
                            id="region-select-all",
                            options=[{"label": "Select All Regions", "value": "All"}],
                            value=[],
                        ),
                    ),
                    html.Div(
                        id="region-select-dropdown-outer",
                        children=dcc.Dropdown(
                            id="region-select",
                            multi=True,
                            searchable=True,
                        ),
                    ),
                ],
            ),
            html.Div(
                id="table-container",
                className="table-container",
                children=[
                    html.Div(
                        id="table-upper",
                        children=[
                            html.P("Hospital Charges Summary"),
                            dcc.Loading(children=html.Div(id="cost-stats-container")),
                        ],
                    ),
                    html.Div(
                        id="table-lower",
                        children=[
                            html.P("Procedure Charges Summary"),
                            dcc.Loading(
                                children=html.Div(id="procedure-stats-container")
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )


def map_card():
    return html.Div(
        id="geo-map-outer",
        className="six columns",
        children=[
            html.P(
                id="map-title",
                children="Medicare Provider Charges in the State of {}".format(
                    state_map[state_list[0]]
                ),
            ),
            html.Div(
                id="geo-map-loading-outer",
                children=[
                    dcc.Loading(
                        id="loading",
                        children=dcc.Graph(
                            id="geo-map",
                            figure={
                                "data": [],
                                "layout": dict(
                                    plot_bgcolor="#171b26",
                                    paper_bgcolor="#171b26",
                                ),
                            },
                        ),
                    )
                ],
            ),
        ],
    )


def procedure_card():
    return html.Div(
        id="lower-container",
        children=[
            dcc.Graph(
                id="procedure-plot",
                figure=generate_procedure_plot(
                    data_dict[state_list[1]], cost_metric[0], init_region, []
                ),
            )
        ],
    )
