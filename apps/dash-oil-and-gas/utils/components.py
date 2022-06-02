from dash import html, dcc
from constants import well_status_options, well_type_options, WELL_STATUSES, WELL_TYPES


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

    logo = html.Img(src=app.get_asset_url("images/plotly-logo-light-theme.png"))
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


def controls_card():
    return html.Div(
        [
            html.P(
                "Filter by construction date (or select range in histogram):",
                className="control_label",
            ),
            dcc.RangeSlider(
                id="year_slider",
                min=1960,
                max=2017,
                value=[1990, 2010],
                step=1,
                marks={
                    1960: "1960",
                    1970: "1970",
                    1980: "1980",
                    1990: "1990",
                    2000: "2000",
                    2010: "2010",
                    2017: "2017",
                },
                className="dcc_control",
            ),
            html.P("Filter by well status:", className="control_label"),
            dcc.RadioItems(
                id="well_status_selector",
                options=[
                    {"label": "All ", "value": "all"},
                    {"label": "Active only ", "value": "active"},
                    {"label": "Customize ", "value": "custom"},
                ],
                value="active",
                labelStyle={"display": "inline-block"},
                className="dcc_control",
            ),
            dcc.Dropdown(
                id="well_statuses",
                options=well_status_options,
                multi=True,
                value=list(WELL_STATUSES.keys()),
                className="dcc_control",
            ),
            dcc.Checklist(
                id="lock_selector",
                options=[{"label": "Lock camera", "value": "locked"}],
                className="dcc_control",
                value=[],
            ),
            html.P("Filter by well type:", className="control_label"),
            dcc.RadioItems(
                id="well_type_selector",
                options=[
                    {"label": "All ", "value": "all"},
                    {"label": "Productive only ", "value": "productive"},
                    {"label": "Customize ", "value": "custom"},
                ],
                value="productive",
                labelStyle={"display": "inline-block"},
                className="dcc_control",
            ),
            dcc.Dropdown(
                id="well_types",
                options=well_type_options,
                multi=True,
                value=list(WELL_TYPES.keys()),
                className="dcc_control",
            ),
        ],
        className="pretty_container four columns",
        id="cross-filter-options",
    )


def top_data_cards():
    return html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [html.H6(id="well_text"), html.P("No. of Wells")],
                        id="wells",
                        className="mini_container",
                    ),
                    html.Div(
                        [html.H6(id="gasText"), html.P("Gas")],
                        id="gas",
                        className="mini_container",
                    ),
                    html.Div(
                        [html.H6(id="oilText"), html.P("Oil")],
                        id="oil",
                        className="mini_container",
                    ),
                    html.Div(
                        [html.H6(id="waterText"), html.P("Water")],
                        id="water",
                        className="mini_container",
                    ),
                ],
                id="info-container",
                className="row container-display",
            ),
            html.Div(
                [dcc.Graph(id="count_graph")],
                id="countGraphContainer",
                className="pretty_container",
            ),
        ],
        id="right-column",
        className="eight columns",
    )


def main_graph():
    return html.Div(
        [dcc.Graph(id="main_graph")],
        className="pretty_container seven columns",
    )


def individual_graph():
    return html.Div(
        [dcc.Graph(id="individual_graph")],
        className="pretty_container five columns",
    )


def pie_graph():
    return html.Div(
        [dcc.Graph(id="pie_graph")],
        className="pretty_container seven columns",
    )


def aggregate_graph():
    return html.Div(
        [dcc.Graph(id="aggregate_graph")],
        className="pretty_container five columns",
    )
