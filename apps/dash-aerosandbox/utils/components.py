from dash import html, dcc
import dash_bootstrap_components as dbc


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

    logo = html.Img(src=app.get_asset_url("assets/images/plotly-logo-light-theme.png"))
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


def key_parameters_card(booms_id, wing_span_id, alpha_id):
    return html.Div(
        [
            html.H5("Key Parameters"),
            html.P("Number of booms:"),
            dcc.Slider(
                id=booms_id,
                min=1,
                max=3,
                step=1,
                value=3,
                marks={
                    1: "1",
                    2: "2",
                    3: "3",
                },
            ),
            html.P("Wing Span [m]:"),
            dcc.Input(id=wing_span_id, value=43, type="number"),
            html.P("Angle of Attack [deg]:"),
            dcc.Input(id=alpha_id, value=7.0, type="number"),
        ]
    )


def commands_card(display_geometry_id, run_ll_analysis_id, run_vlm_analysis_id):
    return html.Div(
        [
            html.H5("Commands"),
            dbc.Button(
                "Display (1s)",
                id="display_geometry",
                color="primary",
                style={"margin": "5px"},
                n_clicks_timestamp="0",
            ),
            dbc.Button(
                "LL Analysis (3s)",
                id="run_ll_analysis",
                color="secondary",
                style={"margin": "5px"},
                n_clicks_timestamp="0",
            ),
            dbc.Button(
                "VLM Analysis (15s)",
                id="run_vlm_analysis",
                color="secondary",
                style={"margin": "5px"},
                n_clicks_timestamp="0",
            ),
        ]
    )


def aerodynamic_performance_card(output_id):
    return html.Div(
        [
            html.H5("Aerodynamic Performance"),
            dbc.Spinner(
                html.P(id=output_id),
                color="primary",
            ),
        ]
    )


def figure_card(display_id):
    return dbc.Col(
        [
            # html.Div(id='display')
            dbc.Spinner(
                dcc.Graph(id=display_id, style={"height": "80vh"}),
                color="primary",
            )
        ],
        width=True,
    )


# def footer_card():
#     return html.P(
#         [
#             html.A(
#                 "Source code",
#                 href="https://github.com/peterdsharpe/AeroSandbox-Interactive-Demo",
#             ),
#             ". Aircraft design tools powered by ",
#             html.A("AeroSandbox", href="https://peterdsharpe.github.com/AeroSandbox"),
#             ". Build beautiful UIs for your scientific computing apps with ",
#             html.A("Plot.ly ", href="https://plotly.com/"),
#             "and ",
#             html.A("Dash", href="https://plotly.com/dash/"),
#             "!",
#         ]
#     )
