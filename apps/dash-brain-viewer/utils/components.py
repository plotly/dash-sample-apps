from dash import html, dcc
import dash_colorscales as dcs

from utils.model import create_mesh_data
from constants import default_colorscale_index, plot_layout


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


def brain_graph(brain_graph_id):
    return dcc.Graph(
        id=brain_graph_id,
        figure={
            "data": create_mesh_data("human_atlas"),
            "layout": plot_layout,
        },
        config={"editable": True, "scrollZoom": False},
    )


def control_and_output():
    return [
        html.Div(
            [
                html.P("Click colorscale to change"),
                dcs.DashColorscales(
                    id="colorscale-picker",
                    colorscale=default_colorscale_index,
                ),
            ]
        ),
        html.Div(
            [
                html.P("Select option"),
                dcc.RadioItems(
                    options=[
                        {"label": "Brain Atlas", "value": "human_atlas"},
                        {"label": "Cortical Thickness", "value": "human"},
                        {"label": "Mouse Brain", "value": "mouse"},
                    ],
                    value="human_atlas",
                    id="radio-options",
                ),
            ],
        ),
        html.Div(
            [
                html.Span("Click data"),
                html.Span("  |  "),
                html.Span("Click on points in the graph."),
                dcc.Loading(
                    html.Pre(id="click-data"),
                    type="dot",
                ),
            ],
        ),
        html.Div(
            [
                html.Span("Relayout data"),
                html.Span("  |  "),
                html.Span("Drag the graph corners to rotate it."),
                dcc.Loading(
                    html.Pre(id="relayout-data"),
                    type="dot",
                ),
            ],
        ),
        html.P(
            [
                "Brain data from Mcgill's ACE Lab ",
                html.A(
                    children="Surface Viewer.",
                    target="_blank",
                    href="https://brainbrowser.cbrain.mcgill.ca/surface-viewer#ct",
                ),
            ]
        ),
    ]
