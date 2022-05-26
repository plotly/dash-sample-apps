from dash import html, dcc
from mni import create_mesh_data
from constants import default_colorscale_index, plot_layout, GITHUB_LINK
import dash_colorscales as dcs


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


def brain_graph(brain_graph_id):
    return html.Div(
        [
            dcc.Graph(
                id=brain_graph_id,
                figure={
                    "data": create_mesh_data("human_atlas"),
                    "layout": plot_layout,
                },
                config={"editable": True, "scrollZoom": False},
            )
        ],
        className="graph__container",
    )


def color_picker(color_picker_id):
    return html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            html.P("Click colorscale to change", className="subheader"),
                            dcs.DashColorscales(
                                id=color_picker_id,
                                colorscale=default_colorscale_index,
                            ),
                        ]
                    )
                ],
                className="colorscale pb-20",
            ),
            html.Div(
                [
                    html.P("Select option", className="subheader"),
                    dcc.RadioItems(
                        options=[
                            {"label": "Brain Atlas", "value": "human_atlas"},
                            {"label": "Cortical Thickness", "value": "human"},
                            {"label": "Mouse Brain", "value": "mouse"},
                        ],
                        value="human_atlas",
                        id="radio-options",
                        labelClassName="label__option",
                        inputClassName="input__option",
                    ),
                ],
                className="pb-20",
            ),
            html.Div(
                [
                    html.Span("Click data", className="subheader"),
                    html.Span("  |  "),
                    html.Span("Click on points in the graph.", className="small-text"),
                    dcc.Loading(
                        html.Pre(id="click-data", className="info__container"),
                        type="dot",
                    ),
                ],
                className="pb-20",
            ),
            html.Div(
                [
                    html.Span("Relayout data", className="subheader"),
                    html.Span("  |  "),
                    html.Span(
                        "Drag the graph corners to rotate it.",
                        className="small-text",
                    ),
                    dcc.Loading(
                        html.Pre(id="relayout-data", className="info__container"),
                        type="dot",
                    ),
                ],
                className="pb-20",
            ),
            html.Div(
                [
                    html.P(
                        [
                            "Dash/Python code on ",
                            html.A(
                                children="GitHub.",
                                target="_blank",
                                href=GITHUB_LINK,
                                className="red-ish",
                            ),
                        ]
                    ),
                    html.P(
                        [
                            "Brain data from Mcgill's ACE Lab ",
                            html.A(
                                children="Surface Viewer.",
                                target="_blank",
                                href="https://brainbrowser.cbrain.mcgill.ca/surface-viewer#ct",
                                className="red-ish",
                            ),
                        ]
                    ),
                ]
            ),
        ],
        className="one-third column app__right__section",
    )
