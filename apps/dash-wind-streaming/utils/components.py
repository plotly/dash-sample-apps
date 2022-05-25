from dash import html, dcc

from constants import app_color

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


def wind_speed_card(graph_id):
    return html.Div([
        html.Div(html.H6("WIND SPEED (MPH)", className="graph__title")),
        dcc.Graph(
            id=graph_id,
            figure=dict(
                layout=dict(
                    plot_bgcolor=app_color["graph_bg"],
                    paper_bgcolor=app_color["graph_bg"],
                )
            ),
        ),
    ], className="two-thirds column wind__speed__container",
    )

def histogram_card(graph_id):
    return html.Div([
        html.Div(html.H6("WIND SPEED HISTOGRAM", className="graph__title")),
        html.Div(
            dcc.Slider(
                id="bin-slider",
                min=1,
                max=60,
                step=1,
                value=20,
                marks={
                    20: {"label": "20"},
                    40: {"label": "40"},
                    60: {"label": "60"},
                },
            ),
            className="slider",
        ),
        html.Div([
            dcc.Checklist(
                id="bin-auto",
                options=["Auto"],
                value=["Auto"],
                inputClassName="auto__checkbox",
                labelClassName="auto__label",
            ),
            html.P(
                "# of Bins: Auto",
                id="bin-size",
                className="auto__p",
            )],
            className="auto__container",
        ),
        dcc.Graph(
            id=graph_id,
            figure=dict(
                layout=dict(
                    plot_bgcolor=app_color["graph_bg"],
                    paper_bgcolor=app_color["graph_bg"],
                )
            ),
        ),
    ], className="graph__container first",
    )

def wind_direction_card(graph_id):
    return html.Div([
        html.Div(html.H6("WIND DIRECTION", className="graph__title")),
        dcc.Graph(
            id=graph_id,
            figure=dict(
                layout=dict(
                    plot_bgcolor=app_color["graph_bg"],
                    paper_bgcolor=app_color["graph_bg"],
                )
            ),
        ),
    ], className="graph__container second",
    )