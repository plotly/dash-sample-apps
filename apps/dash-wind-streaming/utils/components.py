from dash import html, dcc
from utils.helper_functions import app_color
import os

GRAPH_INTERVAL = os.environ.get("GRAPH_INTERVAL", 5000)


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
        "ENTERPRISE DEMO",
        href="https://plotly.com/dash/",
        target="_blank",
        className="demo-button",
    )
    right_logos = html.Div([demo_link, link], className="header-logos")

    return html.Div([left_headers, right_logos], className="header")


slider = [
    dcc.Slider(
        id="bin-slider",
        min=1,
        max=60,
        step=1,
        value=20,
        updatemode="drag",
        marks={
            20: {"label": "20"},
            40: {"label": "40"},
            60: {"label": "60"},
        },
    )
]
checklist = [
    dcc.Checklist(
        id="bin-auto",
        options=[{"label": "Auto", "value": "Auto"}],
        value=["Auto"],
        inputClassName="auto__checkbox",
        labelClassName="auto__label",
    ),
    html.P(
        "# of Bins: Auto",
        id="bin-size",
        className="auto__p",
    ),
]

left_graph = dcc.Graph(
    id="wind-speed",
    figure=dict(
        layout=dict(
            plot_bgcolor=app_color["graph_bg"],
            paper_bgcolor=app_color["graph_bg"],
        )
    ),
)

interval = dcc.Interval(
    id="wind-speed-update",
    interval=int(GRAPH_INTERVAL),
    n_intervals=0,
)

right_graph_one = dcc.Graph(
    id="wind-histogram",
    figure=dict(
        layout=dict(
            plot_bgcolor=app_color["graph_bg"],
            paper_bgcolor=app_color["graph_bg"],
        )
    ),
)

right_graph_two = dcc.Graph(
    id="wind-direction",
    figure=dict(
        layout=dict(
            plot_bgcolor=app_color["graph_bg"],
            paper_bgcolor=app_color["graph_bg"],
        )
    ),
)
