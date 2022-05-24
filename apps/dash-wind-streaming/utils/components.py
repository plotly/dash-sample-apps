from dash import html, dcc
from utils.helper_functions import app_color
import os

GRAPH_INTERVAL = os.environ.get("GRAPH_INTERVAL", 5000)


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
