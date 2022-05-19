from dash import html
import dash_bootstrap_components as dbc
import dash_deck

from constants import CAMERAS, LIDARS, token_list
from utils.helper_functions import unsnake

def Header(name, app):
    title = html.H2(name, style={"margin-top": 5})
    logo = html.Img(src=app.get_asset_url("images/dash-logo.png"), style={"float": "right", "height": 60})
    link = html.A(logo, href="https://plotly.com/dash/", target="_blank") 
    demo_link = html.A("ENTERPRISE DEMO", href="https://plotly.com/get-demo/", target="_blank", className="demo-button")

    return dbc.Row([dbc.Col(title, md=8), dbc.Col([demo_link, link], md=4, className="header-logos")], className="header")


CONTROLS = [
    html.Div(
        [
            dbc.Label("Camera Position"),
            dbc.Select(
                id="camera",
                options=[
                    {"label": unsnake(s.replace("CAM_", "")), "value": s}
                    for s in CAMERAS
                ],
                value=CAMERAS[0],
            ),
        ]
    ),
    html.Div(
        [
            dbc.Label("Image Overlay"),
            dbc.Checklist(
                id="overlay",
                value=["boxes"],
                options=[
                    {"label": x.title(), "value": x} for x in ["pointcloud", "boxes"]
                ],
                inline=True,
                switch=True,
            ),
        ]
    ),
    html.Div(
        [
            dbc.Label("Frame"),
            html.Br(),
            dbc.Spinner(
                dbc.ButtonGroup(
                    [
                        dbc.Button(
                            "Prev", id="prev", n_clicks=0, color="primary", outline=True
                        ),
                        dbc.Button("Next", id="next", n_clicks=0, color="primary"),
                    ],
                    id="button-group",
                    style={"width": "100%"},
                ),
                spinner_style={"margin-top": 0, "margin-bottom": 0},
            ),
        ]
    ),
    html.Div(
        [
            dbc.Label("Progression"),
            dbc.Spinner(
                dbc.Input(id="progression", type="range", min=0, max=len(token_list), value=0, step=1),
                spinner_style={"margin-top": 0, "margin-bottom": 0},
            ),
        ]
    ),
    html.Div(
        [
            dbc.Label("Lidar Position"),
            dbc.Select(
                id="lidar",
                value=LIDARS[0],
                options=[
                    {"label": unsnake(s.replace("LIDAR_", "")), "value": s}
                    for s in LIDARS
                ],
            ),
        ]
    ),
    html.Div(
        [
            dbc.Label("Lidar View Mode"),
            dbc.Select(
                id="view-mode",
                value="map",
                options=[
                    {"label": unsnake(x), "value": x}
                    for x in ["first_person", "orbit", "map"]
                ],
            ),
        ]
    ),
]

DECK_CARD = dbc.Card(
    dash_deck.DeckGL(id="deck-pointcloud", tooltip={"html": "<b>Label:</b> {name}"}),
    body=True,
    style={"height": "calc(95vh - 215px)"},
)