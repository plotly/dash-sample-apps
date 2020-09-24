"""
This demo shows how to interact with event callbacks
like clickInfo, hoverInfo, dragStartInfo, etc.
"""
import os
import json

import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import dash_deck
import pydeck as pdk

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


# 2014 locations of car accidents in the UK
UK_ACCIDENTS_DATA = (
    "https://raw.githubusercontent.com/uber-common/"
    "deck.gl-data/master/examples/3d-heatmap/heatmap-data.csv"
)

# Define a layer to display on a map
layer = pdk.Layer(
    "HexagonLayer",
    UK_ACCIDENTS_DATA,
    get_position=["lng", "lat"],
    auto_highlight=True,
    elevation_scale=50,
    pickable=True,
    elevation_range=[0, 3000],
    extruded=True,
    coverage=1,
)

# Set the viewport location
view_state = pdk.ViewState(
    longitude=-1.415,
    latitude=52.2323,
    zoom=6,
    min_zoom=5,
    max_zoom=15,
    pitch=40.5,
    bearing=-27.36,
)

map_view = pdk.View("MapView", controller=True)

# Render
r = pdk.Deck(
    layers=[layer],
    initial_view_state=view_state,
    mapbox_key=mapbox_api_token,
    views=[map_view],
)

# Start building the layout here
styles = {
    "json-output": {
        "overflowY": "scroll",
        "height": "calc(50% - 25px)",
        "border": "thin lightgrey solid",
    },
    "tab": {"height": "calc(98vh - 115px)"},
}

tabs = [
    dcc.Tab(
        html.Div(
            style=styles["tab"],
            children=[
                html.P("hoverInfo"),
                html.Pre(id="hover-info-json-output", style=styles["json-output"],),
                html.P("hoverEvent"),
                html.Pre(id="hover-event-json-output", style=styles["json-output"],),
            ],
        ),
        label="Hover",
    ),
    dcc.Tab(
        html.Div(
            style=styles["tab"],
            children=[
                html.P("clickInfo"),
                html.Pre(id="click-info-json-output", style=styles["json-output"],),
                html.P("clickEvent"),
                html.Pre(id="click-event-json-output", style=styles["json-output"],),
            ],
        ),
        label="click",
    ),
    dcc.Tab(
        html.Div(
            style=styles["tab"],
            children=[
                html.P("dragStartInfo"),
                html.Pre(id="dragStart-info-json-output", style=styles["json-output"],),
                html.P("dragStartEvent"),
                html.Pre(
                    id="dragStart-event-json-output", style=styles["json-output"],
                ),
            ],
        ),
        label="dragStart",
    ),
    dcc.Tab(
        html.Div(
            style=styles["tab"],
            children=[
                html.P("dragEndInfo"),
                html.Pre(id="dragEnd-info-json-output", style=styles["json-output"],),
                html.P("dragEndEvent"),
                html.Pre(id="dragEnd-event-json-output", style=styles["json-output"],),
            ],
        ),
        label="dragEnd",
    ),
]

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(
    [
        html.Div(
            style={
                "width": "64%",
                "height": "95vh",
                "display": "inline-block",
                "position": "relative",
            },
            children=[
                dash_deck.DeckGL(
                    r.to_json(),
                    id="deck",
                    tooltip=True,
                    enableEvents=True,
                    mapboxKey=r.mapbox_key,
                )
            ],
        ),
        html.Div(
            style={"width": "34%", "float": "right", "display": "inline-block"},
            children=dcc.Tabs(id="tabs", children=tabs),
        ),
    ]
)


def assign_callback(app, out, event):
    @app.callback(Output(f"{out}-json-output", "children"), [Input("deck", event)])
    def dump_json(data):
        return json.dumps(data, indent=2)


assign_callback(app, "click-info", "clickInfo")
assign_callback(app, "click-event", "clickEvent")
assign_callback(app, "hover-info", "hoverInfo")
assign_callback(app, "hover-event", "hoverEvent")
assign_callback(app, "dragStart-info", "dragStartInfo")
assign_callback(app, "dragStart-event", "dragStartEvent")
assign_callback(app, "dragEnd-info", "dragEndInfo")
assign_callback(app, "dragEnd-event", "dragEndEvent")

if __name__ == "__main__":
    app.run_server(debug=True)
