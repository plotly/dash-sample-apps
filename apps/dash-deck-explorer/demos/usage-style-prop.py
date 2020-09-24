"""
Adapted from: https://pypi.org/project/pydeck/

This demo shows how to resize the map using the `style` prop.
"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk
import pandas as pd

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

# Render
r = pdk.Deck(
    layers=[layer], initial_view_state=view_state, mapbox_key=mapbox_api_token,
)

app = dash.Dash(__name__)

app.layout = html.Div(
    [
        dash_deck.DeckGL(
            r.to_json(),
            id="deck-gl",
            tooltip=True,
            mapboxKey=r.mapbox_key,
            style={"width": "55vw", "height": "75vh"},
        ),
        # Random Red Rectangle float on the right side
        html.Div(
            style={
                "width": "40vw",
                "height": "75vh",
                "background-color": "red",
                "float": "right",
            }
        ),
    ]
)


if __name__ == "__main__":
    app.run_server(debug=True)
