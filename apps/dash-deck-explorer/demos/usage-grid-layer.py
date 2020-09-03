"""
Adapted from: https://pydeck.gl/gallery/grid_layer.html

Locations of bike parking within San Francisco.

Adapted from the deck.gl documentation.
"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


CPU_GRID_LAYER_DATA = (
    "https://raw.githubusercontent.com/uber-common/"
    "deck.gl-data/master/website/sf-bike-parking.json"
)
df = pd.read_json(CPU_GRID_LAYER_DATA)

# Define a layer to display on a map

layer = pdk.Layer(
    "GridLayer",
    df,
    pickable=True,
    extruded=True,
    cell_size=200,
    elevation_scale=4,
    get_position="COORDINATES",
)

view_state = pdk.ViewState(
    latitude=37.7749295, longitude=-122.4194155, zoom=11, bearing=0, pitch=45
)

# Render
r = pdk.Deck(layers=[layer], initial_view_state=view_state)


app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(
        r.to_json(),
        id="deck-gl",
        tooltip={"text": "{position}\nCount: {count}"},
        mapboxKey=mapbox_api_token,
    )
)


if __name__ == "__main__":
    app.run_server(debug=True)
