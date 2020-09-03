"""
Adapted from: https://pydeck.gl/gallery/terrain_layer.html

Extruded terrain using AWS Open Data Terrain Tiles and Mapbox Satellite imagery

"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk

MAPBOX_API_KEY = os.getenv("MAPBOX_ACCESS_TOKEN")

# AWS Open Data Terrain Tiles
TERRAIN_IMAGE = (
    '"https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"'
)

# Define how to parse elevation tiles
ELEVATION_DECODER = {"rScaler": 256, "gScaler": 1, "bScaler": 1 / 256, "offset": -32768}

SURFACE_IMAGE = '"https://api.mapbox.com/v4/mapbox.satellite/{{z}}/{{x}}/{{y}}@2x.png?access_token={}"'.format(
    MAPBOX_API_KEY
)

terrain_layer = pdk.Layer(
    "TerrainLayer",
    data=None,
    elevation_decoder=ELEVATION_DECODER,
    texture=SURFACE_IMAGE,
    elevation_data=TERRAIN_IMAGE,
)

view_state = pdk.ViewState(
    latitude=46.24, longitude=-122.18, zoom=11.5, bearing=140, pitch=60
)

r = pdk.Deck(terrain_layer, initial_view_state=view_state)


app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(r.to_json(), id="deck-gl", mapboxKey=MAPBOX_API_KEY)
)


if __name__ == "__main__":
    app.run_server(debug=True)
