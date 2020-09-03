"""
Adapted from: https://pydeck.gl/gallery/geopandas_integration.html

This demos shows how to use the geopandas library in pydeck and Dash Deck.
"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


import geopandas as gpd

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

centroids = gpd.GeoDataFrame()
centroids["geometry"] = world.geometry.centroid
centroids["name"] = world.name

layers = [
    pdk.Layer("GeoJsonLayer", data=world, get_fill_color=[0, 0, 0],),
    pdk.Layer(
        "TextLayer",
        data=centroids,
        get_position="geometry.coordinates",
        get_size=16,
        get_color=[255, 255, 255],
        get_text="name",
    ),
]

r = pdk.Deck(layers, map_provider=None)

app = dash.Dash(__name__)

app.layout = html.Div(dash_deck.DeckGL(r.to_json(), id="deck-gl"))


if __name__ == "__main__":
    app.run_server(debug=True)
