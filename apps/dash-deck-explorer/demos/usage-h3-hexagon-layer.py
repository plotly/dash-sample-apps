"""
Adapted from: https://pydeck.gl/gallery/h3_hexagon_layer.html

Plot of values for a particular hex ID in the H3 geohashing scheme.

This example is adapted from the deck.gl documentation.
"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


H3_HEX_DATA = "https://raw.githubusercontent.com/visgl/deck.gl-data/master/website/sf.h3cells.json"

df = pd.read_json(H3_HEX_DATA)

# Define a layer to display on a map
layer = pdk.Layer(
    "H3HexagonLayer",
    df,
    pickable=True,
    stroked=True,
    filled=True,
    extruded=False,
    get_hexagon="hex",
    get_fill_color="[255 - count, 255, count]",
    get_line_color=[255, 255, 255],
    line_width_min_pixels=2,
)

# Set the viewport location
view_state = pdk.ViewState(
    latitude=37.7749295, longitude=-122.4194155, zoom=14, bearing=0, pitch=30
)


# Render
r = pdk.Deck(layers=[layer], initial_view_state=view_state)


app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(
        r.to_json(),
        id="deck-gl",
        tooltip={"text": "Count: {count}"},
        mapboxKey=mapbox_api_token,
    )
)


if __name__ == "__main__":
    app.run_server(debug=True)
