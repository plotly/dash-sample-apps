"""
Adapted from: https://pydeck.gl/gallery/h3_cluster_layer.html

Data grouped by H3 geohash, as an example of one of the geohash schemes 
supported by pydeck. If you'd simply like to plot a value at an H3 hex 
ID, see the H3HexagonLayer.

"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


H3_CLUSTER_LAYER_DATA = "https://raw.githubusercontent.com/visgl/deck.gl-data/master/website/sf.h3clusters.json"  # noqa

df = pd.read_json(H3_CLUSTER_LAYER_DATA)

# Define a layer to display on a map
layer = pdk.Layer(
    "H3ClusterLayer",
    df,
    pickable=True,
    stroked=True,
    filled=True,
    extruded=False,
    get_hexagons="hexIds",
    get_fill_color="[255, (1 - mean / 500) * 255, 0]",
    get_line_color=[255, 255, 255],
    line_width_min_pixels=2,
)

# Set the viewport location
view_state = pdk.ViewState(
    latitude=37.7749295, longitude=-122.4194155, zoom=11, bearing=0, pitch=30
)


# Render
r = pdk.Deck(layers=[layer], initial_view_state=view_state)


app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(
        r.to_json(),
        id="deck-gl",
        tooltip={"text": "Density: {mean}"},
        mapboxKey=mapbox_api_token,
    )
)


if __name__ == "__main__":
    app.run_server(debug=True)
