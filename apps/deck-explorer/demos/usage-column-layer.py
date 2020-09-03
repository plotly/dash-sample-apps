"""
Adapted from: https://pydeck.gl/gallery/column_layer.html

Real estate values for select properties in Taipei. Data is from 2012-2013.
The height of a column indicates increasing price per unit area, and the color 
indicates distance from a subway stop.

The real estate valuation data set from UC Irvine's Machine Learning repository, 
viewable here:
https://archive.ics.uci.edu/ml/datasets/Real+estate+valuation+data+set
"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")

DATA_URL = (
    "https://raw.githubusercontent.com/ajduberstein/geo_datasets/master/housing.csv"
)
df = pd.read_csv(DATA_URL)

view = pydeck.data_utils.compute_view(df[["lng", "lat"]])
view.pitch = 75
view.bearing = 60

column_layer = pydeck.Layer(
    "ColumnLayer",
    data=df,
    get_position=["lng", "lat"],
    get_elevation="price_per_unit_area",
    elevation_scale=100,
    radius=50,
    get_fill_color=["mrt_distance * 10", "mrt_distance", "mrt_distance * 10", 140],
    pickable=True,
    auto_highlight=True,
)

tooltip = {
    "html": "<b>{mrt_distance}</b> meters away from an MRT station, costs <b>{price_per_unit_area}</b> NTD/sqm",
    "style": {
        "background": "grey",
        "color": "white",
        "font-family": '"Helvetica Neue", Arial',
        "z-index": "10000",
    },
}

r = pydeck.Deck(
    column_layer,
    initial_view_state=view,
    map_style="mapbox://styles/mapbox/satellite-v9",
    mapbox_key=mapbox_api_token,
)

app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(r.to_json(), id="deck-gl", tooltip=tooltip, mapboxKey=r.mapbox_key)
)


if __name__ == "__main__":
    app.run_server(debug=True)
