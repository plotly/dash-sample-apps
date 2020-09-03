"""
Adapted from: https://pydeck.gl/gallery/icon_layer.html

Location of biergartens in Germany listed on OpenStreetMap as of early 2020.

"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


# Data from OpenStreetMap, accessed via osmpy
DATA_URL = "https://raw.githubusercontent.com/ajduberstein/geo_datasets/master/biergartens.json"
ICON_URL = (
    "https://upload.wikimedia.org/wikipedia/commons/c/c4/Projet_bi%C3%A8re_logo_v2.png"
)

icon_data = {
    # Icon from Wikimedia, used the Creative Commons Attribution-Share Alike 3.0
    # Unported, 2.5 Generic, 2.0 Generic and 1.0 Generic licenses
    "url": ICON_URL,
    "width": 242,
    "height": 242,
    "anchorY": 242,
}

data = pd.read_json(DATA_URL)
data["icon_data"] = None
for i in data.index:
    data["icon_data"][i] = icon_data

view_state = pdk.data_utils.compute_view(data[["lon", "lat"]], 0.1)

icon_layer = pdk.Layer(
    type="IconLayer",
    data=data,
    get_icon="icon_data",
    get_size=4,
    size_scale=15,
    get_position=["lon", "lat"],
    pickable=True,
)

r = pdk.Deck(layers=[icon_layer], initial_view_state=view_state)


app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(
        r.to_json(),
        id="deck-gl",
        tooltip={"text": "{tags}"},
        mapboxKey=mapbox_api_token,
    )
)


if __name__ == "__main__":
    app.run_server(debug=True)
