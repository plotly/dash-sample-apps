"""
Adapted from: https://pydeck.gl/gallery/path_layer.html

Locations of the Bay Area Rapid Transit lines.

"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


DATA_URL = "https://raw.githubusercontent.com/visgl/deck.gl-data/master/website/bart-lines.json"
df = pd.read_json(DATA_URL)


def hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i : i + 2], 16) for i in (0, 2, 4))


df["color"] = df["color"].apply(hex_to_rgb)


view_state = pdk.ViewState(latitude=37.782556, longitude=-122.3484867, zoom=10)

layer = pdk.Layer(
    type="PathLayer",
    data=df,
    pickable=True,
    get_color="color",
    width_scale=20,
    width_min_pixels=2,
    get_path="path",
    get_width=5,
)

r = pdk.Deck(layers=[layer], initial_view_state=view_state)


app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(
        r.to_json(),
        id="deck-gl",
        tooltip={"text": "{name}"},
        mapboxKey=mapbox_api_token,
    )
)


if __name__ == "__main__":
    app.run_server(debug=True)
