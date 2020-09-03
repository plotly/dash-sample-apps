"""
Adapted from: https://pydeck.gl/gallery/great_circle_layer.html

Plot of direct flights to and from San Francisco International Airport.

Origin is in green; destinations are in blue.

Adapted from the deck.gl documentation.
"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


GREAT_CIRCLE_LAYER_DATA = "https://raw.githubusercontent.com/visgl/deck.gl-data/master/website/flights.json"  # noqa

df = pd.read_json(GREAT_CIRCLE_LAYER_DATA)

# Use pandas to prepare data for tooltip
df["from_name"] = df["from"].apply(lambda f: f["name"])
df["to_name"] = df["to"].apply(lambda t: t["name"])

# Define a layer to display on a map
layer = pdk.Layer(
    "GreatCircleLayer",
    df,
    pickable=True,
    get_stroke_width=12,
    get_source_position="from.coordinates",
    get_target_position="to.coordinates",
    get_source_color=[64, 255, 0],
    get_target_color=[0, 128, 200],
    auto_highlight=True,
)

# Set the viewport location
view_state = pdk.ViewState(latitude=50, longitude=-40, zoom=1, bearing=0, pitch=0)

# Render
r = pdk.Deck(layers=[layer], initial_view_state=view_state,)
r.picking_radius = 10


app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(
        r.to_json(),
        id="deck-gl",
        mapboxKey=mapbox_api_token,
        tooltip={"text": "{from_name} to {to_name}"},
    )
)


if __name__ == "__main__":
    app.run_server(debug=True)
