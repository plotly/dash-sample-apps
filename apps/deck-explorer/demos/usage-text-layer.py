"""
Adapted from: https://pydeck.gl/gallery/text_layer.html

Names of various public transit stops within San Francisco,
plotted at the location of that stop

"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck as pdk
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


TEXT_LAYER_DATA = "https://raw.githubusercontent.com/visgl/deck.gl-data/master/website/bart-stations.json"  # noqa
df = pd.read_json(TEXT_LAYER_DATA)

# Define a layer to display on a map
layer = pdk.Layer(
    "TextLayer",
    df,
    pickable=True,
    get_position="coordinates",
    get_text="name",
    get_size=16,
    get_color=[255, 255, 255],
    get_angle=0,
    # Note that string constants in pydeck are explicitly passed as strings
    # This distinguishes them from columns in a data set
    get_text_anchor="'middle'",
    get_alignment_baseline="'center'",
)

# Set the viewport location
view_state = pdk.ViewState(
    latitude=37.7749295, longitude=-122.4194155, zoom=10, bearing=0, pitch=45
)

# Render
r = pdk.Deck(
    layers=[layer], initial_view_state=view_state, map_style=pdk.map_styles.SATELLITE,
)


app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(
        r.to_json(),
        id="deck-gl",
        tooltip={"text": "{name}\n{address}"},
        mapboxKey=mapbox_api_token,
    )
)


if __name__ == "__main__":
    app.run_server(debug=True)
