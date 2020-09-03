"""
Adapted from: https://pydeck.gl/gallery/globe_view.html

This demos the experimental Glove View from deck.gl/pydeck by using the GeoJSON
and column layers. The data used contains global plant database and can be found here:
https://github.com/ajduberstein/geo_datasets

Notice that here, we are explicitly convert the r.to_json() into a python dictionary.
This is needed because the data contains NaN, which can't be parsed by the underlying
JavaScript JSON parser, but it can be parsed by Python's JSON engine.

"""
import os
import json

import dash
import dash_deck
import dash_html_components as html
import pydeck
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


COUNTRIES = "https://d2ad6b4ur7yvpq.cloudfront.net/naturalearth-3.3.0/ne_50m_admin_0_scale_rank.geojson"
POWER_PLANTS = "https://raw.githubusercontent.com/ajduberstein/geo_datasets/master/global_power_plant_database.csv"

df = pd.read_csv(POWER_PLANTS)


def is_green(fuel_type):
    if fuel_type.lower() in (
        "nuclear",
        "water",
        "wind",
        "hydro",
        "biomass",
        "solar",
        "geothermal",
    ):
        return [10, 230, 120]
    return [230, 158, 10]


df["color"] = df["primary_fuel"].apply(is_green)

view_state = pydeck.ViewState(latitude=51.47, longitude=0.45, zoom=2)

layers = []
# Set height and width variables
view = pydeck.View(type="_GlobeView", controller=True, width=1000, height=700)


layers = [
    pydeck.Layer(
        "GeoJsonLayer",
        id="base-map",
        data=COUNTRIES,
        stroked=False,
        filled=True,
        get_line_color=[60, 60, 60],
        get_fill_color=[200, 200, 200],
    ),
    pydeck.Layer(
        "ColumnLayer",
        id="power-plant",
        data=df,
        get_elevation="capacity_mw",
        get_position=["longitude", "latitude"],
        elevation_scale=100,
        pickable=True,
        auto_highlight=True,
        radius=20000,
        get_fill_color="color",
    ),
]

r = pydeck.Deck(
    views=[view],
    initial_view_state=view_state,
    layers=layers,
    # Note that this must be set for the globe to be opaque
    parameters={"cull": True},
)


app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(
        json.loads(r.to_json()),
        id="deck-gl",
        style={"background-color": "black"},
        tooltip={"text": "{name}, {primary_fuel} plant, {country}"},
    )
)


if __name__ == "__main__":
    app.run_server(debug=True)
