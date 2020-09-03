"""
This demos shows you how to customize views using the JSON API.
"""
import dash
import dash_deck
import dash_html_components as html


data = {
    "description": "A minimal deck.gl example rendering a circle with text",
    "initialViewState": {"longitude": -122.45, "latitude": 37.8, "zoom": 12},
    "layers": [
        {
            "@@type": "ScatterplotLayer",
            "data": [{"position": [-122.45, 37.8]}],
            "getFillColor": [255, 0, 0, 255],
            "getRadius": 1000,
        },
        {
            "@@type": "TextLayer",
            "data": [{"position": [-122.45, 37.8], "text": "Hello World"}],
        },
    ],
    "views": [{"@@type": "MapView", "width": "75%", "height": "600px"}],
}

app = dash.Dash(__name__)

app.layout = html.Div(dash_deck.DeckGL(data, id="deck-gl"))


if __name__ == "__main__":
    app.run_server(debug=True)
