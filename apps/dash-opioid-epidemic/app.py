import pandas as pd
from dash import dash, html, dcc, Input, Output, State, callback, callback_context
import cufflinks as cf

from utils.figures import display_map, display_selected_data

from utils.components import Header, choropleth_card, slider_graph_card

from constants import (
    YEARS,
    mapbox_access_token,
    mapbox_style,
)

# Initialize app

app = dash.Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.0"}
    ],
)
app.title = "US Opioid Epidemic"
server = app.server


# App layout
app.layout = html.Div(
    [
        Header(
            app,
            "Rate of US Poison-Induced Deaths",
            "† Deaths are classified using the International Classification of Diseases, \
                    Tenth Revision (ICD–10). Drug-poisoning deaths are defined as having ICD–10 underlying \
                    cause-of-death codes X40–X44 (unintentional), X60–X64 (suicide), X85 (homicide), or Y10–Y14 \
                    (undetermined intent).",
        ),
        html.Div(
            id="app-container",
            children=[
                choropleth_card("county-choropleth"),
                slider_graph_card("selected-data"),
            ],
        ),
    ],
)


@callback(
    Output("county-choropleth", "figure"),
    Input("years-slider", "value"),
    State("county-choropleth", "figure"),
)
def return_display_map(year, figure):
    return display_map(year, figure)


@callback(Output("heatmap-title", "children"), Input("years-slider", "value"))
def update_map_title(year):
    return "Heatmap of age adjusted mortality rates \
				from poisonings in year {0}".format(
        year
    )


@callback(
    Output("selected-data", "figure"),
    Input("county-choropleth", "selectedData"),
    Input("chart-dropdown", "value"),
    Input("years-slider", "value"),
)
def return_display_selected_data(selectedData, chart_dropdown, year):
    return display_selected_data(selectedData, chart_dropdown, year)


if __name__ == "__main__":
    app.run_server(debug=True)
