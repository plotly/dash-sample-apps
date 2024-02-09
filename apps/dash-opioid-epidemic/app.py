from dash import dash, html, Input, Output, State, callback
import dash_bootstrap_components as dbc

import utils.figures as figs
from utils.components import header, choropleth_card, slider_graph_card

app = dash.Dash(__name__, title = "US Opioid Epidemic", external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

# App layout
app.layout = dbc.Container(
    [
        header(
            app,
            "inherit",
            "Rate of US Poison-Induced Deaths",
            subheader="† Deaths are classified using the International Classification of Diseases, \
                    Tenth Revision (ICD-10).\n\nDrug-poisoning deaths are defined as having ICD-10 underlying \
                    cause-of-death codes X40-X44 (unintentional), X60-X64 (suicide), X85 (homicide), or Y10-Y14 \
                    (undetermined intent).",
        ),
        dbc.Row([
                dbc.Col(choropleth_card("county-choropleth"), width=7),
                dbc.Col(slider_graph_card("selected-data"),width=5)
            ],
        ),
    ],
    fluid=True
)


@callback(
    Output("county-choropleth", "figure"),
    Input("years-slider", "value"),
    State("county-choropleth", "figure"),
)
def return_display_map(year, figure):
    return figs.display_map(year, figure)


@callback(
    Output("heatmap-title", "children"), 
    Input("years-slider", "value")
)
def update_map_title(year):
    return f"Heatmap of age adjusted mortality rates from poisonings in year {year}"


@callback(
    Output("selected-data", "figure"),
    Input("county-choropleth", "selectedData"),
    Input("chart-dropdown", "value"),
    Input("years-slider", "value"),
)
def return_display_selected_data(selectedData, chart_dropdown, year):
    return figs.display_selected_data(selectedData, chart_dropdown, year)


if __name__ == "__main__":
    app.run_server(debug=True)
