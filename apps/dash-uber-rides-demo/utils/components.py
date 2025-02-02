from datetime import datetime as dt
from dash import html, dcc
import dash_bootstrap_components as dbc

from constants import dict_of_locations

def controls(app):
    return dbc.Col(
    className="div-user-controls",
    md=4,
    children=[
        html.Div([
            html.A(
                html.Img(
                    className="logo",
                    src=app.get_asset_url("images/plotly-logo-dark-theme.png"),
                ),
                href="https://plotly.com/get-demo/",
            ),
            html.A("LEARN MORE", href="https://plotly.com/dash/", target="_blank", className="demo-button")
        ], className="logos"),

        html.H2("DASH - UBER DATA APP"),
        html.P(
            """Select different days using the date picker or by selecting
            different time frames on the histogram."""
        ),
        
        dcc.DatePickerSingle(
            id="date-picker",
            min_date_allowed=dt(2014, 4, 1),
            max_date_allowed=dt(2014, 9, 30),
            initial_visible_month=dt(2014, 4, 1),
            date=dt(2014, 4, 1).date(),
            display_format="MMMM D, YYYY",
            style={"border": "0px solid black"},
            className="div-for-dropdown",
        ),

        dbc.Row([
            # Dropdown for locations on map
            dbc.Col(
                dcc.Dropdown(
                    id="location-dropdown",
                    options=[{'label': k, 'value': k} for k, v in dict_of_locations.items()],
                    placeholder="Select a location",
                    className="div-for-dropdown",
                ),
                sm=6,
            ),

            # Dropdown to select times
            dbc.Col(
                dcc.Dropdown(
                    id="bar-selector",
                    options=[
                        {
                            "label": str(n) + ":00",
                            "value": str(n),
                        }
                        for n in range(24)
                    ],
                    multi=True,
                    placeholder="Select certain hours",
                    className="div-for-dropdown",
                ),
                sm=6,
            ),
        ]),

        html.P(id="total-rides"),
        html.P(id="total-rides-selection"),
        html.P(id="date-value"),
    ],
)