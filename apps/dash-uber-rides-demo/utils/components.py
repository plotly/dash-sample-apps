from datetime import datetime as dt
from dash import html, dcc

from constants import list_of_locations

def controls(app):
    return html.Div(
    className="four columns div-user-controls",
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
        html.Div(
            className="div-for-dropdown",
            children=[
                dcc.DatePickerSingle(
                    id="date-picker",
                    min_date_allowed=dt(2014, 4, 1),
                    max_date_allowed=dt(2014, 9, 30),
                    initial_visible_month=dt(2014, 4, 1),
                    date=dt(2014, 4, 1).date(),
                    display_format="MMMM D, YYYY",
                    style={"border": "0px solid black"},
                )
            ],
        ),
        # Change to side-by-side for mobile layout
        html.Div(
            className="row",
            children=[
                html.Div(
                    className="div-for-dropdown",
                    children=[
                        # Dropdown for locations on map
                        dcc.Dropdown(
                            id="location-dropdown",
                            options=list_of_locations,
                            placeholder="Select a location",
                        )
                    ],
                ),
                html.Div(
                    className="div-for-dropdown",
                    children=[
                        # Dropdown to select times
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
                        )
                    ],
                ),
            ],
        ),
        html.P(id="total-rides"),
        html.P(id="total-rides-selection"),
        html.P(id="date-value"),
    ],
)