import dash
import dash_bootstrap_components as dbc
from dash import Dash, html, dcc, Input, Output, State, callback, callback_context
from dash.dependencies import Input, Output, State

from utils.components import (
    logo,
    predict_button,
    get_new_information_button,
    graphs,
    rul_estimation_indicator,
    info_box,
    active_power_display,
    blade_angle_display,
    active_power_from_wind_display,
    wind_speed_information,
    reactive_power_display,
)
from utils.figures import update_graph, display_click_data

from constants import gauge_size
import pickle

app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    external_stylesheets=[dbc.themes.SUPERHERO],
)
server = app.server
app.title = "Predictive Maintenance Dashboard"

app.layout = dbc.Container(
    fluid=True,
    children=[
        logo(app),
        dbc.Row(
            [
                dbc.Col(graphs, xs=10, md=10, lg=10, width=10),
                dbc.Col(
                    [
                        dbc.Row(
                            dbc.Col(
                                rul_estimation_indicator, xs=12, md=12, lg=12, width=12
                            )
                        ),
                        dbc.Row(dbc.Col(info_box, xs=12, md=12, lg=12, width=12)),
                        dbc.Row(
                            dbc.Col(
                                get_new_information_button,
                                xs=12,
                                md=12,
                                lg=12,
                                width=12,
                            )
                        ),
                        dbc.Row(dbc.Col(predict_button, xs=12, md=12, lg=12, width=12)),
                    ]
                ),
            ],
            justify="start",
            style={"display": "flex", "marginBottom": "-3%"},
        ),
        dbc.Row(
            [
                dbc.Col(
                    active_power_display,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
                dbc.Col(
                    active_power_from_wind_display,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
                dbc.Col(
                    reactive_power_display,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
                dbc.Col(
                    wind_speed_information,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
                dbc.Col(
                    blade_angle_display,
                    xs=gauge_size,
                    md=gauge_size,
                    lg=gauge_size,
                    width=gauge_size,
                ),
            ],
            style={"marginTop": "3%"},
        ),
    ],
)


@callback(
    Output("Main-Graph", "figure"),
    Output("rul-estimation-indicator-led", "value"),
    Output("Info-Textbox", "value"),
    Input("feature-dropdown", "value"),
    Input("date-picker", "start_date"),
    Input("date-picker", "end_date"),
    Input("get-new-info-button", "n_clicks"),
    Input("predict-button", "n_clicks"),
)
def return_updated_graph(selected_column, start_date, end_date, n_get_new_info, n_pred):
    return update_graph(selected_column, start_date, end_date, n_get_new_info, n_pred)


@callback(
    Output("active-power-information-gauge", "value"),
    Output("active-power-from-wind-information-gauge", "value"),
    Output("wind-power-information-gauge", "value"),
    Output("reactive-power-information-gauge", "value"),
    Output("blade-angle-information-gauge", "value"),
    Input("Main-Graph", "clickData"),
)
def _return_displayed_click_data(clickData):
    return display_click_data(clickData)


if __name__ == "__main__":
    app.run_server(debug=True, use_reloader=True)
