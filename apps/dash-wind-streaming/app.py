from dash import Dash, html, dcc, Input, Output, State, callback, callback_context

from constants import GRAPH_INTERVAL
from utils.components import Header, wind_speed_card, histogram_card, wind_direction_card
import utils.figures as figs

app = Dash(__name__, title="Wind Speed Dashboard")
server = app.server

app.layout = html.Div(
    [
        Header(
            app,
            "WIND SPEED STREAMING",
            "This app continually queries a SQL database and displays live charts of wind speed and wind direction.",
        ),
        html.Div(
            [
                wind_speed_card("wind-speed"),
                html.Div(
                    [
                        histogram_card("wind-histogram"),
                        wind_direction_card("wind-direction")
                    ],
                    className="one-third column histogram__direction",
                ),
            ],
            className="app__content",
        ),
        dcc.Interval(id="wind-speed-update", interval=GRAPH_INTERVAL)
    ],
    className="app__container",
)


@callback(
    Output("wind-histogram", "figure"),
    Output("wind-speed", "figure"), 
    Output("wind-direction", "figure"), 
    Input("wind-speed-update", "n_intervals"),
    State("bin-slider", "value"),
    State("bin-auto", "value"),
)
def return_gen_wind_histogram(interval, slider_value, auto_state):
    fig_wind_speed = figs.gen_wind_speed()
    fig_wind_direction =figs.gen_wind_direction()
    fig_wind_histogram = figs.gen_wind_histogram(fig_wind_speed, slider_value, auto_state)
    return fig_wind_histogram, fig_wind_speed, fig_wind_direction


@callback(
    Output("bin-auto", "value"),
    Output("bin-size", "children"),
    Input("bin-slider", "value"),
    Input("bin-auto", "value"),
    prevent_initial_call=True,
)
def deselect_auto(slider_value, bin_auto):
    triggered_id = callback_context.triggered[0]["prop_id"].split(".")[0]
    if triggered_id == "bin-slider":
        return [], "# of Bins: " + str(int(slider_value))
    else:
        return ["Auto"], "# of Bins: Auto"


if __name__ == "__main__":
    app.run_server(debug=True)
