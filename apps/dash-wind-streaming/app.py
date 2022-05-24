import os
import pathlib
import numpy as np

import dash
from dash import html

from dash.dependencies import Input, Output, State
from db.api import get_wind_data, get_wind_data_by_id

from dash.exceptions import PreventUpdate

from utils.components import (
    Header,
    slider,
    checklist,
    left_graph,
    interval,
    right_graph_one,
    right_graph_two,
)

from utils.helper_functions import (
    gen_wind_speed,
    gen_wind_direction,
    gen_wind_histogram,
)


app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)
app.title = "Wind Speed Dashboard"

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
                # wind speed
                html.Div(
                    [
                        html.Div(
                            [html.H6("WIND SPEED (MPH)", className="graph__title")]
                        ),
                        left_graph,
                        interval,
                    ],
                    className="two-thirds column wind__speed__container",
                ),
                html.Div(
                    [
                        # histogram
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H6(
                                            "WIND SPEED HISTOGRAM",
                                            className="graph__title",
                                        )
                                    ]
                                ),
                                html.Div(
                                    slider,
                                    className="slider",
                                ),
                                html.Div(
                                    checklist,
                                    className="auto__container",
                                ),
                                right_graph_one,
                            ],
                            className="graph__container first",
                        ),
                        # wind direction
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H6(
                                            "WIND DIRECTION", className="graph__title"
                                        )
                                    ]
                                ),
                                right_graph_two,
                            ],
                            className="graph__container second",
                        ),
                    ],
                    className="one-third column histogram__direction",
                ),
            ],
            className="app__content",
        ),
    ],
    className="app__container",
)


@app.callback(
    Output("wind-speed", "figure"), [Input("wind-speed-update", "n_intervals")]
)
def return_gen_wind_speed(interval):
    """
    Generate the wind speed graph.

    :params interval: update the graph based on an interval
    """

    return gen_wind_speed()


@app.callback(
    Output("wind-direction", "figure"), Input("wind-speed-update", "n_intervals")
)
def return_gen_wind_direction(interval):
    """
    Generate the wind direction graph.

    :params interval: update the graph based on an interval
    """

    return gen_wind_direction()


@app.callback(
    Output("wind-histogram", "figure"),
    [Input("wind-speed-update", "n_intervals")],
    [
        State("wind-speed", "figure"),
        State("bin-slider", "value"),
        State("bin-auto", "value"),
    ],
)
def return_gen_wind_histogram(interval, wind_speed_figure, slider_value, auto_state):
    """
    Genererate wind histogram graph.

    :params interval: upadte the graph based on an interval
    :params wind_speed_figure: current wind speed graph
    :params slider_value: current slider value
    :params auto_state: current auto state
    """
    return gen_wind_histogram(wind_speed_figure, slider_value, auto_state)


@app.callback(
    Output("bin-auto", "value"),
    [Input("bin-slider", "value")],
    [State("wind-speed", "figure")],
)
def deselect_auto(slider_value, wind_speed_figure):
    """Toggle the auto checkbox."""

    # prevent update if graph has no data
    if "data" not in wind_speed_figure:
        raise PreventUpdate
    if not len(wind_speed_figure["data"]):
        raise PreventUpdate

    if wind_speed_figure is not None and len(wind_speed_figure["data"][0]["y"]) > 5:
        return [""]
    return ["Auto"]


@app.callback(
    Output("bin-size", "children"),
    [Input("bin-auto", "value")],
    [State("bin-slider", "value")],
)
def show_num_bins(autoValue, slider_value):
    """Display the number of bins."""

    if "Auto" in autoValue:
        return "# of Bins: Auto"
    return "# of Bins: " + str(int(slider_value))


if __name__ == "__main__":
    app.run_server(debug=True)
