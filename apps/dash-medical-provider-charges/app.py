import dash
from dash import (
    html,
    Input,
    Output,
    State,
    callback,
)
import plotly.graph_objs as go


import pandas as pd
import os

from utils.helper_functions import (
    generate_aggregation,
    get_lat_lon_add,
    region_dropdown,
    checklist,
    procedure_stats,
)

from utils.components import build_upper_left_panel, map_card, procedure_card

from utils.figures import (
    generate_geo_map,
    generate_procedure_plot,
    hospital_datatable,
    update_geo_map,
    update_procedure_plot,
)

app = dash.Dash(
    __name__,
    meta_tags=[
        {
            "name": "viewport",
            "content": "width=device-width, initial-scale=1, maximum-scale=1.0, user-scalable=no",
        }
    ],
)
app.title = "Medical Provider Charges"
server = app.server

app.config["suppress_callback_exceptions"] = True

app.layout = html.Div(
    className="container scalable",
    children=[
        html.Div(
            id="banner",
            className="banner",
            children=[
                html.H6("Dash Clinical Analytics"),
                html.Img(src=app.get_asset_url("plotly_logo_white.png")),
            ],
        ),
        html.Div(
            id="upper-container",
            className="row",
            children=[
                build_upper_left_panel(),
                map_card(),
            ],
        ),
        procedure_card(),
    ],
)


@callback(
    Output("region-select", "value"),
    Output("region-select", "options"),
    Output("map-title", "children"),
    Input("region-select-all", "value"),
    Input("state-select", "value"),
)
def update_region_dropdown(select_all, state_select):
    return region_dropdown(select_all, state_select)


@callback(
    Output("checklist-container", "children"),
    Input("region-select", "value"),
    State("region-select", "options"),
    State("region-select-all", "value"),
)
def update_checklist(selected, select_options, checked):
    return checklist(selected, select_options, checked)


@callback(
    Output("cost-stats-container", "children"),
    Input("geo-map", "selectedData"),
    Input("procedure-plot", "selectedData"),
    Input("metric-select", "value"),
    Input("state-select", "value"),
)
def update_hospital_datatable(geo_select, procedure_select, cost_select, state_select):
    return hospital_datatable(geo_select, procedure_select, cost_select, state_select)


@callback(
    Output("procedure-stats-container", "children"),
    Input("procedure-plot", "selectedData"),
    Input("geo-map", "selectedData"),
    Input("metric-select", "value"),
    State("state-select", "value"),
)
def update_procedure_stats(procedure_select, geo_select, cost_select, state_select):
    return procedure_stats(procedure_select, geo_select, cost_select, state_select)


@callback(
    Output("geo-map", "figure"),
    Input("metric-select", "value"),
    Input("region-select", "value"),
    Input("procedure-plot", "selectedData"),
    Input("state-select", "value"),
)
def return_updated_geo_map(cost_select, region_select, procedure_select, state_select):
    return update_geo_map(cost_select, region_select, procedure_select, state_select)


@callback(
    Output("procedure-plot", "figure"),
    Input("metric-select", "value"),
    Input("region-select", "value"),
    Input("geo-map", "selectedData"),
    Input("state-select", "value"),
)
def return_updated_procedure_plot(cost_select, region_select, geo_select, state_select):
    return update_procedure_plot(cost_select, region_select, geo_select, state_select)


if __name__ == "__main__":
    app.run_server(debug=True)
