import copy
import dash
import math
import datetime as dt

from dash import (
    Dash,
    html,
    dcc,
    Input,
    Output,
    State,
    callback,
    callback_context,
    ClientsideFunction,
)

# Multi-dropdown options
from constants import (
    WELL_STATUSES,
    WELL_TYPES,
)

from utils.data import df

from utils.figures import (
    make_main_figure,
    make_individual_figure,
    make_aggregate_figure,
    make_pie_figure,
    make_count_figure,
)
from utils.helper_functions import (
    human_format,
    filter_dataframe,
    produce_individual,
    produce_aggregate,
)

from utils.components import (
    header,
    controls_card,
    top_data_cards,
    main_graph,
    individual_graph,
    pie_graph,
    aggregate_graph,
)

app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width"}],
)
app.title = "Oil & Gas Wells"
server = app.server

# Create app layout
app.layout = html.Div(
    [
        dcc.Store(id="aggregate_data"),
        # empty Div to trigger javascript file for graph resizing
        html.Div(id="output-clientside"),
        # Header to be replaced
        header(app, "black", "New York Oil and Gas", "Production Overview"),
        html.Div(
            [
                controls_card(),
                # top data cards
                top_data_cards(),
            ],
            className="row flex-display",
        ),
        html.Div(
            [
                main_graph(),
                individual_graph(),
            ],
            className="row flex-display",
        ),
        html.Div(
            [
                pie_graph(),
                aggregate_graph(),
            ],
            className="row flex-display",
        ),
    ],
    id="mainContainer",
    style={"display": "flex", "flex-direction": "column"},
)

# Create callbacks
app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="resize"),
    Output("output-clientside", "children"),
    [Input("count_graph", "figure")],
)


@callback(
    Output("aggregate_data", "data"),
    Input("well_statuses", "value"),
    Input("well_types", "value"),
    Input("year_slider", "value"),
)
def update_production_text(well_statuses, well_types, year_slider):

    dff = filter_dataframe(df, well_statuses, well_types, year_slider)
    selected = dff["API_WellNo"].values
    index, gas, oil, water = produce_aggregate(selected, year_slider)
    return [human_format(sum(gas)), human_format(sum(oil)), human_format(sum(water))]


# Radio -> multi
@callback(Output("well_statuses", "value"), Input("well_status_selector", "value"))
def display_status(selector):
    if selector == "all":
        return list(WELL_STATUSES.keys())
    elif selector == "active":
        return ["AC"]
    return []


# Radio -> multi
@callback(Output("well_types", "value"), Input("well_type_selector", "value"))
def display_type(selector):
    if selector == "all":
        return list(WELL_TYPES.keys())
    elif selector == "productive":
        return ["GD", "GE", "GW", "IG", "IW", "OD", "OE", "OW"]
    return []


# Slider -> count graph
@callback(Output("year_slider", "value"), Input("count_graph", "selectedData"))
def update_year_slider(count_graph_selected):

    if count_graph_selected is None:
        return [1990, 2010]

    nums = [int(point["pointNumber"]) for point in count_graph_selected["points"]]
    return [min(nums) + 1960, max(nums) + 1961]


# Selectors -> well text
@callback(
    Output("well_text", "children"),
    [
        Input("well_statuses", "value"),
        Input("well_types", "value"),
        Input("year_slider", "value"),
    ],
)
def update_well_text(well_statuses, well_types, year_slider):

    dff = filter_dataframe(df, well_statuses, well_types, year_slider)
    return dff.shape[0]


@app.callback(
    Output("gasText", "children"),
    Output("oilText", "children"),
    Output("waterText", "children"),
    Input("aggregate_data", "data"),
)
def update_text(data):
    return data[0] + " mcf", data[1] + " bbl", data[2] + " bbl"


# Selectors -> main graph
@app.callback(
    Output("main_graph", "figure"),
    Input("well_statuses", "value"),
    Input("well_types", "value"),
    Input("year_slider", "value"),
    State("lock_selector", "value"),
    State("main_graph", "relayoutData"),
)
def return_make_main_figure(
    well_statuses, well_types, year_slider, selector, main_graph_layout
):

    return make_main_figure(
        well_statuses, well_types, year_slider, selector, main_graph_layout
    )


# Main graph -> individual graph
@callback(Output("individual_graph", "figure"), Input("main_graph", "hoverData"))
def return_make_individual_figure(main_graph_hover):
    return make_individual_figure(main_graph_hover)


# Selectors, main graph -> aggregate graph
@app.callback(
    Output("aggregate_graph", "figure"),
    Input("well_statuses", "value"),
    Input("well_types", "value"),
    Input("year_slider", "value"),
    Input("main_graph", "hoverData"),
)
def return_make_aggregate_figure(
    well_statuses, well_types, year_slider, main_graph_hover
):
    return make_aggregate_figure(
        well_statuses, well_types, year_slider, main_graph_hover
    )


# Selectors, main graph -> pie graph
@callback(
    Output("pie_graph", "figure"),
    Input("well_statuses", "value"),
    Input("well_types", "value"),
    Input("year_slider", "value"),
)
def return_make_pie_figure(well_statuses, well_types, year_slider):
    return make_pie_figure(well_statuses, well_types, year_slider)


# Selectors -> count graph
@callback(
    Output("count_graph", "figure"),
    Input("well_statuses", "value"),
    Input("well_types", "value"),
    Input("year_slider", "value"),
)
def return_make_count_figure(well_statuses, well_types, year_slider):
    return make_count_figure(well_statuses, well_types, year_slider)


# Main
if __name__ == "__main__":
    app.run_server(debug=True)
