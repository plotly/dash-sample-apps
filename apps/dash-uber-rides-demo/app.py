from dash import Dash, dcc, html, Input, Output, callback

import utils.figures as figs
from constants import totalList
from utils.helper_functions import total_rides_calculation
from utils.components import controls

app = Dash(__name__, title = "New York Uber Rides")
server = app.server


# Layout of Dash App
app.layout = html.Div(
    children=[
        html.Div(
            className="row",
            children=[
                # Column for user controls
                controls(app),
                # Column for app graphs and plots
                html.Div(
                    className="eight columns div-for-charts bg-grey",
                    children=[
                        dcc.Graph(id="map-graph"),
                        html.Div(
                            className="text-padding",
                            children=[
                                "Select any of the bars on the histogram to section data by time."
                            ],
                        ),
                        dcc.Graph(id="histogram"),
                    ],
                ),
            ],
        )
    ]
)


@callback(
    Output("bar-selector", "value"),
    Input("histogram", "selectedData"), 
    Input("histogram", "clickData"),
)
def update_bar_selector(value, clickData):
    " Selected Data in the Histogram updates the Values in the Hours selection dropdown menu "
    holder = []
    if clickData:
        holder.append(str(int(clickData["points"][0]["x"])))
    if value:
        for x in value["points"]:
            holder.append(str(int(x["x"])))
    return list(set(holder))


@callback(
    Output("histogram", "selectedData"), 
    Input("histogram", "clickData"),
)
def update_selected_data(clickData):
    " Clear Selected Data if Click Data is used "
    if clickData:
        return {"points": []}


@callback(
    Output("total-rides", "children"), 
    Input("date-picker", "date"),
)
def update_total_rides(datePicked):
    " Update the total number of rides Tag "
    date_picked = dt.strptime(datePicked, "%Y-%m-%d")
    total_rides = len(totalList[date_picked.month - 4][date_picked.day - 1])
    return f"Total Number of rides: {total_rides:,d}"


@callback(
    Output("total-rides-selection", "children"), 
    Output("date-value", "children"),
    Input("date-picker", "date"), 
    Input("bar-selector", "value"),
)
def update_total_rides_selection(date_picked, bars_selected):
    " Update the total number of rides in selected times "
    return total_rides_calculation(date_picked, bars_selected)


@callback(
    Output("histogram", "figure"),
    Input("date-picker", "date"), 
    Input("bar-selector", "value"),
)
def update_histogram(date_picked, bars_selected):
    " Update Histogram Figure based on Month, Day and Times Chosen "
    return figs.histogram(date_picked, bars_selected)


@callback(
    Output("map-graph", "figure"),
    Input("date-picker", "date"),
    Input("bar-selector", "value"),
    Input("location-dropdown", "value"),
)
def update_graph(date_picked, bars_selected, location):
    " Update Map Graph based on date-picker, selected data on histogram and location dropdown "
    return figs.map(date_picked, bars_selected, location)


if __name__ == "__main__":
    app.run_server(debug=True)
