import pathlib
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.plotly as py
import pandas as pd
import numpy as np

from dash.dependencies import Input, Output, State
from plotly import graph_objs as go
from plotly.graph_objs import *

app = dash.Dash(
    "UberApp",
    meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)
server = app.server


# Plotly mapbox public token
mapbox_access_token = "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNqdnBvNDMyaTAxYzkzeW5ubWdpZ2VjbmMifQ.TXcBE-xg9BFdV2ocecc_7g"

# Dictionary of important locations in New York
list_of_locations = {
    "Madison Square Garden": {"lat": 40.7505, "lon": -73.9934},
    "Yankee Stadium": {"lat": 40.8296, "lon": -73.9262},
    "Empire State Building": {"lat": 40.7484, "lon": -73.9857},
    "New York Stock Exchange": {"lat": 40.7069, "lon": -74.0113},
    "JFK Airport": {"lat": 40.644987, "lon": -73.785607},
    "Grand Central Station": {"lat": 40.7527, "lon": -73.9772},
    "Times Square": {"lat": 40.7589, "lon": -73.9851},
    "Columbia University": {"lat": 40.8075, "lon": -73.9626},
    "United Nations HQ": {"lat": 40.7489, "lon": -73.9680},
}

# Initialize Dataframes that the app will use
def initialize():
    df = pd.read_csv("https://www.dropbox.com/s/vxe7623o7eqbe6n/output.csv?dl=1")
    df.drop("Unnamed: 0", 1, inplace=True)
    df["Date/Time"] = pd.to_datetime(df["Date/Time"], format="%Y-%m-%d %H:%M:%S")
    df.index = df["Date/Time"]
    df.drop("Date/Time", 1, inplace=True)
    df.drop("Base", 1, inplace=True)
    totalList = []
    for month in df.groupby(df.index.month):
        dailyList = []
        for day in month[1].groupby(month[1].index.day):
            dailyList.append(day[1])
        totalList.append(dailyList)
    return np.array(totalList)


# Layout of Dash App
app.layout = html.Div(
    className="",
    children=[
        html.Div(
            className="row",
            children=[
                # Column for user controls
                html.Div(
                    className="four columns div-user-controls",
                    children=[
                        html.Img(className="logo", src=app.get_asset_url('dash-logo-stripe.png')),
                        html.H2("DASH - UBER DATA APP"),
                        html.P(
                            """Select different days using the dropdown and the slider
                    below or by selecting different time frames on the
                    histogram."""
                        ),
                        html.P(
                            "Select any of the bars on the histogram to section data by time."
                        ),
                        html.Div(
                            className="div-for-dropdown",
                            children=[
                                # Dropdown for Months
                                dcc.Dropdown(
                                    id="my-dropdown",
                                    options=[
                                        {"label": "April 2014", "value": "Apr"},
                                        {"label": "May 2014", "value": "May"},
                                        {"label": "June 2014", "value": "June"},
                                        {"label": "July 2014", "value": "July"},
                                        {"label": "Aug 2014", "value": "Aug"},
                                        {"label": "Sept 2014", "value": "Sept"},
                                    ],
                                    value="Apr",
                                    clearable=False,
                                )
                            ],
                        ),
                        html.Div(
                            className="div-for-dropdown",
                            children=[
                                # Dropdown for days in a month
                                dcc.Dropdown(
                                    id="day-dropdown", clearable=False, value=1
                                )
                            ],
                        ),
                        html.Div(
                            className="div-for-dropdown",
                            children=[
                                # Dropdown for locations on map
                                dcc.Dropdown(
                                    id="location-dropdown",
                                    options=[
                                        {"label": i, "value": i}
                                        for i in list_of_locations
                                    ],
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
                                        {"label": "0:00", "value": "0"},
                                        {"label": "1:00", "value": "1"},
                                        {"label": "2:00", "value": "2"},
                                        {"label": "3:00", "value": "3"},
                                        {"label": "4:00", "value": "4"},
                                        {"label": "5:00", "value": "5"},
                                        {"label": "6:00", "value": "6"},
                                        {"label": "7:00", "value": "7"},
                                        {"label": "8:00", "value": "8"},
                                        {"label": "9:00", "value": "9"},
                                        {"label": "10:00", "value": "10"},
                                        {"label": "11:00", "value": "11"},
                                        {"label": "12:00", "value": "12"},
                                        {"label": "13:00", "value": "13"},
                                        {"label": "14:00", "value": "14"},
                                        {"label": "15:00", "value": "15"},
                                        {"label": "16:00", "value": "16"},
                                        {"label": "17:00", "value": "17"},
                                        {"label": "18:00", "value": "18"},
                                        {"label": "19:00", "value": "19"},
                                        {"label": "20:00", "value": "20"},
                                        {"label": "21:00", "value": "21"},
                                        {"label": "22:00", "value": "22"},
                                        {"label": "23:00", "value": "23"},
                                    ],
                                    multi=True,
                                    placeholder="Select certain hours",
                                )
                            ],
                        ),
                        html.P(id="total-rides"),
                        html.P(id="total-rides-selection"),
                        html.P(id="date-value"),
                        dcc.Markdown(
                            children=[
                                "Source: [FiveThirtyEight](https://github.com/fivethirtyeight/uber-tlc-foil-response/tree/master/uber-trip-data)"
                            ]
                        ),
                    ],
                ),
                # Column for app graphs and plots
                html.Div(
                    className="eight columns div-for-charts",
                    children=[dcc.Graph(id="map-graph"), dcc.Graph(id="histogram")],
                ),
            ],
        )
    ],
)

# Gets the amount of days in the specified month
def getValue(value):
    val = {"Apr": 30, "May": 31, "June": 30, "July": 31, "Aug": 31, "Sept": 30}[value]
    return val


# Get index for the specified month in the dataframe
def getIndex(value):
    if value == None:
        return 0
    val = {"Apr": 0, "May": 1, "June": 2, "July": 3, "Aug": 4, "Sept": 5}[value]
    return val


# Get the amount of rides per hour based on the time selected
# This also higlights the color of the histogram bars based on
# if the hours are selected
def get_selection(value, day_value, selection):
    xVal = []
    yVal = []
    xSelected = []
    colorVal = [
        "#F4EC15",
        "#DAF017",
        "#BBEC19",
        "#9DE81B",
        "#80E41D",
        "#66E01F",
        "#4CDC20",
        "#34D822",
        "#24D249",
        "#25D042",
        "#26CC58",
        "#28C86D",
        "#29C481",
        "#2AC093",
        "#2BBCA4",
        "#2BB5B8",
        "#2C99B4",
        "#2D7EB0",
        "#2D65AC",
        "#2E4EA4",
        "#2E38A4",
        "#3B2FA0",
        "#4E2F9C",
        "#603099",
    ]
    # Put selected times into a list of numbers - xSelected
    if selection is not None:
        for x in selection:
            xSelected.append(int(x))
    for i in range(0, 24, 1):
        # If bar is selected then color it white
        if i in xSelected and len(xSelected) < 24:
            colorVal[i] = "#FFFFFF"
        xVal.append(i)
        # Get the number of rides at a particular time
        yVal.append(
            len(
                totalList[getIndex(value)][day_value - 1][
                    totalList[getIndex(value)][day_value - 1].index.hour == i
                ]
            )
        )
    return [np.array(xVal), np.array(yVal), np.array(colorVal)]


# Update day-dropdown based off month in my-dropdown
@app.callback(
    [Output("day-dropdown", "options"), Output("day-dropdown", "value")],
    [Input("my-dropdown", "value")],
)
def update_day_dropdown(month):
    opts = []
    for i in range(1, getValue(month) + 1, 1):
        opts.append({"label": i, "value": i})
    return opts, 1


# Selected Data in the Histogram updates the Values in the Hours Dropdown menu
@app.callback(Output("bar-selector", "value"), 
              [Input("histogram", "selectedData"),
              Input("histogram", "clickData")])
def update_bar_selector(value, clickData):
    holder = []
    if clickData:
        holder.append(str(int(clickData["points"][0]["x"])))
    if value:
        for x in value["points"]:
            holder.append(str(int(x["x"])))
    return list(set(holder))

# Clear Selected Data if Click Data is used
@app.callback(Output("histogram", "selectedData"),
              [Input("histogram", "clickData")])
def update_selected_data(clickData):
    if(clickData):
        return {"points":[]}


# Update the total number of rides Tag
@app.callback(
    Output("total-rides", "children"),
    [Input("my-dropdown", "value"), Input("day-dropdown", "value")],
)
def update_total_rides(month_value, day_value):
    return "Total Number of rides: {:,d}".format(
        len(totalList[getIndex(month_value)][day_value - 1])
    )


# Update the total number of rides in selected times
@app.callback(
    Output("total-rides-selection", "children"),
    [
        Input("my-dropdown", "value"),
        Input("day-dropdown", "value"),
        Input("bar-selector", "value"),
    ],
)
def update_total_rides_selection(month_value, day_value, selection):
    if selection is None or len(selection) is 0:
        return ""
    totalInSelction = 0
    for x in selection:
        totalInSelction += len(
            totalList[getIndex(month_value)][day_value - 1][
                totalList[getIndex(month_value)][day_value - 1].index.hour == int(x)
            ]
        )
    return "Total rides in selection: {:,d}".format(totalInSelction)


# Update Range of List
@app.callback(
    Output("date-value", "children"),
    [
        Input("my-dropdown", "value"),
        Input("day-dropdown", "value"),
        Input("bar-selector", "value"),
    ],
)
def update_date(month_value, day_value, selection):
    holder = []
    if (
        month_value is None
        or selection is None
        or len(selection) is 24
        or len(selection) is 0
    ):
        return (month_value, " ", day_value, " - showing hour(s): All")

    for x in selection:
        holder.append(int(x))
    holder.sort()

    if holder[len(holder) - 1] - holder[0] + 2 == len(holder) + 1 and len(holder) > 2:
        return (
            month_value,
            " ",
            day_value,
            " - showing hour(s): ",
            holder[0],
            "-",
            holder[len(holder) - 1],
        )

    x = ""
    for h in holder:
        if holder.index(h) == (len(holder) - 1):
            x += str(h)
        else:
            x += str(h) + ", "
    return (month_value, " ", day_value, " - showing hour(s): ", x)


# Update Histogram Figure based on Month, Day and Times Chosen
@app.callback(
    Output("histogram", "figure"),
    [
        Input("my-dropdown", "value"),
        Input("day-dropdown", "value"),
        Input("bar-selector", "value"),
    ],
)
def update_histogram(value, day_value, selection):

    [xVal, yVal, colorVal] = get_selection(value, day_value, selection)

    layout = go.Layout(
        bargap=0.01,
        bargroupgap=0,
        barmode="group",
        margin=go.layout.Margin(l=10, r=0, t=0, b=30),
        showlegend=False,
        plot_bgcolor="#323130",
        paper_bgcolor="#323130",
        height=250,
        dragmode="select",
        font=dict(color="white"),
        xaxis=dict(
            range=[-0.5, 23.5],
            showgrid=False,
            nticks=25,
            fixedrange=True,
            ticksuffix=":00",
        ),
        yaxis=dict(
            range=[0, max(yVal) + max(yVal) / 4],
            showticklabels=False,
            showgrid=False,
            fixedrange=True,
            rangemode="nonnegative",
            zeroline=False,
        ),
        annotations=[
            dict(
                x=xi,
                y=yi,
                text=str(yi),
                xanchor="center",
                yanchor="bottom",
                showarrow=False,
                font=dict(color="white"),
            )
            for xi, yi in zip(xVal, yVal)
        ],
    )

    return go.Figure(
        data=[
            go.Bar(x=xVal, y=yVal, marker=dict(color=colorVal), hoverinfo="x"),
            go.Scatter(
                opacity=0,
                x=xVal,
                y=yVal / 2,
                hoverinfo="none",
                mode="markers",
                marker=dict(color="rgb(66, 134, 244, 0)", symbol="square", size=40),
                visible=True,
            ),
        ],
        layout=layout,
    )


# Get Color of point on scattermapbox based on date and time
def get_lat_lon_color(selectedData, month_value, day_value):
    listStr = "totalList[getIndex(month_value)][day_value-1]"
    if selectedData is None or len(selectedData) is 0:
        return listStr
    elif (
        int(selectedData[len(selectedData) - 1]) - int(selectedData[0]) + 2
        == len(selectedData) + 1
        and len(selectedData) > 2
    ):
        listStr += (
            "[(totalList[getIndex(month_value)][day_value-1].index.hour>"
            + str(int(selectedData[0]))
            + ") & \
                    (totalList[getIndex(month_value)][day_value-1].index.hour<"
            + str(int(selectedData[len(selectedData) - 1]))
            + ")]"
        )
    else:
        listStr += "["
        for point in selectedData:
            if selectedData.index(point) is not len(selectedData) - 1:
                listStr += (
                    "(totalList[getIndex(month_value)][day_value-1].index.hour=="
                    + str(int(point))
                    + ") | "
                )
            else:
                listStr += (
                    "(totalList[getIndex(month_value)][day_value-1].index.hour=="
                    + str(int(point))
                    + ")]"
                )

    return listStr


# Update Map Graph based on dropdown, slider, selected data on histogram and location dropdown
@app.callback(
    Output("map-graph", "figure"),
    [
        Input("my-dropdown", "value"),
        Input("day-dropdown", "value"),
        Input("bar-selector", "value"),
        Input("location-dropdown", "value"),
    ],
)
def update_graph(month_value, day_value, selectedData, selectedLocation):
    zoom = 12.0
    latInitial = 40.7272
    lonInitial = -73.991251
    bearing = 0

    if selectedLocation:
        zoom = 15.0
        latInitial = list_of_locations[selectedLocation]["lat"]
        lonInitial = list_of_locations[selectedLocation]["lon"]

    listStr = get_lat_lon_color(selectedData, month_value, day_value)

    return go.Figure(
        data=[
            # Data for all rides based on date and time
            Scattermapbox(
                lat=eval(listStr)["Lat"],
                lon=eval(listStr)["Lon"],
                mode="markers",
                hoverinfo="lat+lon+text",
                text=eval(listStr).index.hour,
                marker=dict(
                    showscale=True,
                    color=np.append(np.insert(eval(listStr).index.hour, 0, 0), 23),
                    opacity=0.5,
                    size=5,
                    colorscale=[
                        [0, "#F4EC15"],
                        [0.04167, "#DAF017"],
                        [0.0833, "#BBEC19"],
                        [0.125, "#9DE81B"],
                        [0.1667, "#80E41D"],
                        [0.2083, "#66E01F"],
                        [0.25, "#4CDC20"],
                        [0.292, "#34D822"],
                        [0.333, "#24D249"],
                        [0.375, "#25D042"],
                        [0.4167, "#26CC58"],
                        [0.4583, "#28C86D"],
                        [0.50, "#29C481"],
                        [0.54167, "#2AC093"],
                        [0.5833, "#2BBCA4"],
                        [1.0, "#613099"],
                    ],
                    colorbar=dict(
                        title="Time of<br>Day",
                        x=0.93,
                        xpad=0,
                        nticks=24,
                        tickfont=dict(color="#d8d8d8"),
                        titlefont=dict(color="#d8d8d8"),
                        thicknessmode="pixels",
                    ),
                ),
            ),
            # Plot of important locations on the map
            Scattermapbox(
                lat=[list_of_locations[i]["lat"] for i in list_of_locations],
                lon=[list_of_locations[i]["lon"] for i in list_of_locations],
                mode="markers",
                hoverinfo="text",
                text=[i for i in list_of_locations],
                marker=dict(size=8, color="#ffa0a0"),
            ),
        ],
        layout=Layout(
            autosize=True,
            height=500,
            margin=go.layout.Margin(l=0, r=0, t=0, b=0),
            showlegend=False,
            mapbox=dict(
                accesstoken=mapbox_access_token,
                center=dict(lat=latInitial, lon=lonInitial),  # 40.7272  # -73.991251
                style="dark",
                bearing=bearing,
                zoom=zoom,
            ),
            updatemenus=[
                dict(
                    buttons=(
                        [
                            dict(
                                args=[
                                    {
                                        "mapbox.zoom": 12,
                                        "mapbox.center.lon": "-73.991251",
                                        "mapbox.center.lat": "40.7272",
                                        "mapbox.bearing": 0,
                                        "mapbox.style": "dark",
                                    }
                                ],
                                label="Reset Zoom",
                                method="relayout",
                            )
                        ]
                    ),
                    direction="left",
                    pad={"r": 0, "t": 0, "b": 0, "l": 0},
                    showactive=False,
                    type="buttons",
                    x=0.45,
                    y=0.02,
                    xanchor="left",
                    yanchor="bottom",
                    bgcolor="#323130",
                    borderwidth=1,
                    bordercolor="#6d6d6d",
                    font=dict(color="#FFFFFF"),
                )
            ],
        ),
    )


@app.server.before_first_request
def defineTotalList():
    global totalList
    totalList = initialize()


if __name__ == "__main__":
    app.run_server(debug=True)
