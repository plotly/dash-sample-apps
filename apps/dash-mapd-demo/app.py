import dash
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import numpy as np
import pandas as pd
from dash.dependencies import State, Input, Output
from dash.exceptions import PreventUpdate
import pymapd
import datetime
from datetime import datetime as dt
import os
import time

wk_map = {
    "1": "Mon",
    "2": "Tues",
    "3": "Wed",
    "4": "Thu",
    "5": "Fri",
    "6": "Sat",
    "7": "Sun",
}
wk_map_rev_map = {"Mon": 1, "Tues": 2, "Wed": 3, "Thu": 4, "Fri": 5, "Sat": 6, "Sun": 7}

app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

server = app.server

app.config.suppress_callback_exceptions = True

user = "mapd"
password = "HyperInteractive"
db_name = "mapd"

# Get credentials from env
if "DB_HOST" in os.environ:
    host = os.environ.get("DB_HOST")
else:
    host = "localhost"

table = "flights_2008_10k"


# Connect to omnisci server
def db_connect():
    """Return a valid connection. """
    try:
        connection = pymapd.connect(
            user=user, password=password, host=host, dbname=db_name
        )

        if table not in connection.get_tables():
            print(
                "Table {} not found in this database, please load sample data.".format(
                    table
                )
            )
    except Exception as e:
        print("Error connecting to OmniSci database: {}, retry".format(e))
        time.sleep(1)
        # Make retry connection
        connection = pymapd.connect(
            user=user, password=password, host="localhost", dbname=db_name
        )

    return connection


# First-time-connection upon deployment. if fails, raise attention in app log
INIT_SELECT_QUERY = (
    f"SELECT flight_dayofweek, depdelay, arrdelay, dest_state, origin_state, dep_timestamp, arr_timestamp, deptime, airtime, "
    f"carrier_name, uniquecarrier, flightnum, origin_city, dest_city FROM {table} LIMIT 50"
)  # 14 columns

init_con = None
try:
    init_con = db_connect()
    init_df = pd.read_sql(INIT_SELECT_QUERY, init_con).dropna()
    print(init_df.head(5))
except Exception as e:
    print(
        "Initial test query failed, check your Omnisci database and re-deploy this app"
    )
finally:
    if init_con is not None:
        init_con.close()


def generate_dest_choro(dd_select, start, end):
    """
    :param dd_select: dropdown select value.
    :param start: start date from date-picker.
    :param end: end date from date-picker.

    :return: Choropleth map displaying average delay time.
    """
    start_f = f"{start} 00:00:00"
    end_f = f"{end} 00:00:00"

    state_col = "dest_state"
    if dd_select == "dep":
        state_col = "origin_state"

    choro_query = f"SELECT AVG(depdelay) AS avg_delay, {state_col} AS state FROM {table} WHERE dep_timestamp BETWEEN '{start_f}' AND '{end_f}' GROUP BY {state_col}"

    # Create new connection in callback for querying from remote, this will prevent unexpected rollback exception
    # caused by global db_connect()

    new_con = None
    try:
        new_con = db_connect()
        dest_df = pd.read_sql(choro_query, new_con).dropna()
    except Exception as e:
        print("Error querying for choropleth: ", e)
        return {}
    finally:
        if new_con:
            new_con.close()

    zmin, zmax = np.min(dest_df["avg_delay"]), np.max(dest_df["avg_delay"])

    state_data = [
        go.Choropleth(
            colorbar=dict(
                tickvals=[zmin, zmax],
                tickformat=".2f",
                ticks="",
                title="Average Delay(Min)",
                titlefont=dict(family="Open Sans", color="#515151"),
                thickness=15,
                len=0.7,
                tickcolor="#515151",
                tickfont=dict(family="Open Sans", color="#515151"),
            ),
            colorscale=[[0, "#71cde4"], [1, "#ecae50"]],
            reversescale=True,
            locations=dest_df["state"],
            z=dest_df["avg_delay"],
            locationmode="USA-states",
            marker=dict(line={"color": "rgb(255,255,255)"}),
            customdata=dest_df["state"],
        )
    ]
    title = (
        "Average Departure Delay <b>(Minutes)</b> <br> By Original State"
        if dd_select == "dep"
        else "Average Arrival Delay <b>(Minutes)</b> <br> By Destination State"
    )

    layout = dict(
        title=dict(
            text=title,
            font=dict(family="Open Sans, sans-serif", size=15, color="#515151"),
        ),
        margin=dict(l=20, r=20, b=20, pad=5),
        automargin=False,
        clickmode="event+select",
        geo=go.layout.Geo(
            scope="usa", projection=go.layout.geo.Projection(type="albers usa")
        ),
    )

    return {"data": state_data, "layout": layout}


def generate_flights_hm(state, dd_select, start, end, select=False):
    # Total flight count by Day of Week / Hour
    hm = []

    state_query = ""
    start_f = f"{start} 00:00:00"
    end_f = f"{end} 00:00:00"

    if select:
        state_query = f"origin_state = '{state}' AND "

    for i in range(0, 24):
        hm_query = (
            f"SELECT flight_dayofweek, COUNT(*) AS Time_{i} FROM {table} WHERE {state_query}{dd_select}time >= {i}00 AND {dd_select}time < {i + 1}00 "
            f"AND {dd_select}_timestamp BETWEEN '{start_f}' AND '{end_f}' group by flight_dayofweek"
        )

        con = None

        try:
            con = db_connect()
            hm_df = pd.read_sql(hm_query, con).dropna()
            hm_df = hm_df.set_index("flight_dayofweek")

            hm.append(hm_df)
        except Exception as e:
            print("Error querying for heatmap: ", e)
            return {}
        finally:
            if con is not None:
                con.close()

    hm_df = pd.concat(hm, axis=1)
    hm_df = hm_df.fillna(0).reset_index()

    y = []
    for dayofweek in hm_df["flight_dayofweek"]:
        if str(dayofweek) in wk_map:
            y.append(wk_map[str(dayofweek)])

    zmin, zmax = np.min(hm_df.to_numpy()), np.max(hm_df.to_numpy())

    trace = dict(
        type="heatmap",
        z=hm_df.to_numpy(),
        x=list("{}:00".format(i) for i in range(24)),
        y=y,
        colorscale=[[0, "#71cde4"], [1, "#ecae50"]],
        reversescale=True,
        showscale=True,
        xgap=2,
        ygap=2,
        colorbar=dict(
            len=0.7,
            ticks="",
            title="Delays",
            titlefont=dict(family="Gravitas One", color="#515151"),
            thickness=15,
            tickcolor="#515151",
            tickfont=dict(family="Open Sans, sans serif", color="#515151"),
            tickvals=[zmin, zmax],
        ),
    )

    title = "Arrival Flights by days/hours State <b>{}</b>".format(state)
    if dd_select == "dep":
        title = "Departure Flights by days/hours State <b>{}</b>".format(state)

    layout = dict(
        title=dict(
            text=title,
            font=dict(family="Open Sans, sans-serif", size=15, color="#515151"),
        ),
        font=dict(family="Open Sans, sans-serif", size=13),
        automargin=True,
    )

    return {"data": [trace], "layout": layout}


def generate_time_series_chart(state, start, end, dd_select):
    """
    :param state: state selection from choropleth.
    :param start: start date from date-picker.
    :param end: end date from date-picker.
    :param dd_select: dropdown select value.
    :return: Binned departure/arrival flight record chart.
    """
    # bin at each hour

    start_f = "'" + start + " 00:00:00" + "'"
    end_f = "'" + end + " 00:00:00" + "'"

    title = (
        f"Flights by Departure Time State {state} <br> {start} - {end}"
        if dd_select == "dep"
        else f"Flights by Arrival Time State {state} <br> {start} - {end}"
    )
    state_col = "origin_state" if dd_select == "dep" else "dest_state"
    state_query = f"{state_col} = '{state}' AND " if len(state) > 0 else ""

    ts_query_x = f"SELECT {dd_select}_timestamp AS ts_timestamp, {dd_select}time, {state_col} FROM {table} WHERE {state_query}{dd_select}_timestamp BETWEEN {start_f} AND {end_f}"

    new_con = None
    try:
        new_con = db_connect()
        df_ts = pd.read_sql(ts_query_x, new_con).dropna()
    except Exception as e:
        print("Error querying for time-series", e)
        return {}
    finally:
        if new_con:
            new_con.close()

    start_time = datetime.datetime(2008, 1, 1)

    count = {}
    for x_val in df_ts["ts_timestamp"]:
        if x_val is pd.NaT:
            continue
        elapsed_seconds = (x_val - start_time).total_seconds()

        hour = int(elapsed_seconds / 3600)

        if hour not in count.keys():
            count[hour] = 1
        else:
            count[hour] += 1

    x = []
    y = []
    for key in sorted(count.keys()):
        x_timestamp = start_time + datetime.timedelta(hours=key)
        x.append(x_timestamp)
        y.append(count[key])

    data = [go.Scatter(x=x, y=y, mode="lines", line=dict(color="#71cde4"))]
    layout = dict(
        title=dict(
            text=title,
            font=dict(family="Open Sans, sans-serif", size=15, color="#515151"),
        ),
        font=dict(family="Open Sans, sans-serif", size=13),
        hovermode="closest",
        xaxis=dict(rangeslider=dict(visible=True), yaxis=dict(title="Records")),
    )
    return {"data": data, "layout": layout}


def generate_count_chart(state, dd_select, start, end):
    """
    :param state: state selection from choropleth.
    :param dd_select: dropdown select value.
    :param start: start date from date-picker.
    :param end: end date from date-picker.

    :return: Flight count sum graph by dayofweek
    """
    select_f = "origin_state" if dd_select == "dep" else "dest_state"
    start_f = f"{start} 00:00:00"
    end_f = f"{end} 00:00:00"
    state_query = f"{select_f} = '{state}' AND " if len(state) > 0 else ""

    count_query = (
        f"SELECT flight_dayofweek, COUNT(*) AS total_count FROM {table} WHERE {state_query}{dd_select}_timestamp BETWEEN '{start_f}' AND '{end_f}' "
        f"group by flight_dayofweek"
    )

    new_con = None
    try:
        new_con = db_connect()
        df_count = pd.read_sql(count_query, new_con).dropna()
    except Exception as e:
        print("Error querying for count_chart : ", e)
        return {}
    finally:
        if new_con:
            new_con.close()

    y = []
    for dayofweek in df_count["flight_dayofweek"]:
        if str(dayofweek) in wk_map:
            y.append(wk_map[str(dayofweek)])

    data = [
        go.Bar(
            x=df_count["total_count"],
            y=y,
            orientation="h",
            marker=dict(color="#71cde4"),
        )
    ]

    layout = dict(
        title=dict(
            text="Flight counts by Days {0} State <b>{1}<b>".format(
                select_f.split("_")[0].capitalize(), state
            ),
            font=dict(size=15),
        ),
        font=dict(family="Open Sans, sans-serif", size=13),
        xaxis=dict(
            title=dict(text="Total Flight Counts", font=dict(size=12)),
            zerolinecolor="#999999",
        ),
        clickmode="event+select",
    )

    return {"data": data, "layout": layout}


def generate_city_graph(state_select, dd_select, start, end):
    """
    :param end: end date from date-picker.
    :param start: start date from date-picker.
    :param dd_select: dropdown select value.
    :param state_select: state selection from choropleth.

    :return: city delay scatter graph.
    """
    start_f = f"{start} 00:00:00"
    end_f = f"{end} 00:00:00"

    days = ["Mon", "Tues", "Wed", "Thu", "Fri", "Sat", "Sun"]

    count_df = []
    state_f = state_select if len(state_select) > 0 else "NY"
    title = (
        f"Avg Arrival Delays By Cities for State <b>{state_f}</b>"
        if dd_select == "arr"
        else f"Avg Departure Delays By Cities for State <b>{state_f}</b>"
    )
    state_col = "origin_state" if dd_select == "arr" else "dest_state"
    city_col = "origin_city" if dd_select == "arr" else "dest_city"

    new_con = None

    for i, day in enumerate(days):
        count_query = (
            f"SELECT AVG(arrdelay) AS {day}, {city_col} AS city FROM {table} WHERE {state_col} ='{state_f}' AND flight_dayofweek = {i + 1}"
            f" AND {dd_select}_timestamp BETWEEN '{start_f}' AND '{end_f}' group by {city_col}"
        )

        try:
            new_con = db_connect()
            df_city_count = pd.read_sql(count_query, new_con).dropna().set_index("city")
        except Exception as e:
            print("Error reading count queries", e)
            return {}
        finally:
            if new_con:
                new_con.close()

        count_df.append(df_city_count)

    count_df = pd.concat(count_df, axis=1, sort=True)
    data = []
    for city in count_df.index:
        customdata = list(city for _ in range(7))
        trace = go.Scatter(
            y=days,
            x=count_df.loc[city, :],
            name=city,
            customdata=customdata,
            mode="markers",
        )
        data.append(trace)

    layout = dict(
        title=dict(
            text=title,
            font=dict(family="Open Sans, sans-serif", size=15, color="#515151"),
        ),
        automargin=True,
        font=dict(family="Open Sans, sans-serif", size=13),
        xaxis=dict(
            title=dict(text="Minutes", font=dict(size=12)), zerolinecolor="#999999"
        ),
        hovermode="closest",
        clickmode="event+select",
        dragmode="select",
        showlegend=True,
    )

    return {"data": data, "layout": layout}


def query_helper(state_query, dd_select, start, end, weekday_query):
    add_and = ""
    if state_query:
        add_and = "AND"
    query = (
        f"SELECT uniquecarrier AS carrier, flightnum, dep_timestamp, arr_timestamp, origin_city, dest_city "
        f"FROM {table} WHERE {state_query} {add_and} {dd_select}_timestamp BETWEEN '{start}' AND '{end}' {weekday_query} limit 100"
    )
    new_con = None
    try:
        new_con = db_connect()
        dff = pd.read_sql(query, new_con).dropna()
        dff["flightnum"] = dff["carrier"] + dff["flightnum"].map(str)
        dff.drop(["carrier"], axis=1)
        return dff.to_dict("rows")
    except Exception as e:
        print(f"Error querying {query} for making table", e)
        raise PreventUpdate
    finally:
        if new_con:
            new_con.close()


app.layout = html.Div(
    className="container scalable",
    children=[
        html.Div(
            id="banner",
            className="banner",
            children=[
                html.H6("US Flights"),
                html.Img(src=app.get_asset_url("plotly_logo.png")),
            ],
        ),
        html.Div(
            className="app_main_content",
            children=[
                html.Div(
                    id="dropdown-select-outer",
                    children=[
                        html.Div(
                            [
                                html.P("Select Departure/Arrival"),
                                dcc.Dropdown(
                                    id="dropdown-select",
                                    options=[
                                        {"label": "Departure", "value": "dep"},
                                        {"label": "Arrival", "value": "arr"},
                                    ],
                                    value="dep",
                                ),
                            ],
                            className="selector",
                        ),
                        html.Div(
                            [
                                html.P("Select Date Range"),
                                dcc.DatePickerRange(
                                    id="date-picker-range",
                                    min_date_allowed=dt(2008, 1, 1),
                                    max_date_allowed=dt(
                                        2008, 1, 7
                                    ),  # set maximum limit according to local casting
                                    initial_visible_month=dt(2008, 1, 1),
                                    minimum_nights=3,
                                    display_format="MMM Do, YY",
                                    start_date=dt(2008, 1, 1),
                                    end_date=dt(2008, 1, 7),
                                ),
                            ],
                            id="date-picker-outer",
                            className="selector",
                        ),
                    ],
                ),
                html.Div(
                    id="top-row",
                    className="row",
                    children=[
                        html.Div(
                            id="map_geo_outer",
                            className="seven columns",
                            # avg arrival/dep delay by destination state
                            children=dcc.Graph(id="choropleth"),
                        ),
                        html.Div(
                            id="flights_by_day_hm_outer",
                            className="five columns",
                            children=dcc.Loading(children=dcc.Graph(id="flights_hm")),
                        ),
                    ],
                ),
                html.Div(
                    id="middle-row",
                    className="row",
                    children=[
                        html.Div(
                            id="Flights-by-city-outer",
                            className="six columns",
                            children=dcc.Loading(
                                children=dcc.Graph(id="value_by_city_graph")
                            ),
                        ),
                        html.Div(
                            id="time-series-outer",
                            className="six columns",
                            children=dcc.Loading(
                                children=dcc.Graph(
                                    id="flights_time_series",
                                    figure=generate_time_series_chart(
                                        "",
                                        "2018-01-01 00:00:00",
                                        "2018-01-08 00:00:00",
                                        "dep",
                                    ),
                                )
                            ),
                        ),
                    ],
                ),
                html.Div(
                    id="bottom-row",
                    className="row",
                    children=[
                        html.Div(
                            id="Count_by_days_outer",
                            className="four columns",
                            children=dcc.Loading(
                                children=dcc.Graph(id="count_by_day_graph")
                            ),
                        ),
                        html.Div(
                            id="flight_info_table_outer",
                            className="eight columns",
                            children=dcc.Loading(
                                id="table-loading",
                                children=dash_table.DataTable(
                                    id="flights-table",
                                    columns=[
                                        {"name": i, "id": i}
                                        for i in [
                                            "flightnum",
                                            "dep_timestamp",
                                            "arr_timestamp",
                                            "origin_city",
                                            "dest_city",
                                        ]
                                    ],
                                    filter_action="native",
                                    fill_width=True,
                                    data=[],
                                    style_as_list_view=True,
                                    style_header={
                                        "textTransform": "Uppercase",
                                        "fontWeight": "bold",
                                        "backgroundColor": "#ffffff",
                                        "padding": "10px 0px",
                                    },
                                    style_data_conditional=[
                                        {
                                            "if": {"row_index": "even"},
                                            "backgroundColor": "#f5f6f7",
                                        },
                                        {
                                            "if": {"row_index": "odd"},
                                            "backgroundColor": "#ffffff",
                                        },
                                    ],
                                ),
                            ),
                        ),
                    ],
                ),
            ],
        ),
    ],
)


@app.callback(
    Output("choropleth", "figure"),
    [
        Input("dropdown-select", "value"),
        Input("date-picker-range", "start_date"),
        Input("date-picker-range", "end_date"),
    ],
)
def update_choro(dd_select, start, end):
    # Update choropleth when dropdown or date-picker change
    if dd_select is None:
        dd_select = "dep"

    start, end = start.replace("T", " "), end.replace("T", " ")
    return generate_dest_choro(dd_select, start, end)


@app.callback(
    Output("flights-table", "data"),
    [
        Input("flights_time_series", "relayoutData"),
        Input("count_by_day_graph", "clickData"),
        Input("value_by_city_graph", "selectedData"),
        Input("choropleth", "figure"),
    ],
    [
        State("dropdown-select", "value"),
        State("date-picker-range", "start_date"),
        State("date-picker-range", "end_date"),
        State("choropleth", "clickData"),
    ],
)
def update_sel_for_table(
    ts_select, count_click, city_select, choro_fig, dd_select, start, end, choro_click
):
    """
    :return: Data for generating flight info datatable.
    """
    if dd_select is None:
        dd_select = "dep"

    start, end = start.replace("T", " "), end.replace("T", " ")
    start_f = f"{start} 00:00:00"
    end_f = f"{end} 00:00:00"

    ctx = dash.callback_context
    inputs = ctx.inputs
    prop_id = ctx.triggered[0]["prop_id"].split(".")[0]

    state_query = ""
    if choro_click is not None:
        state = choro_click["points"][0]["location"]
        state_query = f" origin_state = '{state}'"

    try:
        if prop_id == "choropleth":
            print("table triggered by choropleth figure update")
            return query_helper(state_query, dd_select, start_f, end_f, "")

        elif prop_id == "flights_time_series":
            if "xaxis.range[0]" not in ts_select:
                print("table triggered by time-series figure update but autorange-true")
                raise PreventUpdate

            range_min, range_max = (
                inputs["flights_time_series.relayoutData"]["xaxis.range[0]"],
                inputs["flights_time_series.relayoutData"]["xaxis.range[1]"],
            )

            print("table triggered by choropleth figure")
            return query_helper(
                state_query, dd_select, range_min[:-5], range_max[:-5], ""
            )

        elif prop_id == "count_by_day_graph":

            wk_day = wk_map_rev_map[
                inputs["count_by_day_graph.clickData"]["points"][0]["y"]
            ]
            wk_day_query = f"AND flight_dayofweek = {wk_day}"
            print("table triggered by count_by_day barchart")
            return query_helper(state_query, dd_select, start_f, end_f, wk_day_query)

        elif prop_id == "value_by_city_graph":
            print("table triggered by city scatterplot")
            wk_days = []
            cities = []
            for selected_point in city_select["points"]:
                city = selected_point["customdata"]
                wk_day = selected_point["y"]
                if city not in cities:
                    cities.append(city)
                if wk_day not in wk_days:
                    wk_days.append(wk_day)
            # Only update table by querying selected city.
            frames = []
            q_template = (
                "SELECT uniquecarrier AS carrier, flightnum, dep_timestamp, arr_timestamp, origin_city, dest_city "
                "FROM {} WHERE {}_timestamp BETWEEN '{}' AND '{}' AND flight_dayofweek = {} AND {} = '{}' limit 25"
            )

            city_col = "dest_city"
            if dd_select == "dep":
                city_col = "origin_city"

            new_con = None

            for wk_day in wk_days:
                for city in cities:
                    q = q_template.format(
                        table,
                        dd_select,
                        start_f,
                        end_f,
                        wk_map_rev_map[wk_day],
                        city_col,
                        city,
                    )
                    try:
                        new_con = db_connect()
                        dff = pd.read_sql(q, new_con).dropna()
                        dff["flightnum"] = dff["carrier"] + dff["flightnum"].map(str)
                        dff.drop(["carrier"], axis=1)
                        frames.append(dff)
                    except Exception as e:
                        print("Error querying for updating datatable {}".format(e))
                        return []
                    finally:
                        if new_con:
                            new_con.close()
            if len(frames) == 0:
                return []
            return pd.concat(frames).to_dict("rows")
        else:  # should not reach
            print("table triggered by none of above ids")
            raise PreventUpdate
    except Exception as e:
        raise PreventUpdate


@app.callback(
    Output("flights_hm", "figure"),
    [Input("choropleth", "clickData"), Input("choropleth", "figure")],
    [
        State("dropdown-select", "value"),
        State("date-picker-range", "end_date"),
        State("date-picker-range", "start_date"),
    ],
)
def update_hm(choro_click, choro_figure, dd_select, end, start):
    if dd_select is None:
        dd_select = "dep"

    start, end = start.replace("T", " "), end.replace("T", " ")

    if choro_click is not None:
        state = []
        for point in choro_click["points"]:
            state.append(point["location"])

        return generate_flights_hm(state[0], dd_select, start, end, select=True)
    else:
        return generate_flights_hm("", dd_select, start, end, select=False)


@app.callback(
    Output("flights_time_series", "figure"),
    [Input("choropleth", "clickData"), Input("choropleth", "figure")],
    [
        State("dropdown-select", "value"),
        State("date-picker-range", "end_date"),
        State("date-picker-range", "start_date"),
    ],
)
def update_time_series(choro_click, choro_figure, dd_select, end, start):
    # Update time-series chart based on state select
    if dd_select is None:
        dd_select = "dep"
    start, end = start.replace("T", " "), end.replace("T", " ")

    if choro_click is not None:
        state = []
        for point in choro_click["points"]:
            state.append(point["location"])

        return generate_time_series_chart(state[0], start, end, dd_select)
    else:
        return generate_time_series_chart("", start, end, dd_select)


@app.callback(
    [Output("count_by_day_graph", "figure"), Output("value_by_city_graph", "figure")],
    [Input("choropleth", "clickData"), Input("choropleth", "figure")],
    [
        State("dropdown-select", "value"),
        State("date-picker-range", "end_date"),
        State("date-picker-range", "start_date"),
    ],
)
def update_state_click(choro_click, choro_fig, dd_select, end, start):
    # Update count graph/city graph based on state select
    if dd_select is None:
        dd_select = "dep"
    start, end = start.replace("T", " "), end.replace("T", " ")

    if choro_click is not None:
        state = []
        for point in choro_click["points"]:
            state.append(point["location"])

        return (
            generate_count_chart(state[0], dd_select, start, end),
            generate_city_graph(state[0], dd_select, start, end),
        )
    else:
        return (
            generate_count_chart("", dd_select, start, end),
            generate_city_graph("", dd_select, start, end),
        )


# Run the server
if __name__ == "__main__":
    app.run_server(
        debug=True, port=8050, dev_tools_hot_reload=False, use_reloader=False
    )
