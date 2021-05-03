import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from dash.dependencies import ClientsideFunction, Input, Output, State
import plotly.express as px
from ipywidgets import widgets
from helper.plots import custom_dims_plot

import urllib.request, json

# name mapper
with open("mapping_dict_final_binary.json", "r") as f:
    dimension_mapper_binary = json.load(f)

with open("mapping_dict_final.json", "r") as f:
    dimension_mapper = json.load(f)

city_info_bin = pd.read_csv("data_bool_geo_final.csv")

filters = []
for k, v in dimension_mapper_binary.items():
    filters.append(dict(label=v, value=k))

temperature_df = pd.read_csv("city_temperature.csv")

# datasets needed for plots
data = pd.read_csv("df_final.csv")
city_info_num = data.copy()

colnames_to_lower = dict(
    zip(
        city_info_num.drop(columns=["City", "Country", "Lat", "Long"]).columns,
        map(
            str.lower,
            city_info_num.drop(columns=["City", "Country", "Lat", "Long"]).columns,
        ),
    )
)
city_info_num.rename(columns=colnames_to_lower, inplace=True)

city_info_num_agg = city_info_num.drop(columns=["City", "Country"]).apply(np.median)

filters_layout = html.Div(
    [
        html.Div(
            [
                html.H3("Select your filters", style={"display": "inline"}),
                html.Span(
                    [html.Span(className="Select-arrow", title="is_open")],
                    className="Select-arrow-zone",
                    id="select_filters_arrow",
                ),
            ],
        ),
        html.Div(
            [
                html.P("Applied filters:", id="preferencesText"),
                dcc.Dropdown(
                    placeholder="Select Filters",
                    id="filters_drop",
                    options=filters,
                    clearable=False,
                    className="dropdownMenu",
                    multi=True,
                ),
            ],
            id="dropdown_menu_applied_filters",
        ),
    ],
    id="filters_container",
    style={"display": "block"},
    className="stack-top col-3",
)

initial_popup_layout = html.Div(
    [
        html.H1("Where to live!", className="title"),
        html.H3(
            "In this dashboard you can indicate your preferences and navigate the map to find the perfect city \
    for you to live."
        ),
        html.H3("Instructions:"),
        html.P(
            "First, select your preferences using the filters on the top left corner."
        ),
        html.P(
            "Then navigate through the map and click on the locations to see more details."
        ),
        html.H3("Click anywhere to start!"),
        html.Div(
            [
                html.Div(
                    [
                        html.H6("Authors:"),
                        html.P("Mario Rodríguez Ibáñez", className="author_name"),
                        html.P("Diogo Acabado", className="author_name"),
                        html.P("Doris Macean", className="author_name"),
                        html.P("Daniel Philippi", className="author_name"),
                    ]
                ),
                html.Div(
                    [
                        html.H6("Sources:"),
                        html.P(
                            "https://www.kaggle.com/stephenofarrell/cost-of-living",
                            className="source",
                        ),
                        html.P(
                            "https://www.economist.com/big-mac-index",
                            className="source",
                        ),
                        html.P("https://worldhappiness.report/", className="source"),
                        html.P(
                            "https://www.nestpick.com/millennial-city-ranking-2018/",
                            className="source",
                        ),
                        html.P(
                            "https://www.numbeo.com/quality-of-life/rankings.jsp",
                            className="source",
                        ),
                        html.P(
                            "https://www.kaggle.com/sudalairajkumar/daily-temperature-of-major-cities",
                            className="source",
                        ),
                        html.P(
                            "https://simplemaps.com/data/world-cities",
                            className="source",
                        ),
                    ]
                ),
            ],
            style={"display": "flex", "align": "right", "bottom": "10%"},
        ),
    ],
    id="initial_popup",
)

info_bar_layout = html.Div(
    [
        html.H1("Where to live!", className="title"),
        html.H3(
            "In this dashboard you can indicate your preferences and navigate the map to find the perfect city \
    for you to live",
            className="subtitle",
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.H6("Authors:"),
                        html.P(
                            "Mario Rodríguez Ibáñez - Diogo Acabado - Doris Macean - Daniel Philippi",
                            className="author_name",
                        ),
                    ]
                ),
            ],
            style={"display": "flex", "align": "right"},
        ),
    ],
    className="stack-top info_bar row",
    style={"display": "block"},
    id="info_bar",
)

selected_location_layout = html.Div(
    [
        html.Div(
            [
                html.H3("Insert selected location", id="title_selected_location"),
                html.Span("X", id="x_close_selection"),
            ]
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.H4("... Common indices and Quartiles analysis"),
                        html.H6(
                            "Bubble size is GDP per capita. Brush to highlight cities."
                        ),
                        dcc.Graph(id="bubble"),
                    ],
                    className="plot_container_child",
                ),
                html.Div(
                    [
                        html.H4("... Selected preferences"),
                        dcc.Graph(id="custom_dims_plot"),
                    ],
                    className="plot_container_child",
                ),
            ],
            className="plots_container",
        ),
        html.Div(
            [
                html.H4("... Average temperature of the city"),
                dcc.Graph(id="temperature_plot"),
            ],
            className="",
        ),
    ],
    id="selected_location",
    style={"display": "none"},
)

hovered_location_layout = html.Div(
    [html.Div([html.H3("city", id="hover_title"), dcc.Graph("radar")]),],
    id="hovered_location",
    style={"display": "none"},
)

print(dcc.__version__)  # 0.6.0 or above is required

app = dash.Dash(__name__, external_stylesheets="")

suppress_callback_exceptions = True

app.layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        initial_popup_layout,
        html.Div(
            [
                html.Div(
                    id="width", style={"display": "none"}
                ),  # Just to retrieve the width of the window
                html.Div(
                    id="height", style={"display": "none"}
                ),  # Just to retrieve the height of the window
                html.Div(
                    [
                        dcc.Graph(
                            id="map",
                            clear_on_unhover=True,
                            config={"doubleClick": "reset"},
                        )
                    ],
                    style={"width": "100%", "height": "100%"},
                    className="background-map-container",
                ),
            ],
            id="map_container",
            style={"display": "flex"},
        ),
        filters_layout,
        info_bar_layout,
        selected_location_layout,
        hovered_location_layout,
    ],
    id="page-content",
    style={"position": "relative"},
)

#################
#   Figures     #
#################
selections = set()


@app.callback(Output("initial_popup", "style"), Input("initial_popup", "n_clicks"))
def close_initial_popup(n_clicks):
    show_block = {"display": "block"}
    hide = {"display": "none"}
    if n_clicks is not None:
        return hide
    else:
        return show_block


@app.callback(
    Output("dropdown_menu_applied_filters", "style"),
    Output("select_filters_arrow", "title"),
    Input("select_filters_arrow", "n_clicks"),
    State("select_filters_arrow", "title"),
)
def toggle_applied_filters(n_clicks, state):
    style = {"display": "none"}
    if n_clicks is not None:
        if state == "is_open":
            style = {"display": "none"}
            state = "is_closed"
        else:
            style = {"display": "block"}
            state = "is_open"

    return style, state


selected_location = ""
x_close_selection_clicks = -1


@app.callback(Output("bubble", "clickData"), [Input("map", "clickData")])
def update_bubble_selection(click_map):
    point = click_map
    return point


@app.callback(
    Output("selected_location", "style"),
    Output("title_selected_location", "children"),
    Output("custom_dims_plot", "figure"),
    Output("bubble", "figure"),
    Output("temperature_plot", "figure"),
    [Input("map", "clickData")],
    Input("x_close_selection", "n_clicks"),
    [Input("filters_drop", "value")],
    [Input("bubble", "selectedData"), Input("bubble", "clickData")],
    Input("width", "n_clicks"),
    Input("height", "n_clicks"),
    [State("bubble", "figure")],
)
def update_selected_location(
    clickData,
    n_clicks,
    dims_selected,
    bubbleSelect,
    bubbleClick,
    width,
    height,
    bubbleState,
):

    global selected_location
    global x_close_selection_clicks
    location = ""

    if clickData is not None or dims_selected is not None:
        if clickData is not None:
            location = clickData["points"][0]["text"]
        if len(location) != 0:
            selected_location = location
            style = {"display": "block"}
        else:
            selected_location = ""
            location = selected_location
            style = {"display": "none"}
    else:
        style = {"display": "none"}

    if n_clicks != x_close_selection_clicks:
        style = {"display": "none"}
        selected_location = ""
        x_close_selection_clicks = n_clicks

    if bubbleSelect is not None or bubbleClick is not None or bubbleState is not None:
        bubble_fig = update_color(bubbleSelect, bubbleClick, bubbleState, width, height)
    else:
        bubble_fig = build_bubble_figure(width, height)
    return (
        style,
        "Compare " + location + " with other cities using...",
        update_custom_dims_plot(location, dims_selected, width, height),
        bubble_fig,
        update_temperature(location, width, height),
    )


def update_temperature(city, width, height):
    df = temperature_df
    row = df[df["City"] == city]

    avg = df.groupby("Month").mean().reset_index()

    trace0 = go.Scatter(
        x=row["Month"], y=row["AvgTemperature"], name=city, marker=dict(color="#f3d576")
    )

    trace1 = go.Scatter(
        x=avg["Month"],
        y=avg["AvgTemperature"],
        name="Average Temperature",
        marker=dict(color="#d1d1cf"),
    )

    layout = go.Layout(
        xaxis=dict(
            showline=True,
            linecolor="white",
            showgrid=False,
            tickmode="array",
            tickvals=[i for i in range(1, 13)],
            ticktext=[
                "Jan",
                "Feb",
                "Mar",
                "Apr",
                "May",
                "Jun",
                "Jul",
                "Aug",
                "Sep",
                "Oct",
                "Nov",
                "Dec",
            ],
        ),
        yaxis=dict(showline=True, linecolor="white", showgrid=False),
    )

    data = [trace0, trace1]

    fig = go.Figure(data=data, layout=layout)

    fig.update_yaxes(ticksuffix="°C")

    fig.update_layout(
        height=int(height * 0.2),
        width=int(width * 0.7),
        margin=dict(
            l=120,  # left margin
            r=120,  # right margin
            b=0,  # bottom margin
            t=0,  # top margin
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="Open sans", size=12, color="White"),
    )

    return fig


def update_custom_dims_plot(location, dims_selected, width, height):
    if dims_selected is None or len(dims_selected) == 0:
        dims_selected = ["tourism"]
    if len(location) == 0:
        return go.Figure()
    fig = custom_dims_plot(
        location,
        dims_selected,
        city_info_num,
        city_info_num_agg,
        dimension_mapper,
        width,
        height,
    )
    return fig


hovered_location = ""


@app.callback(
    Output("hovered_location", "style"),
    Output("radar", "figure"),
    Output("hover_title", "children"),
    [Input("map", "hoverData")],
)
def update_hovered_location(hoverData):
    global hovered_location
    location = ""
    if hoverData is not None:
        location = hoverData["points"][0]["text"]
        if location != hovered_location:
            hovered_location = location
            style = {"display": "block"}
        else:
            hovered_location = ""
            location = ""
            style = {"display": "none"}
    else:
        hovered_location = ""
        location = ""
        style = {"display": "none"}

    return style, update_radar(location), location


# radar plot to compare index values
def update_radar(city):
    # creating a subset dataframe
    df = data[
        [
            "City",
            "Cost of Living Index",
            "Purchasing Power Index",
            "Safety Index",
            "Health Care Index",
            "Pollution Index",
        ]
    ]

    # categories
    cat = df.columns[1:].tolist()

    select_df = df[df["City"] == city]

    Row_list = []
    r = []
    # Iterate over each row
    for index, rows in select_df.iterrows():
        for i in range(len(cat)):
            # Create list for the current
            r.append(rows[cat[i]])

            # append the list to the final list
        Row_list.append(r)
        Row_list = list(np.concatenate(Row_list).flat)

    fig = go.Figure()

    fig.add_trace(
        go.Barpolar(
            r=Row_list,
            theta=cat,
            name=city,
            marker_color=["rgb(243,203,70)"] * 6,
            marker_line_color="white",
            hoverinfo=["theta"] * 9,
            opacity=0.7,
            base=0,
        )
    )

    fig.add_trace(
        go.Barpolar(
            r=df.mean(axis=0).tolist(),
            theta=cat,
            name="Average",
            marker_color=["#986EA8"] * 6,
            marker_line_color="white",
            hoverinfo=["theta"] * 9,
            opacity=0.7,
            base=0,
        )
    )

    fig.update_layout(
        title="",
        font_size=12,
        margin=dict(
            l=110,  # left margin
            r=120,  # right margin
            b=0,  # bottom margin
            t=0,  # top margin
        ),
        height=150,
        width=300,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font_color="white",
        legend=dict(orientation="h",),
        polar=dict(
            bgcolor="rgba(0,0,0,0)",
            angularaxis=dict(linewidth=3, showline=False, showticklabels=True),
            radialaxis=dict(
                showline=False,
                showticklabels=False,
                linewidth=2,
                gridcolor="rgba(0,0,0,0)",
                gridwidth=2,
            ),
        ),
    )

    return fig


@app.callback(
    Output("page-content", "style"),
    Input("width", "n_clicks"),
    Input("height", "n_clicks"),
)
def set_page_size(width, height):
    return {"width": width, "height": height}


@app.callback(
    Output("map", "figure"),
    [Input("filters_drop", "value")],
    Input("width", "n_clicks"),
    Input("height", "n_clicks"),
)
def update_map(filter_list, width, height):
    fig = go.Figure()

    if filter_list is not None and len(filter_list) != 0:

        filters = []
        for f in filter_list:
            filters.append(city_info_bin[f])
        highlighted = city_info_bin.loc[
            np.all(filters, 0), ["City", "Country", "Lat", "Long"]
        ]
        not_highlighted = city_info_bin.loc[
            ~np.all(filters, 0), ["City", "Country", "Lat", "Long"]
        ]

        # Highlighted
        fig.add_trace(
            go.Scattermapbox(
                lat=highlighted.Lat,
                lon=highlighted.Long,
                text=highlighted.City,
                name="Compatible location",
                mode="markers",
                marker=go.scattermapbox.Marker(size=15, opacity=0.9, color="#F3D576",),
                hovertemplate="<extra></extra>",
            )
        )
    else:
        not_highlighted = city_info_bin

    # Not highlighted
    fig.add_trace(
        go.Scattermapbox(
            lat=not_highlighted.Lat,
            lon=not_highlighted.Long,
            text=not_highlighted.City,
            name="Incompatible location",
            mode="markers",
            marker=go.scattermapbox.Marker(size=10, opacity=0.9, color="#333333",),
            hovertemplate="<extra></extra>",
        )
    )

    mapbox_token = "pk.eyJ1IjoiZmFya2l0ZXMiLCJhIjoiY2ttaHYwZnQzMGI0cDJvazVubzEzc2lncyJ9.fczsOA4Hfgdf8_bAAZkdYQ"
    all_plots_layout = dict(
        mapbox=dict(
            style="mapbox://styles/farkites/ckn0lwfm319ae17o5jmk3ckvu",
            accesstoken=mapbox_token,
        ),
        legend=dict(
            bgcolor="rgba(51,51,51,0.6)",
            yanchor="top",
            y=0.35,
            xanchor="left",
            x=0,
            font=dict(family="Open Sans", size=15, color="white",),
        ),
        autosize=False,
        width=width,
        height=height,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        geo_bgcolor="rgba(0,0,0,0)",
    )
    fig.layout = all_plots_layout

    return fig


md = data[
    [
        "City",
        "Employment",
        "Startup",
        "Tourism",
        "Housing",
        "Logged GDP per capita",
        "Transport",
        "Health",
        "Food",
        "Internet Speed",
        "Generosity",
        "Freedom to make life choices",
        "Access to Contraception",
        "Gender Equality",
        "Immigration Tolerance",
        "LGBT Friendly",
        "Nightscene",
        "Beer",
        "Festival",
    ]
].copy()

for column in md.columns.tolist()[1:]:
    md["{column}.".format(column=column)] = pd.qcut(
        md[column].rank(method="first"), 4, labels=False
    )


# Build parcats dimensions
quartiles = [
    "Startup.",
    "Internet Speed.",
    "Gender Equality.",
    "Immigration Tolerance.",
    "LGBT Friendly.",
    "Nightscene.",
]


dimensions = []
for label in quartiles:
    tmp = go.parcats.Dimension(
        values=md[label], categoryorder="category descending", label=label
    )
    dimensions.append(tmp)

# Build colorscale
color = np.zeros(len(md), dtype="uint8")
colorscale = [[0, "gray"], [1, "rgb(243,203,70)"]]


gdp = 10 ** (md["Logged GDP per capita"] / 10)

sizeref = 2.0 * max(gdp) / (20 ** 2)

customdata = np.stack((pd.Series(md.index), md["City"], gdp.round(5) * 1000), axis=-1)

# bubble plot related indicators
def build_bubble_figure(width, height):
    # Build figure as FigureWidget
    fig = go.Figure(
        data=[
            go.Scatter(
                x=md["Freedom to make life choices"],
                y=md["Generosity"],
                text=md["City"],
                customdata=customdata,
                hovertemplate="""<extra></extra>
        <em>%{customdata[1]}</em><br>
        GDP per capita = %{customdata[2]} €""",
                marker={
                    "color": "#986EA8",
                    "sizeref": sizeref,
                    "sizemin": 0.005,
                    "sizemode": "area",
                    "size": gdp,
                },
                mode="markers",
                selected={"marker": {"color": "rgb(243,203,70)"}},
                unselected={"marker": {"opacity": 0.3}},
            ),
            go.Parcats(
                domain={"y": [0, 0.4]},
                dimensions=dimensions,
                line={
                    "colorscale": colorscale,
                    "cmin": 0,
                    "cmax": 1,
                    "color": color,
                    "shape": "hspline",
                },
            ),
        ]
    )

    fig.update_layout(
        font_color="white",
        font_size=9,
        xaxis={"title": "Freedom to make life choices"},
        yaxis={"title": "Generosity", "domain": [0.6, 1]},
        dragmode="lasso",
        hovermode="closest",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        autosize=False,
        bargap=0.35,
        height=int(height * 0.35),
        width=int(width * 0.45),
        margin=dict(
            l=15,  # left margin
            r=15,  # right margin
            b=0,  # bottom margin
            t=0,  # top margin
        ),
    )

    return fig


# Update color callback
def update_color(selectedData, clickData, fig, width, height):
    selection = None

    # Update selection based on which event triggered the update.
    trigger = dash.callback_context.triggered[0]["prop_id"]
    if trigger == "bubble.clickData":
        selection = [point["pointNumber"] for point in clickData["points"]]
    if trigger == "bubble.selectedData":
        selection = [point["pointIndex"] for point in selectedData["points"]]
    # Update scatter selection
    fig["data"][0]["selectedpoints"] = selection
    # Update parcats colors
    new_color = np.zeros(len(md), dtype="uint8")
    new_color[selection] = 1
    fig["data"][1]["line"]["color"] = new_color
    return fig


# Get window size function
app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="get_window_width"),
    Output("width", "n_clicks"),
    [Input("url", "href")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="get_window_height"),
    Output("height", "n_clicks"),
    [Input("url", "href")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="move_hover"),
    Output("hovered_location", "title"),
    [Input("map", "hoverData")],
)

server = app.server

if __name__ == "__main__":
    app.run_server(debug=True)
