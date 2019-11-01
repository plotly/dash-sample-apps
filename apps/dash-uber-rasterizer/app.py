import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import numpy as np

from dash.dependencies import Input, Output
from plotly import graph_objs as go

import datashader as ds
from datashader import transfer_functions as tf
from datashader import reductions
from datashader.utils import lnglat_to_meters
from couleurs import colorscales, format_colorscale

# Initialize data frame
# df1 = pd.read_csv(
#     "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data1.csv",
#     dtype=object,
# )
# df2 = pd.read_csv(
#     "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data2.csv",
#     dtype=object,
# )
# df3 = pd.read_csv(
#     "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data3.csv",
#     dtype=object,
# )

df1 = pd.read_csv("/Users/Caner/Desktop/tmp/dashr-uber-rasterizer/data/df1.csv")
df2 = pd.read_csv("/Users/Caner/Desktop/tmp/dashr-uber-rasterizer/data/df2.csv")
df3 = pd.read_csv("/Users/Caner/Desktop/tmp/dashr-uber-rasterizer/data/df3.csv")

rides_df = pd.concat([df1, df2, df3], axis=0)
# Convert Lat & Lon to projected points in meters from the origin (Web Mercator coordinates)
rides_df.loc[:, 'LonM'], rides_df.loc[:, 'LatM'] = lnglat_to_meters(rides_df.Lon,rides_df.Lat)

default_x_range = (-74.929, -72.5)
default_y_range = (39.9, 42.1166)

default_colorscale = format_colorscale(colorscales["plasma"][::-1].reset_index())
default_colorbar = {"title": "Log(No. of Rides)"}
default_plot_bg = "black"

# Helper Function 2 for generating initial plot
def tf_to_plotly(
        df = rides_df,
        x_range=default_x_range,
        y_range=default_y_range,
        cs = default_colorscale,
        cb = default_colorbar,
        bg = default_plot_bg,
        scale = "log"):



    cvs = ds.Canvas(x_range=x_range, y_range=y_range)
    agg_heatmap = cvs.points(df, "Lon", "Lat")#, agg=reductions.sum("z"))
    img = tf.shade(agg_heatmap)
    arr_rides = np.array(img)
    z_rides = arr_rides.tolist()

    dims = len(z_rides[0]), len(z_rides)

    if scale == "log":
        z_rides = np.log(z_rides)

    data = [dict(
        z=z_rides,
        x=np.linspace(x_range[0], x_range[1], dims[0]),
        y=np.linspace(y_range[0], y_range[1], dims[1]),
        colorscale=cs,
        colorbar=cb,
        marker=dict(size=20),
        # showscale=False,
        # reversescale = True,
        type='heatmap')]

    layout = dict(
        font={"color": "rgb(226, 239, 250)"},
        paper_bgcolor="rgb(38, 43, 61)",
        plot_bgcolor=bg,
        xaxis={"title": "Longitude",
               "constrain": "domain",
               "automarging": True,
               "scaleanchor": "y",
               "scaleratio": np.cos(40.8 * np.pi / 180)},
        yaxis={"title": "Latitude",
               "automarging": True,
               "constrain": "domain"},
               #"fixedrange": True # Expands to xaxis},
        margin=dict(t=15, b=1, l=1, r=1, pad=0, autoexpand=True)
        #margin=dict(t=0, b=0) # Expands yaxis
    )

    fig = dict(data=data, layout=layout)
    return fig


# App Start ------------------------------

# Initiate application
app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)
server = app.server


# Create Layout Variables ------------------------------

header = html.Div(
    id="app-page-header",
    style={"width": "100%",
           "background": "#262B3D",
           "color": "#E2EFFA"},
    children=[
        html.A(
            id="dash-logo",
            children=[
                html.Img(src="assets/plotly-dash-logo.png",
                         height="36",
                         width="180",
                         style={"top": "10", "margin": "10px"})
            ],
            href="/Portal"
        ),
        html.H2("Uber NYC Rasterizer"),
        html.A(
            id="gh-link",
            children=["View on GitHub"],
            href="https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-uber-rasterizer",
            style={"color": "white", "border": "solid 1px white"}
        ),
        html.Img(src="assets/GitHub-Mark-Light-64px.png")
    ]
)

tabs = html.Div(
    dcc.Tabs(id="circos-control-tabs", value="what-is", children=[
        dcc.Tab(
            label="About",
            value="what-is",
            children=html.Div(
                id="control-tab", children=[
                    html.H4("What is Uber NYC Rasterizer?",
                            style={"font-size": "24pt",
                                   "font-weight": "200",
                                   "letter-spacing": "1px"}),
                    dcc.Markdown(
                        '''
                    This Dash app is a simple demonstration of the rasterizing capabilities of the _rasterly_ package.
                    The dataset consists of over 4.5 million observations, representing Uber rides taken in New York City in 2014.
                    In CSV format, the source data are over 165 MB in size. _rasterly_ is capable of processing datasets an order
                    of magnitude larger in similarly brisk fashion. The raster data required to produce the aggregation layers and color
                    gradients displayed here are computed efficiently enough to maintain the interactive feel of the application.
                        ''',
                        style={"padding": "5px"}),
                    dcc.Markdown(
                        "Visit the _rasterly_ package repository [here](https://github.com/plotly/rasterly) to learn more.",
                        style={"padding": "5px"}),
            ])
        ),
        dcc.Tab(
            label="Options",
            value="data",
            children=html.Div(
                className="circos-tab",
                children=[
                    html.Div(className="app-controls-block", children=[
                        html.H4("Colorscale",
                                style={"font-size": "18pt",
                                       "font-weight": "200",
                                       "letter-spacing": "1px"}),
                        html.Div(dcc.Dropdown(
                            id="cmap",
                            value="plasma",
                            options=[
                                {"label": "Plasma", "value": "plasma"},
                                {"label": "Viridis", "value": "viridis"},
                                {"label": "Blues", "value": "blue"},
                                {"label": "Magma", "value": "fire"}
                            ]
                        ), style = {"color": "white"}),
                        html.H4("Background",
                                style={"font-size": "18pt", "font-weight": "200", "letter-spacing": "1px"}),
                        html.Div(dcc.Dropdown(
                            id="background",
                            value="black",
                            options=[
                                {"label": "Black", "value": "black"},
                                {"label": "Grey", "value": "grey"},
                                {"label": "White", "value": "white"}
                            ]
                        ), style={"color": "white"}),
                        html.H4("Point Scaling", style={"font-size": "18pt", "font-weight": "200", "letter-spacing": "1px"}),
                        dcc.Dropdown(
                            id="scaling",
                            value="log",
                            options=[
                                {"label": "Log", "value": "log"},
                                {"label": "Origin", "value": "origin"}
                            ]
                        ),
                        html.H4("Reduction method", style={"font-size": "18pt", "font-weight": "200", "letter-spacing": "1px"}),
                        dcc.Dropdown(
                            id="reduc",
                            value="sum",
                            options=[
                                {"label": "sum", "value": "sum"},
                                {"label": "any", "value": "any"},
                                {"label": "mean", "value": "mean"}
                            ]
                        ),
                        html.H4("Pixel Size",
                                style={"font-size": "18pt", "font-weight": "200", "letter-spacing": "1px"}),
                        dcc.Slider(
                            id="point-size",
                            min=0,
                            max=10,
                            step=1,
                            value=0,
                            marks={str(i):str(i) for i in range(1,11)}
                        ),
                        html.Br(),
                        html.Br(),
                        html.Button(
                            id="reset-button",
                            children="Reset Graph"
                        )
                    ])

                ]

            )
        )

    ])
)

options = html.Div([tabs], className="item-a")


# Create Layout ------------------------------

app.layout = html.Div(
    children=[
        header,
        html.Div([
            options,
            html.Div(
                children=dcc.Graph(
                    id="rasterizer-output",
                    figure=tf_to_plotly(),
                    style={"height": "88vh"}),
                className="item-b"
            ),
            dcc.Store(id="store")
        ], className="container")
    ]
)

# Callbacks Start ------------------------------





# Callback to subset and return data ranges based on the zoom of the plot.

@app.callback(
    Output("store", "data"),
    [Input("rasterizer-output", "relayoutData"), Input("reset-button", "n_clicks")],
)
def update_stored_data(relayout, n_clicks):
    if n_clicks is None:
        n_clicks = 0

    if n_clicks > 0:
        x_range = (min(rides_df.loc[:, "Lon"]), -72.5)
        y_range = (39.9, max(rides_df.loc[:, "Lat"]))
    else:
        if len(relayout) == 4:
            x_range = (relayout['xaxis.range[0]'], relayout['xaxis.range[1]'])
            y_range = (relayout['yaxis.range[0]'], relayout['yaxis.range[1]'])
        else:
            x_range = (min(rides_df.loc[:, "Lon"]), -72.5)
            y_range = (39.9, max(rides_df.loc[:, "Lat"]))
    return [x_range, y_range]


# Callback to generate rasterization plot.

@app.callback(
    Output("rasterizer-output", "figure"),
    [Input("store", "data"),
     Input("cmap", "value"),
     Input("background", "value"),
     Input("reduc", "value"),
     Input("scaling", "value"),
     Input("point-size", "value")],
)


# Callback to generate rasterization plot.

def update_graph(stored_ranges, cmap, background, reduc, scale, point_size):
    if cmap == "blue":
        color = colorscales["blue"]
    elif cmap == "viridis":
        color = colorscales["viridis"]
    elif cmap == "plasma":
        color = colorscales["plasma"]
    elif cmap == "fire":
        color = colorscales["fire"]

    if background == "black":
        color = color[::-1].reset_index()
    print(color)

    color = format_colorscale(color)

    x_min = stored_ranges[0][0]
    x_max = stored_ranges[0][1]
    y_min = stored_ranges[1][0]
    y_max = stored_ranges[1][1]


    filtered_df_lat = rides_df[(rides_df["Lat"] > y_min) & (rides_df["Lat"] < y_max)]
    filtered_df_lon = filtered_df_lat[(filtered_df_lat["Lon"] > x_min) & (filtered_df_lat["Lat"] > x_max)]

    colorbar_title = {"title": "Log(No. of Rides)"} if scale == "log" else {"title": "No. of Rides"}

    return tf_to_plotly(
        df=filtered_df_lon,
        x_range=(x_min, x_max),
        y_range=(y_min, y_max),
        cs=color,
        cb=colorbar_title,
        bg=background,
        scale=scale)













# @app.callback(
#     Output("bar-selector", "value"),
#     [Input("histogram", "selectedData"), Input("histogram", "clickData")],
# )
# def update_bar_selector(value, clickData):
#     holder = []
#     if clickData:
#         holder.append(str(int(clickData["points"][0]["x"])))
#     if value:
#         for x in value["points"]:
#             holder.append(str(int(x["x"])))
#     return list(set(holder))


if __name__ == "__main__":
    app.run_server(port = 8051)
