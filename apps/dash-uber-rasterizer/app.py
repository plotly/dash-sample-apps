import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import numpy as np

from dash.dependencies import Input, Output
from plotly import graph_objs as go

import datashader as ds
from datashader import transfer_functions as tf
from couleurs import colorscales

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

x_range = (-74.929, -72.5)
y_range = (39.9, 42.1166)

# Helper Function 1 for generating initial plot
def gen_ds_image():
    cvs = ds.Canvas(x_range=x_range, y_range=y_range)
    agg_heatmap = cvs.points(rides_df, "Lon", "Lat")

    img = tf.shade(agg_heatmap)
    #img = tf.dynspread(img, threshold=0.95, max_px=5, shape='circle')

    return img #img.to_pil()

# Helper Function 2 for generating initial plot
def tf_to_plotly():
    rides_img = gen_ds_image()
    arr_rides = np.array(rides_img)
    z_rides = arr_rides.tolist()

    dims = len(z_rides[0]), len(z_rides)

    data = [dict(
        z=z_rides,
        x=np.linspace(x_range[0], x_range[1], dims[0]),
        y=np.linspace(y_range[0], y_range[1], dims[1]),
        colorscale=colorscales["plasma"],
        colorbar={"title": "Log(No. of Rides)"},
        # showscale=False,
        # reversescale = True,
        type='heatmap')]

    layout = dict(
        font={"color": "rgb(226, 239, 250)"},
        paper_bgcolor="rgb(38, 43, 61)",
        plot_bgcolor="rgb(38, 43, 61)",
        xaxis={"title": "Longitude",
               "constrain": "domain",
               "scaleanchor": "y",
               "scaleratio": np.cos(40.8 * np.pi / 180)},
        yaxis={"title": "Latitude",
               "constrain": "domain"},
               #"fixedrange": True # Expands to xaxis},
        margin=dict(t=0, b=0) # Expands yaxis
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
                            n_clicks=0,
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
                    style={"height": "88vh"},
                    className="item-b")
            ),
            dcc.Store(id="store")
        ], className="container")
    ]
)

# Callbacks Start ------------------------------

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
    app.run_server(debug=True, port = 8051)
