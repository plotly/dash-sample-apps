# -*- coding: utf-8 -*-
# Import required libraries
import pandas as pd
import numpy as np
import dash
import pathlib
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html


# Setup the app
app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)
app.title = "Yield Curve Analysis"
server = app.server

# get relative data folder
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()


app.layout = html.Div(
    [
        dcc.Store(id="click-output"),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.Img(
                                    src=app.get_asset_url("dash-logo.png"),
                                    className="plotly-logo",
                                )
                            ]
                        ),
                        dcc.Markdown(
                            """
                ### A View of a Chart That Predicts The Economic Future:
                
                ### The Yield Curve
                """.replace(
                                "  ", ""
                            ),
                            className="title",
                        ),
                        dcc.Markdown(
                            """This interactive report is a rendition of a
                [New York Times original](https://www.nytimes.com/interactive/2015/03/19/upshot/3d-yield-curve-economic-growth.html).""".replace(
                                "  ", ""
                            ),
                            className="subtitle",
                        ),
                    ]
                ),
                html.Div(
                    [
                        html.A(
                            html.Button("Learn More", className="learn-more-button"),
                            href="https://plot.ly/dash/pricing/",
                            target="_blank",
                        )
                    ],
                    className="info-button",
                ),
                html.Div(
                    [
                        dcc.Slider(
                            min=0,
                            max=5,
                            value=0,
                            marks={i: "".format(i + 1) for i in range(6)},
                            id="slider",
                        )
                    ],
                    className="timeline-slider",
                ),
                html.Div(
                    [
                        html.Button(
                            "Back",
                            id="back",
                            style={"display": "inline-block"},
                            n_clicks=0,
                        ),
                        html.Button(
                            "Next",
                            id="next",
                            style={"display": "inline-block", "marginLeft": "10px"},
                            n_clicks=0,
                        ),
                    ],
                    className="page-buttons",
                ),
            ],
            className="four columns sidebar",
        ),
        html.Div(
            [
                html.Div([dcc.Markdown(id="text")], className="text-box"),
                dcc.Graph(id="graph", style={"margin": "0px 20px", "height": "45vh"}),
            ],
            id="page",
            className="eight columns",
        ),
    ],
    className="row flex-display",
    style={"height": "100vh"},
)

df = pd.read_csv(DATA_PATH.joinpath("yield_curve.csv"))

xlist = list(df["x"].dropna())
ylist = list(df["y"].dropna())

del df["x"]
del df["y"]

zlist = []
for row in df.iterrows():
    index, data = row
    zlist.append(data.tolist())

UPS = {
    0: dict(x=0, y=0, z=1),
    1: dict(x=0, y=0, z=1),
    2: dict(x=0, y=0, z=1),
    3: dict(x=0, y=0, z=1),
    4: dict(x=0, y=0, z=1),
    5: dict(x=0, y=0, z=1),
}

CENTERS = {
    0: dict(x=0.3, y=0.8, z=-0.5),
    1: dict(x=0, y=0, z=-0.37),
    2: dict(x=0, y=1.1, z=-1.3),
    3: dict(x=0, y=-0.7, z=0),
    4: dict(x=0, y=-0.2, z=0),
    5: dict(x=-0.11, y=-0.5, z=0),
}

EYES = {
    0: dict(x=2.7, y=2.7, z=0.3),
    1: dict(x=0.01, y=3.8, z=-0.37),
    2: dict(x=1.3, y=3, z=0),
    3: dict(x=2.6, y=-1.6, z=0),
    4: dict(x=3, y=-0.2, z=0),
    5: dict(x=-0.1, y=-0.5, z=2.66),
}

TEXTS = {
    0: """
    ##### Yield curve 101
    The yield curve shows how much it costs the federal government to borrow
    money for a given amount of time, revealing the relationship between long-
    and short-term interest rates.
    
    It is, inherently, a forecast for what the economy holds in the future —
    how much inflation there will be, for example, and how healthy growth will
    be over the years ahead — all embodied in the price of money today,
    tomorrow and many years from now.
    """.replace(
        "  ", ""
    ),
    1: """
    ##### Where we stand
    On Wednesday, both short-term and long-term rates were lower than they have
    been for most of history – a reflection of the continuing hangover
    from the financial crisis.
    
    The yield curve is fairly flat, which is a sign that investors expect
    mediocre growth in the years ahead.
    """.replace(
        "  ", ""
    ),
    2: """
    ##### Deep in the valley
    In response to the last recession, the Federal Reserve has kept short-term
    rates very low — near zero — since 2008. (Lower interest rates stimulate
    the economy, by making it cheaper for people to borrow money, but also
    spark inflation.)
    
    Now, the Fed is getting ready to raise rates again, possibly as early as
    June.
    """.replace(
        "  ", ""
    ),
    3: """
    ##### Last time, a puzzle
    The last time the Fed started raising rates was in 2004. From 2004 to 2006,
    short-term rates rose steadily.
    
    But long-term rates didn't rise very much.
    
    The Federal Reserve chairman called this phenomenon a “conundrum," and it
    raised questions about the ability of the Fed to guide the economy.
    Part of the reason long-term rates failed to rise was because of strong
    foreign demand.
    """.replace(
        "  ", ""
    ),
    4: """
    ##### Long-term rates are low now, too
    Foreign buyers have helped keep long-term rates low recently, too — as have
    new rules encouraging banks to hold government debt and expectations that
    economic growth could be weak for a long time.
    
    The 10-year Treasury yield was as low as it has ever been in July 2012 and
    has risen only modestly since.
    Some economists refer to the economic pessimism as “the new normal.”
    """.replace(
        "  ", ""
    ),
    5: """
    ##### Long-term rates are low now, too
    Here is the same chart viewed from above.
    """.replace(
        "  ", ""
    ),
}

ANNOTATIONS = {
    0: [],
    1: [
        dict(
            showarrow=False,
            x="1-month",
            y="2015-03-18",
            z=0.046,
            text="Short-term rates basically <br>follow the interest rates set <br>by the Federal Reserve.",
            xref="x",
            yref="y",
            zref="z",
            xanchor="left",
            yanchor="auto",
        )
    ],
    2: [],
    3: [],
    4: [],
    5: [],
}


# Make 3d graph
@app.callback(Output("graph", "figure"), [Input("slider", "value")])
def make_graph(value):

    if value is None:
        value = 0

    if value in [0, 2, 3]:
        z_secondary_beginning = [z[1] for z in zlist if z[0] == "None"]
        z_secondary_end = [z[0] for z in zlist if z[0] != "None"]
        z_secondary = z_secondary_beginning + z_secondary_end
        x_secondary = ["3-month"] * len(z_secondary_beginning) + ["1-month"] * len(
            z_secondary_end
        )
        y_secondary = ylist
        opacity = 0.7

    elif value == 1:
        x_secondary = xlist
        y_secondary = [ylist[-1] for i in xlist]
        z_secondary = zlist[-1]
        opacity = 0.7

    elif value == 4:
        z_secondary = [z[8] for z in zlist]
        x_secondary = ["10-year" for i in z_secondary]
        y_secondary = ylist
        opacity = 0.25

    if value in range(0, 5):

        trace1 = dict(
            type="surface",
            x=xlist,
            y=ylist,
            z=zlist,
            hoverinfo="x+y+z",
            lighting={
                "ambient": 0.95,
                "diffuse": 0.99,
                "fresnel": 0.01,
                "roughness": 0.01,
                "specular": 0.01,
            },
            colorscale=[
                [0, "rgb(230,245,254)"],
                [0.4, "rgb(123,171,203)"],
                [0.8, "rgb(40,119,174)"],
                [1, "rgb(37,61,81)"],
            ],
            opacity=opacity,
            showscale=False,
            zmax=9.18,
            zmin=0,
            scene="scene",
        )

        trace2 = dict(
            type="scatter3d",
            mode="lines",
            x=x_secondary,
            y=y_secondary,
            z=z_secondary,
            hoverinfo="x+y+z",
            line=dict(color="#444444"),
        )

        data = [trace1, trace2]

    else:

        trace1 = dict(
            type="contour",
            x=ylist,
            y=xlist,
            z=np.array(zlist).T,
            colorscale=[
                [0, "rgb(230,245,254)"],
                [0.4, "rgb(123,171,203)"],
                [0.8, "rgb(40,119,174)"],
                [1, "rgb(37,61,81)"],
            ],
            showscale=False,
            zmax=9.18,
            zmin=0,
            line=dict(smoothing=1, color="rgba(40,40,40,0.15)"),
            contours=dict(coloring="heatmap"),
        )

        data = [trace1]

    layout = dict(
        autosize=True,
        font=dict(size=12, color="#CCCCCC"),
        margin=dict(t=5, l=5, b=5, r=5),
        showlegend=False,
        hovermode="closest",
        scene=dict(
            aspectmode="manual",
            aspectratio=dict(x=2, y=5, z=1.5),
            camera=dict(up=UPS[value], center=CENTERS[value], eye=EYES[value]),
            annotations=[
                dict(
                    showarrow=False,
                    y="2015-03-18",
                    x="1-month",
                    z=0.046,
                    text="Point 1",
                    xanchor="left",
                    xshift=10,
                    opacity=0.7,
                ),
                dict(
                    y="2015-03-18",
                    x="3-month",
                    z=0.048,
                    text="Point 2",
                    textangle=0,
                    ax=0,
                    ay=-75,
                    font=dict(color="black", size=12),
                    arrowcolor="black",
                    arrowsize=3,
                    arrowwidth=1,
                    arrowhead=1,
                ),
            ],
            xaxis={
                "showgrid": True,
                "title": "",
                "type": "category",
                "zeroline": False,
                "categoryorder": "array",
                "categoryarray": list(reversed(xlist)),
            },
            yaxis={"showgrid": True, "title": "", "type": "date", "zeroline": False},
        ),
    )

    figure = dict(data=data, layout=layout)
    return figure


# Make annotations
@app.callback(Output("text", "children"), [Input("slider", "value")])
def make_text(value):
    if value is None:
        value = 0

    return TEXTS[value]


# Button controls
@app.callback(
    [Output("slider", "value"), Output("click-output", "data")],
    [Input("back", "n_clicks"), Input("next", "n_clicks")],
    [State("slider", "value"), State("click-output", "data")],
)
def advance_slider(back, nxt, slider, last_history):

    try:
        if back > last_history["back"]:
            last_history["back"] = back
            return max(0, slider - 1), last_history

        if nxt > last_history["next"]:
            last_history["next"] = nxt
            return min(5, slider + 1), last_history

    # if last_history store is None
    except Exception as error:
        last_history = {"back": 0, "next": 0}
        return slider, last_history


# Run the Dash app
if __name__ == "__main__":
    app.run_server(debug=True)
