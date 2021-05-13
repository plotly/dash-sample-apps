import os

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import geopandas as gpd
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats

from peaky_finders.predictor import (
    create_load_duration,
    ISO_LIST,
    get_peak_data,
    get_iso_map,
    get_forecasts,
)

app_name = os.getenv("APP_NAME", "dash-peaky-finders")

iso_map = get_iso_map()
peak_data = get_peak_data(ISO_LIST)
predictions, load, temperature = get_forecasts(ISO_LIST)
load_duration_curves = create_load_duration(peak_data)


TEMPLATE = "plotly_white"

app = dash.Dash(
    external_stylesheets=[dbc.themes.LUX], suppress_callback_exceptions=True
)
app.title = "US Electric Grid Forecasting"
server = app.server

"""Homepage"""
app.layout = html.Div(
    [dcc.Location(id="url", refresh=False), html.Div(id="page-content"),]
)

index_page = html.Div(
    [
        html.Br(),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(html.H1(children="Welcome to Peaky Finders"), width=5),
                dbc.Col(width=5),
            ],
            justify="center",
        ),
        html.Br(),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            html.H4(
                                children="To what extent do weather and weekday determine total electricity demand on the grid? Click an ISO button below to find out."
                            ),
                            html.Div(
                                [
                                    dcc.Link(
                                        html.Button(
                                            "HOME", id="home-button", className="mr-1"
                                        ),
                                        href=f"/{app_name}/",
                                    ),
                                    dcc.Link(
                                        html.Button(
                                            "CAISO", id="caiso-button", className="mr-1"
                                        ),
                                        href=f"/{app_name}/caiso",
                                    ),
                                    dcc.Link(
                                        html.Button(
                                            "MISO", id="miso-button", className="mr-1"
                                        ),
                                        href=f"/{app_name}/miso",
                                    ),
                                    dcc.Link(
                                        html.Button(
                                            "PJM", id="pjm-button", className="mr-1"
                                        ),
                                        href=f"/{app_name}/pjm",
                                    ),
                                    dcc.Link(
                                        html.Button(
                                            "NYISO", id="nyiso-button", className="mr-1"
                                        ),
                                        href=f"/{app_name}/nyiso",
                                    ),
                                    dcc.Link(
                                        html.Button(
                                            "ISONE", id="isone-button", className="mr-1"
                                        ),
                                        href=f"/{app_name}/isone",
                                    ),
                                ]
                            ),
                        ]
                    ),
                    width=7,
                ),
                dbc.Col(width=3),
            ],
            justify="center",
        ),
        html.Br(),
        html.Br(),
        html.Br(),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H4(children="ISO Territory Map"), width=4), dbc.Col(width=4)],
            justify="center",
        ),
        html.Div(
            [
                dcc.Graph(
                    figure=px.choropleth(
                        iso_map,
                        geojson=iso_map.geometry,
                        locations=iso_map.index,
                        color="NAME",
                        projection="mercator",
                    )
                    .update_geos(fitbounds="locations", visible=False)
                    .update_layout(height=600, margin={"r": 0, "t": 0, "l": 0, "b": 0},)
                )
            ],
            style={"display": "inline-block", "width": "90%"},
        ),
    ]
)


"""NYISO LAYOUT"""
nyiso_layout = html.Div(
    [
        html.Div(id="nyiso-content"),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Link(
                                html.Button("HOME", id="home-button", className="mr-1"),
                                href=f"/{app_name}/",
                            ),
                            dcc.Link(
                                html.Button(
                                    "CAISO", id="caiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/caiso",
                            ),
                            dcc.Link(
                                html.Button("MISO", id="miso-button", className="mr-1"),
                                href=f"/{app_name}/miso",
                            ),
                            dcc.Link(
                                html.Button("PJM", id="pjm-button", className="mr-1"),
                                href=f"/{app_name}/pjm",
                            ),
                            dcc.Link(
                                html.Button(
                                    "NYISO", id="nyiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/nyiso",
                            ),
                            dcc.Link(
                                html.Button(
                                    "ISONE", id="isone-button", className="mr-1"
                                ),
                                href=f"/{app_name}/isone",
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(width=7),
            ],
            justify="center",
        ),
        html.Br(),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.H1("New York Independent System Operator (NYISO)"), width=9
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
            "The NYISO is the New York Independent System Operator — the organization responsible for
            managing New York’s electric grid and its competitive wholesale electric marketplace." For more information,
            visit https://www.nyiso.com/.
        """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Model Performance"), width=9), dbc.Col(width=2),],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""Mean Absolute Error (MAE) for February, 2021: 347.62 (pretty good)"""
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id="nyiso-dropdown",
                        options=[
                            {"label": "Actual", "value": "Actual"},
                            {"label": "Predicted", "value": "Predicted"},
                        ],
                        value=["Actual", "Predicted"],
                        multi=True,
                    ),
                    width=6,
                ),
                dbc.Col(width=5),
            ],
            justify="center",
        ),
        dcc.Graph(id="nyiso-graph"),
        html.Br(),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Training Data"), width=9), dbc.Col(width=2)],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
                    The NYISO forecasting model was trained on historical load and weather data
                    from 2018-2021. Temperature readings are from New York City.
                """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=px.histogram(
                                    peak_data["NYISO"],
                                    x=peak_data["NYISO"]["load_MW"],
                                    nbins=75,
                                    marginal="rug",
                                    title=f"Distribution of NYISO Daily Peaks",
                                    color_discrete_sequence=["darkturquoise"],
                                ).update_layout(
                                    template=TEMPLATE,
                                    xaxis_title="Peak Load (MW)",
                                    yaxis_title="Number of Days",
                                )
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=go.Figure()
                                .add_trace(
                                    go.Scatter(
                                        x=load_duration_curves["NYISO"]
                                        .reset_index()
                                        .index,
                                        y=load_duration_curves["NYISO"].values,
                                        mode="lines",
                                        fill="tozeroy",
                                        line=dict(color="maroon", width=3),
                                    )
                                )
                                .update_layout(
                                    title="Peak Load Sorted by Day (Highest to Lowest)",
                                    xaxis_title="Number of Days",
                                    yaxis_title="Load (MW)",
                                    template=TEMPLATE,
                                ),
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Dropdown(
                                id="nyiso-scatter-dropdown",
                                options=[
                                    {"label": "Day of Week", "value": "weekday"},
                                    {"label": "Season", "value": "season"},
                                ],
                                value="season",
                                multi=False,
                            ),
                            dcc.Graph(id="nyiso-scatter"),
                        ]
                    ),
                    width=4,
                ),
            ]
        ),
    ]
)


@app.callback(
    dash.dependencies.Output("nyiso-content", "children"),
    [dash.dependencies.Input("nyiso-button", "value")],
)
@app.callback(
    dash.dependencies.Output("nyiso-graph", "figure"),
    [dash.dependencies.Input("nyiso-dropdown", "value")],
)
def plot_nyiso_load_(value):
    fig = go.Figure()
    if "Actual" in value:
        fig.add_trace(
            go.Scatter(
                x=load["NYISO"].index,
                y=load["NYISO"].values,
                name="Actual Load",
                line=dict(color="maroon", width=3),
            )
        )
    if "Predicted" in value:
        fig.add_trace(
            go.Scatter(
                x=predictions["NYISO"].index,
                y=predictions["NYISO"].values,
                name="Forecasted Load",
                line=dict(color="darkturquoise", width=3, dash="dash"),
            )
        )
    return fig.update_layout(
        title="System Load: Actual vs. Predicted",
        xaxis_title="Date",
        yaxis_title="Load (MW)",
        template=TEMPLATE,
    )


@app.callback(
    dash.dependencies.Output("nyiso-scatter", "figure"),
    [dash.dependencies.Input("nyiso-scatter-dropdown", "value")],
)
def nyiso_scatter_plot(value):
    fig = px.scatter(peak_data["NYISO"], x="load_MW", y="temperature", color=value)
    return fig.update_layout(template=TEMPLATE, title="Peak Load vs. Temperature")


app_name = os.getenv("APP_NAME", "dash-peaky-finders")

"""PJM LAYOUT"""
pjm_layout = html.Div(
    [
        html.Div(id="pjm-content"),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Link(
                                html.Button("HOME", id="home-button", className="mr-1"),
                                href=f"/{app_name}/",
                            ),
                            dcc.Link(
                                html.Button(
                                    "CAISO", id="caiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/caiso",
                            ),
                            dcc.Link(
                                html.Button("MISO", id="miso-button", className="mr-1"),
                                href=f"/{app_name}/miso",
                            ),
                            dcc.Link(
                                html.Button("PJM", id="pjm-button", className="mr-1"),
                                href=f"/{app_name}/pjm",
                            ),
                            dcc.Link(
                                html.Button(
                                    "NYISO", id="nyiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/nyiso",
                            ),
                            dcc.Link(
                                html.Button(
                                    "ISONE", id="isone-button", className="mr-1"
                                ),
                                href=f"/{app_name}/isone",
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(width=7),
            ],
            justify="center",
        ),
        html.Br(),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.H1("Pennsylvania, Jersey, Maryland Power Pool (PJM)"), width=9
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
            "PJM is a regional transmission organization (RTO) that coordinates the
            movement of wholesale electricity in all or parts of 13 states and
            the District of Columbia." For more information, visit https://www.pjm.com.
        """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Model Performance"), width=9), dbc.Col(width=2),],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""Mean Absolute Error (MAE) for February, 2021: 2,886.66 (not great)"""
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id="pjm-dropdown",
                        options=[
                            {"label": "Actual", "value": "Actual"},
                            {"label": "Predicted", "value": "Predicted"},
                        ],
                        value=["Actual", "Predicted"],
                        multi=True,
                    ),
                    width=6,
                ),
                dbc.Col(width=5),
            ],
            justify="center",
        ),
        dcc.Graph(id="pjm-graph"),
        html.Br(),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Training Data"), width=9), dbc.Col(width=2)],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
                    The PJM forecasting model was trained on historical load and weather data
                    from 2018-2021. Temperature readings are from Philadelphia.
                """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=px.histogram(
                                    peak_data["PJM"],
                                    x=peak_data["PJM"]["load_MW"],
                                    nbins=75,
                                    marginal="rug",
                                    title=f"Distribution of PJM Daily Peaks",
                                    color_discrete_sequence=["darkturquoise"],
                                ).update_layout(
                                    template=TEMPLATE,
                                    xaxis_title="Peak Load (MW)",
                                    yaxis_title="Number of Days",
                                )
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=go.Figure()
                                .add_trace(
                                    go.Scatter(
                                        x=load_duration_curves["PJM"]
                                        .reset_index()
                                        .index,
                                        y=load_duration_curves["PJM"].values,
                                        mode="lines",
                                        fill="tozeroy",
                                        line=dict(color="maroon", width=3),
                                    )
                                )
                                .update_layout(
                                    title="Peak Load Sorted by Day (Highest to Lowest)",
                                    xaxis_title="Number of Days",
                                    yaxis_title="Load (MW)",
                                    template=TEMPLATE,
                                ),
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Dropdown(
                                id="pjm-scatter-dropdown",
                                options=[
                                    {"label": "Day of Week", "value": "weekday"},
                                    {"label": "Season", "value": "season"},
                                ],
                                value="season",
                                multi=False,
                            ),
                            dcc.Graph(id="pjm-scatter"),
                        ]
                    ),
                    width=4,
                ),
            ]
        ),
    ]
)


@app.callback(
    dash.dependencies.Output("pjm-content", "children"),
    [dash.dependencies.Input("pjm-button", "value")],
)
@app.callback(
    dash.dependencies.Output("pjm-graph", "figure"),
    [dash.dependencies.Input("pjm-dropdown", "value")],
)
def plot_pjm_load_(value):
    fig = go.Figure()
    if "Actual" in value:
        fig.add_trace(
            go.Scatter(
                x=load["PJM"].index,
                y=load["PJM"].values,
                name="Actual Load",
                line=dict(color="maroon", width=3),
            )
        )
    if "Predicted" in value:
        fig.add_trace(
            go.Scatter(
                x=predictions["PJM"].index,
                y=predictions["PJM"].values,
                name="Forecasted Load",
                line=dict(color="darkturquoise", width=3, dash="dash"),
            )
        )
    return fig.update_layout(
        title="System Load: Actual vs. Predicted",
        xaxis_title="Date",
        yaxis_title="Load (MW)",
        template=TEMPLATE,
    )


@app.callback(
    dash.dependencies.Output("pjm-scatter", "figure"),
    [dash.dependencies.Input("pjm-scatter-dropdown", "value")],
)
def pjm_scatter_plot(value):
    fig = px.scatter(peak_data["PJM"], x="load_MW", y="temperature", color=value)
    return fig.update_layout(template=TEMPLATE, title="Peak Load vs. Temperature")


"""MISO LAYOUT"""
miso_layout = html.Div(
    [
        html.Div(id="miso-content"),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Link(
                                html.Button("HOME", id="home-button", className="mr-1"),
                                href=f"/{app_name}/",
                            ),
                            dcc.Link(
                                html.Button(
                                    "CAISO", id="caiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/caiso",
                            ),
                            dcc.Link(
                                html.Button("MISO", id="miso-button", className="mr-1"),
                                href=f"/{app_name}/miso",
                            ),
                            dcc.Link(
                                html.Button("PJM", id="pjm-button", className="mr-1"),
                                href=f"/{app_name}/pjm",
                            ),
                            dcc.Link(
                                html.Button(
                                    "NYISO", id="nyiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/nyiso",
                            ),
                            dcc.Link(
                                html.Button(
                                    "ISONE", id="isone-button", className="mr-1"
                                ),
                                href=f"/{app_name}/isone",
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(width=7),
            ],
            justify="center",
        ),
        html.Br(),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.H1("Midcontinent Independent System Operator (MISO)"), width=9
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
                        "Midcontinent Independent System Operator (MISO) is an independent,
                        not-for-profit organization that delivers safe, cost-effective
                        electric power across 15 U.S. states and the Canadian province of
                        Manitoba." For more information,
                        visit www.misoenergy.org.
                        """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Model Performance"), width=9), dbc.Col(width=2),],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""Mean Absolute Error (MAE) for February, 2021: 2382.66 (not great)"""
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id="miso-dropdown",
                        options=[
                            {"label": "Actual", "value": "Actual"},
                            {"label": "Predicted", "value": "Predicted"},
                        ],
                        value=["Actual", "Predicted"],
                        multi=True,
                    ),
                    width=6,
                ),
                dbc.Col(width=5),
            ],
            justify="center",
        ),
        dcc.Graph(id="miso-graph"),
        html.Br(),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Training Data"), width=9), dbc.Col(width=2)],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
                                    The MISO forecasting model was trained on historical load and weather data
                                    from 2018-2021. Temperature readings are from Minneapolis.
                                """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=px.histogram(
                                    peak_data["MISO"],
                                    x=peak_data["MISO"]["load_MW"],
                                    nbins=75,
                                    marginal="rug",
                                    title=f"Distribution of MISO Daily Peaks",
                                    color_discrete_sequence=["darkturquoise"],
                                ).update_layout(
                                    template=TEMPLATE,
                                    xaxis_title="Peak Load (MW)",
                                    yaxis_title="Number of Days",
                                )
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=go.Figure()
                                .add_trace(
                                    go.Scatter(
                                        x=load_duration_curves["MISO"]
                                        .reset_index()
                                        .index,
                                        y=load_duration_curves["MISO"].values,
                                        mode="lines",
                                        fill="tozeroy",
                                        line=dict(color="maroon", width=3),
                                    )
                                )
                                .update_layout(
                                    title="Peak Load Sorted by Day (Highest to Lowest)",
                                    xaxis_title="Number of Days",
                                    yaxis_title="Load (MW)",
                                    template=TEMPLATE,
                                ),
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Dropdown(
                                id="miso-scatter-dropdown",
                                options=[
                                    {"label": "Day of Week", "value": "weekday"},
                                    {"label": "Season", "value": "season"},
                                ],
                                value="season",
                                multi=False,
                            ),
                            dcc.Graph(id="miso-scatter"),
                        ]
                    ),
                    width=4,
                ),
            ]
        ),
    ]
)


@app.callback(
    dash.dependencies.Output("miso-content", "children"),
    [dash.dependencies.Input("miso-button", "value")],
)
@app.callback(
    dash.dependencies.Output("miso-graph", "figure"),
    [dash.dependencies.Input("miso-dropdown", "value")],
)
def plot_miso_load_(value):
    fig = go.Figure()
    if "Actual" in value:
        fig.add_trace(
            go.Scatter(
                x=load["MISO"].index,
                y=load["MISO"].values,
                name="Actual Load",
                line=dict(color="maroon", width=3),
            )
        )
    if "Predicted" in value:
        fig.add_trace(
            go.Scatter(
                x=predictions["MISO"].index,
                y=predictions["MISO"].values,
                name="Forecasted Load",
                line=dict(color="darkturquoise", width=3, dash="dash"),
            )
        )
    return fig.update_layout(
        title="System Load: Actual vs. Predicted",
        xaxis_title="Date",
        yaxis_title="Load (MW)",
        template=TEMPLATE,
    )


@app.callback(
    dash.dependencies.Output("miso-scatter", "figure"),
    [dash.dependencies.Input("miso-scatter-dropdown", "value")],
)
def miso_scatter_plot(value):
    fig = px.scatter(peak_data["MISO"], x="load_MW", y="temperature", color=value)
    return fig.update_layout(template=TEMPLATE, title="Peak Load vs. Temperature")


"""ISONE LAYOUT"""
isone_layout = html.Div(
    [
        html.Div(id="isone-content"),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Link(
                                html.Button("HOME", id="home-button", className="mr-1"),
                                href=f"/{app_name}/",
                            ),
                            dcc.Link(
                                html.Button(
                                    "CAISO", id="caiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/caiso",
                            ),
                            dcc.Link(
                                html.Button("MISO", id="miso-button", className="mr-1"),
                                href=f"/{app_name}/miso",
                            ),
                            dcc.Link(
                                html.Button("PJM", id="pjm-button", className="mr-1"),
                                href=f"/{app_name}/pjm",
                            ),
                            dcc.Link(
                                html.Button(
                                    "NYISO", id="nyiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/nyiso",
                            ),
                            dcc.Link(
                                html.Button(
                                    "ISONE", id="isone-button", className="mr-1"
                                ),
                                href=f"/{app_name}/isone",
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(width=7),
            ],
            justify="center",
        ),
        html.Br(),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.H1("Independent System Operator of New England (ISONE)"),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
            ISONE is the "independent, not-for-profit corporation responsible
            for keeping electricity flowing across the six New England states
            and ensuring that the region has reliable, competitively priced
            wholesale electricity today and into the future." For more information,
            visit www.iso-ne.com.
        """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Model Performance"), width=9), dbc.Col(width=2),],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""Mean Absolute Error (MAE) for February, 2021: 522.43 (pretty good)"""
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id="isone-dropdown",
                        options=[
                            {"label": "Actual", "value": "Actual"},
                            {"label": "Predicted", "value": "Predicted"},
                        ],
                        value=["Actual", "Predicted"],
                        multi=True,
                    ),
                    width=6,
                ),
                dbc.Col(width=5),
            ],
            justify="center",
        ),
        dcc.Graph(id="isone-graph"),
        html.Br(),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Training Data"), width=9), dbc.Col(width=2)],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
                    The ISONE model was trained on historical load and weather data
                    from 2018-2021. Temperature readings are from Boston.
                """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=px.histogram(
                                    peak_data["ISONE"],
                                    x=peak_data["ISONE"]["load_MW"],
                                    nbins=75,
                                    marginal="rug",
                                    title=f"Distribution of ISONE Daily Peaks",
                                    color_discrete_sequence=["darkturquoise"],
                                ).update_layout(
                                    template=TEMPLATE,
                                    xaxis_title="Peak Load (MW)",
                                    yaxis_title="Number of Days",
                                )
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=go.Figure()
                                .add_trace(
                                    go.Scatter(
                                        x=load_duration_curves["ISONE"]
                                        .reset_index()
                                        .index,
                                        y=load_duration_curves["ISONE"].values,
                                        mode="lines",
                                        fill="tozeroy",
                                        line=dict(color="maroon", width=3),
                                    )
                                )
                                .update_layout(
                                    title="Peak Load Sorted by Day (Highest to Lowest)",
                                    xaxis_title="Number of Days",
                                    yaxis_title="Load (MW)",
                                    template=TEMPLATE,
                                ),
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Dropdown(
                                id="isone-scatter-dropdown",
                                options=[
                                    {"label": "Day of Week", "value": "weekday"},
                                    {"label": "Season", "value": "season"},
                                ],
                                value="season",
                                multi=False,
                            ),
                            dcc.Graph(id="isone-scatter"),
                        ]
                    ),
                    width=4,
                ),
            ]
        ),
    ]
)


@app.callback(
    dash.dependencies.Output("isone-content", "children"),
    [dash.dependencies.Input("isone-button", "value")],
)
@app.callback(
    dash.dependencies.Output("isone-graph", "figure"),
    [dash.dependencies.Input("isone-dropdown", "value")],
)
def plot_isone_load_(value):
    fig = go.Figure()
    if "Actual" in value:
        fig.add_trace(
            go.Scatter(
                x=load["ISONE"].index,
                y=load["ISONE"].values,
                name="Actual Load",
                line=dict(color="maroon", width=3),
            )
        )
    if "Predicted" in value:
        fig.add_trace(
            go.Scatter(
                x=predictions["ISONE"].index,
                y=predictions["ISONE"].values,
                name="Forecasted Load",
                line=dict(color="darkturquoise", width=3, dash="dash"),
            )
        )
    return fig.update_layout(
        title="System Load: Actual vs. Predicted",
        xaxis_title="Date",
        yaxis_title="Load (MW)",
        template=TEMPLATE,
    )


@app.callback(
    dash.dependencies.Output("isone-scatter", "figure"),
    [dash.dependencies.Input("isone-scatter-dropdown", "value")],
)
def isone_scatter_plot(value):
    fig = px.scatter(peak_data["ISONE"], x="load_MW", y="temperature", color=value)
    return fig.update_layout(template=TEMPLATE, title="Peak Load vs. Temperature")


"""CAISO LAYOUT"""
caiso_layout = html.Div(
    [
        html.Div(id="caiso-content"),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Link(
                                html.Button("HOME", id="home-button", className="mr-1"),
                                href=f"/{app_name}/",
                            ),
                            dcc.Link(
                                html.Button(
                                    "CAISO", id="caiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/caiso",
                            ),
                            dcc.Link(
                                html.Button("MISO", id="miso-button", className="mr-1"),
                                href=f"/{app_name}/miso",
                            ),
                            dcc.Link(
                                html.Button("PJM", id="pjm-button", className="mr-1"),
                                href=f"/{app_name}/pjm",
                            ),
                            dcc.Link(
                                html.Button(
                                    "NYISO", id="nyiso-button", className="mr-1"
                                ),
                                href=f"/{app_name}/nyiso",
                            ),
                            dcc.Link(
                                html.Button(
                                    "ISONE", id="isone-button", className="mr-1"
                                ),
                                href=f"/{app_name}/isone",
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(width=7),
            ],
            justify="center",
        ),
        html.Br(),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.H1("California Independent System Operator (CAISO)"), width=9
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
            "The California Independent System Operator (ISO) maintains
            reliability on one of the largest and most modern power grids in
            the world, and operates a transparent, accessible wholesale energy
            market."  For more information,
            visit http://www.caiso.com/.
        """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Model Performance"), width=9), dbc.Col(width=2),],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""Mean Absolute Error (MAE) for February, 2021: 455.91 (pretty good)"""
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id="caiso-dropdown",
                        options=[
                            {"label": "Actual", "value": "Actual"},
                            {"label": "Predicted", "value": "Predicted"},
                        ],
                        value=["Actual", "Predicted"],
                        multi=True,
                    ),
                    width=6,
                ),
                dbc.Col(width=5),
            ],
            justify="center",
        ),
        dcc.Graph(id="caiso-graph"),
        html.Br(),
        html.Br(),
        dbc.Row(
            [dbc.Col(html.H3("Training Data"), width=9), dbc.Col(width=2)],
            justify="center",
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        children="""
                    The CAISO forecasting model was trained on historical load and weather data
                    from 2018-2021. Temperature readings were from Los Angeles.
                """
                    ),
                    width=9,
                ),
                dbc.Col(width=2),
            ],
            justify="center",
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=px.histogram(
                                    peak_data["CAISO"],
                                    x=peak_data["CAISO"]["load_MW"],
                                    nbins=75,
                                    marginal="rug",
                                    title=f"Distribution of CAISO Daily Peaks",
                                    color_discrete_sequence=["darkturquoise"],
                                ).update_layout(
                                    template=TEMPLATE,
                                    xaxis_title="Peak Load (MW)",
                                    yaxis_title="Number of Days",
                                )
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Graph(
                                figure=go.Figure()
                                .add_trace(
                                    go.Scatter(
                                        x=load_duration_curves["CAISO"]
                                        .reset_index()
                                        .index,
                                        y=load_duration_curves["CAISO"].values,
                                        mode="lines",
                                        fill="tozeroy",
                                        line=dict(color="maroon", width=3),
                                    )
                                )
                                .update_layout(
                                    title="Peak Load Sorted by Day (Highest to Lowest)",
                                    xaxis_title="Number of Days",
                                    yaxis_title="Load (MW)",
                                    template=TEMPLATE,
                                ),
                            ),
                        ]
                    ),
                    width=4,
                ),
                dbc.Col(
                    html.Div(
                        [
                            dcc.Dropdown(
                                id="caiso-scatter-dropdown",
                                options=[
                                    {"label": "Day of Week", "value": "weekday"},
                                    {"label": "Season", "value": "season"},
                                ],
                                value="season",
                                multi=False,
                            ),
                            dcc.Graph(id="caiso-scatter"),
                        ]
                    ),
                    width=4,
                ),
            ]
        ),
    ]
)


@app.callback(
    dash.dependencies.Output("caiso-content", "children"),
    [dash.dependencies.Input("caiso-button", "value")],
)
@app.callback(
    dash.dependencies.Output("caiso-graph", "figure"),
    [dash.dependencies.Input("caiso-dropdown", "value")],
)
def plot_caiso_load_(value):
    fig = go.Figure()
    if "Actual" in value:
        fig.add_trace(
            go.Scatter(
                x=load["CAISO"].index,
                y=load["CAISO"].values,
                name="Actual Load",
                line=dict(color="maroon", width=3),
            )
        )
    if "Predicted" in value:
        fig.add_trace(
            go.Scatter(
                x=predictions["CAISO"].index,
                y=predictions["CAISO"].values,
                name="Forecasted Load",
                line=dict(color="darkturquoise", width=3, dash="dash"),
            )
        )
    return fig.update_layout(
        title="System Load: Actual vs. Predicted",
        xaxis_title="Date",
        yaxis_title="Load (MW)",
        template=TEMPLATE,
    )


@app.callback(
    dash.dependencies.Output("caiso-scatter", "figure"),
    [dash.dependencies.Input("caiso-scatter-dropdown", "value")],
)
def caiso_scatter_plot(value):
    fig = px.scatter(
        peak_data["CAISO"].dropna(), x="load_MW", y="temperature", color=value
    )
    return fig.update_layout(template=TEMPLATE, title="Peak Load vs. Temperature")


# Update the index
@app.callback(
    dash.dependencies.Output("page-content", "children"),
    [dash.dependencies.Input("url", "pathname")],
)
def display_page(pathname):
    if pathname.endswith("/nyiso"):
        return nyiso_layout
    elif pathname.endswith("/pjm"):
        return pjm_layout
    elif pathname.endswith("/isone"):
        return isone_layout
    elif pathname.endswith("/miso"):
        return miso_layout
    elif pathname.endswith("/caiso"):
        return caiso_layout
    else:
        return index_page


if __name__ == "__main__":
    app.run_server(debug=False)
