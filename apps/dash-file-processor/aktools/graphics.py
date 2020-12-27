import dash
from dash.dependencies import Output, Input
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from aktools.conf import hosts,dashboard_refresh


external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

colors = {"background": "#111111", "text": "#7FDBFF"}


app.layout = html.Div(
    style={"backgroundColor": colors["background"]},
    children=[
        html.H1(
            children="Network log analytics",
            style={"textAlign": "center", "color": colors["text"]},
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.H2(
                            """Select a host:""",
                            style={"margin-right": "1em", "color": colors["text"]},
                        )
                    ],
                ),
                dcc.Dropdown(
                    id="hosts_dropdown",
                    options=[{"label": i, "value": i} for i in hosts],
                    value="Hannibal",  # default value
                    style=dict(width="40%", display="inline-block"),
                ),
            ],
            style={"display": "flex", "align-items": "center"},
        ),
        dcc.Graph(id="live-graphs_host"),
        dcc.Interval(id="graph-update", interval=dashboard_refresh), # graph updates in ms
    ],
)
