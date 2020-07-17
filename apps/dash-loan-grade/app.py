import base64
import os
import pickle
import time
from textwrap import dedent
import xml

import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from flask_caching import Cache
import joblib
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import pydot
from snowflake.sqlalchemy import URL
from sqlalchemy import create_engine
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier, plot_tree, export_graphviz
from sklearn import metrics


def svg_to_fig(svg_bytes, title=None, plot_bgcolor='white'):
    svg_enc = base64.b64encode(svg_bytes)
    svg = f'data:image/svg+xml;base64, {svg_enc.decode()}'
    
    # Get the width and height
    xml_tree = xml.etree.ElementTree.fromstring(svg_bytes.decode())
    img_width = int(xml_tree.attrib['width'].strip('pt'))
    img_height = int(xml_tree.attrib['height'].strip('pt'))

    fig = go.Figure()
    # Add invisible scatter trace.
    # This trace is added to help the autoresize logic work.
    fig.add_trace(
        go.Scatter(
            x=[0, img_width],
            y=[img_height, 0],
            mode="markers",
            marker_opacity=0,
            hoverinfo="none",
        )
    )
    fig.add_layout_image(
        dict(
            source=svg,
            x=0,
            y=0,
            xref="x",
            yref="y",
            sizex=img_width,
            sizey=img_height,
            sizing="stretch",
            opacity=1,
            layer="below",
        )
    )

    # Adapt axes to the right width and height, lock aspect ratio
    fig.update_xaxes(
        showgrid=False, 
        visible=False, 
        constrain="domain", 
        range=[0, img_width]
    )
    fig.update_yaxes(
        showgrid=False,
        visible=False,
        scaleanchor="x",
        scaleratio=1,
        range=[img_height, 0],
    )
    
    fig.update_layout(plot_bgcolor=plot_bgcolor)

    if title:
        fig.update_layout(title=title)

    return fig


# Models
MODELS = (
    'model-shallow',
    'model-deep',
    'model-with-entropy',
    'model-random-split'
)

# Snowflake vars
FLAKE_ACCOUNT = os.getenv("FLAKE_ACCOUNT")
FLAKE_USER = os.getenv("FLAKE_USER")
FLAKE_PW = os.getenv("FLAKE_PW")

flake_warehouse = "snowflake_demos"
flake_db = "LOANS"
stage = "LOAN_STAGE"


# Load feature names
feature_names = pickle.load(open('feature_names.pickle', 'rb'))


# Create Engine and connect to DB
engine = create_engine(
    URL(
        account=FLAKE_ACCOUNT,
        user=FLAKE_USER,
        password=FLAKE_PW,
        database=flake_db,
        schema="public",
        warehouse=flake_warehouse,
        role="sysadmin",
    ),
    pool_size=5,
    pool_recycle=1800,
    pool_pre_ping=True,
)

# connection = engine.connect()


# Define app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server
cache = Cache(
    app.server,
    config={"CACHE_TYPE": "redis", "CACHE_REDIS_URL": os.environ.get("REDIS_URL", "")},
)


# Cache functions
@cache.memoize(timeout=300)
def connect_read_sql(query, engine):
    connection = engine.connect()
    result = pd.read_sql(query, connection)
    connection.close()
    return result



model_selection = dbc.InputGroup(
    [
        dbc.InputGroupAddon("Select Model", addon_type="prepend"),
        dbc.Select(
            id='model-selection',
            options=[
                {
                    'label': m.replace("-", " ").capitalize(), 
                    'value': 'assets/' + m + '.joblib'
                }
                for m in MODELS
            ],
            value='assets/' + MODELS[0] + '.joblib'
        )
    ]
)

controls = [
    dbc.Col(model_selection)
]


# Define Layout
app.layout = dbc.Container(
    fluid=True,
    children=[
        html.H1("Dash Loan Grade Classification with Snowflake"),
        html.Hr(),
        dbc.Card(dbc.Row(controls), body=True),
        dbc.Row(
            [
                dbc.Col(dcc.Graph(id='graph-tree', style={'height': '700px'}))
            ]
        )
    ],
    style={"margin": "auto"},
)


@app.callback(
    Output('graph-tree', 'figure'), 
    [Input('model-selection', 'value')]
)
def visualize_tree(path):
    print(path)
    model = joblib.load(open(path, 'rb'))
    dot_data = export_graphviz(
        model, 
        out_file=None, 
        filled=True, 
        rounded=True, 
        feature_names=feature_names,
        class_names=model.classes_,
        proportion=True,
        rotate=True,
        precision=2
    )

    pydot_graph = pydot.graph_from_dot_data(dot_data)[0]
    svg_bytes = pydot_graph.create_svg()
    fig = svg_to_fig(svg_bytes, title='Decision Tree Explanation')

    return fig



if __name__ == "__main__":
    app.run_server(debug=True)
