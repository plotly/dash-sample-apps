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
from dash_table import DataTable
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

os.environ["REDIS_URL"] = os.getenv("REDIS_URL", os.getenv("EXTERNAL_REDIS_URL"))
# os.environ['DATABASE_URL'] = os.getenv('DATABASE_URL', os.getenv('EXTERNAL_DATABASE_URL'))


def connect_read_sql(query, engine):
    connection = engine.connect()
    result = pd.read_sql(query, connection)
    connection.close()
    return result


def svg_to_fig(svg_bytes, title=None, plot_bgcolor="white", x_lock=False, y_lock=True):
    svg_enc = base64.b64encode(svg_bytes)
    svg = f"data:image/svg+xml;base64, {svg_enc.decode()}"

    # Get the width and height
    xml_tree = xml.etree.ElementTree.fromstring(svg_bytes.decode())
    img_width = int(xml_tree.attrib["width"].strip("pt"))
    img_height = int(xml_tree.attrib["height"].strip("pt"))

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
            opacity=1,
            layer="below",
        )
    )

    # Adapt axes to the right width and height, lock aspect ratio
    fig.update_xaxes(showgrid=False, visible=False, range=[0, img_width])
    fig.update_yaxes(showgrid=False, visible=False, range=[img_height, 0])

    if x_lock is True:
        fig.update_xaxes(constrain="domain")
    if y_lock is True:
        fig.update_yaxes(scaleanchor="x", scaleratio=1)

    fig.update_layout(plot_bgcolor=plot_bgcolor, margin=dict(r=5, l=5, b=5))

    if title:
        fig.update_layout(title=title)

    return fig


# Models
MODELS = ("model-shallow", "model-deep", "model-with-entropy", "model-random-split")

# Snowflake vars
FLAKE_ACCOUNT = os.getenv("FLAKE_ACCOUNT")
FLAKE_USER = os.getenv("FLAKE_USER")
FLAKE_PW = os.getenv("FLAKE_PW")

flake_warehouse = "snowflake_demos"
flake_db = "NEW_CLIENTS"
stage = "NEW_CLIENT_STAGE"


# Load feature names
feature_names = pickle.load(open("feature_names.pickle", "rb"))


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
app.title = "Loan Grade Classification with Snowflake"
server = app.server


query_card = dbc.Card(
    [
        dbc.CardHeader("Dynamic SnowSQL Query"),
        dbc.CardBody(dcc.Markdown(id="sql-query")),
    ]
)

output_card = dbc.Card(
    [
        dbc.CardHeader("Real-time Prediction"),
        dbc.CardBody(html.H2(id="predicted-grade", style={"text-align": "center"})),
    ]
)

model_selection = dbc.InputGroup(
    [
        dbc.InputGroupAddon("Select Model", addon_type="prepend"),
        dbc.Select(
            id="model-selection",
            options=[
                {
                    "label": m.replace("-", " ").capitalize(),
                    "value": "assets/" + m + ".joblib",
                }
                for m in MODELS
            ],
            value="assets/" + MODELS[0] + ".joblib",
        ),
    ]
)


table_selection = dbc.InputGroup(
    [
        dbc.InputGroupAddon("Query Table", addon_type="prepend"),
        dbc.Select(
            id="table-selection",
            options=[{"label": t, "value": t} for t in ["RENT", "MORTGAGE", "OWNER"]],
            value="RENT",
        ),
    ]
)

sample_controls = [
    dbc.Col(table_selection),
    dbc.Col(
        dbc.ButtonGroup(
            [
                dbc.Button(
                    "Prev. Sample",
                    id="prev-sample",
                    color="success",
                    n_clicks=0,
                    outline=True,
                ),
                dbc.Button(
                    "Next Sample", id="next-sample", n_clicks=0, color="success"
                ),
            ],
            style={"width": "100%"},
        )
    ),
]

controls = [dbc.Col(dbc.Row(sample_controls), md=7), dbc.Col(model_selection, md=5)]


# Define Layout
app.layout = dbc.Container(
    fluid=True,
    children=[
        dcc.Store(id="query-store"),
        html.H1("Dash Loan Grade Classification with Snowflake"),
        html.Hr(),
        dbc.Row(controls, style={"padding": "20px 0px"}),
        dbc.Row(
            [
                dbc.Col(
                    children=[
                        dbc.CardDeck([query_card, output_card]),
                        # dbc.Row(
                        #     [dbc.Col(query_card), dbc.Col(output_card)]
                        # ),
                        DataTable(
                            id="table-sample",
                            style_table={
                                "height": "500px",
                                "overflowY": "auto",
                                "padding": "20px",
                            },
                        ),
                    ],
                    md=7,
                ),
                dbc.Col(dcc.Graph(id="graph-tree", style={"height": "700px"}), md=5),
            ]
        ),
    ],
    style={"margin": "auto"},
)


@app.callback(Output("graph-tree", "figure"), [Input("model-selection", "value")])
def visualize_tree(path):
    model = joblib.load(open(path, "rb"))
    dot_data = export_graphviz(
        model,
        out_file=None,
        filled=True,
        rounded=True,
        feature_names=feature_names,
        class_names=model.classes_,
        proportion=True,
        rotate=True,
        precision=2,
    )

    pydot_graph = pydot.graph_from_dot_data(dot_data)[0]
    svg_bytes = pydot_graph.create_svg()
    fig = svg_to_fig(svg_bytes, title="Decision Tree Explanation")

    return fig


@app.callback(
    [Output("query-store", "data"), Output("sql-query", "children")],
    [Input("table-selection", "value")],
)
def query_and_store(table_name):
    query = f"""
SELECT *
FROM NEW_CLIENTS.PUBLIC.{table_name}
LIMIT 100;
    """
    query_df = connect_read_sql(query, engine)

    return query_df.to_json(), f"```\n{query}\n```"


@app.callback(
    [
        Output("table-sample", "data"),
        Output("table-sample", "columns"),
        Output("predicted-grade", "children"),
    ],
    [
        Input("query-store", "data"),
        Input("prev-sample", "n_clicks"),
        Input("next-sample", "n_clicks"),
        Input("model-selection", "value"),
    ],
)
def generate_table(query_json, prev_clicks, next_clicks, model_path):

    query_df = pd.read_json(query_json)
    cat_cols = query_df.columns[query_df.dtypes == "object"]
    dummies_df = pd.get_dummies(query_df, columns=cat_cols, prefix_sep=": ")

    # Create the missing columns that the model needs to process input
    missing_cols = set(feature_names) - set(query_df.columns)
    for col in missing_cols:
        dummies_df[col] = 0

    # Build the sample table
    i = max(0, next_clicks - prev_clicks)
    table_df = query_df.loc[i:i].T.reset_index()
    table_df.columns = ["Feature", "Current Value"]

    # Load model and make prediction
    model = joblib.load(open(model_path, "rb"))
    model_input = dummies_df.loc[i:i]
    pred = model.predict(model_input)[0]

    columns = [{"name": i, "id": i} for i in table_df.columns]

    return table_df.to_dict("records"), columns, f"Grade = {pred}"


if __name__ == "__main__":
    app.run_server(debug=True)
