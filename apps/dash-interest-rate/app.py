import os
import time
from textwrap import dedent

import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from flask_caching import Cache
import plotly.express as px
import pandas as pd
from snowflake.sqlalchemy import URL
from sqlalchemy import create_engine
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split
from sklearn import metrics

from utils import *

os.environ["REDIS_URL"] = os.getenv("REDIS_URL", os.getenv("EXTERNAL_REDIS_URL"))
# os.environ['DATABASE_URL'] = os.getenv('DATABASE_URL', os.getenv('EXTERNAL_DATABASE_URL'))

# Snowflake vars
FLAKE_ACCOUNT = os.getenv("FLAKE_ACCOUNT")
FLAKE_USER = os.getenv("FLAKE_USER")
FLAKE_PW = os.getenv("FLAKE_PW")

flake_warehouse = "snowflake_demos"
flake_db = "LOANS"
stage = "LOAN_STAGE"


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

connection = engine.connect()


# Get distincts and range
loan_min, loan_max = get_range(connection, flake_db, "LOAN_AMNT")
inc_min, inc_max = get_range(connection, flake_db, "ANNUAL_INC")
app_types = get_unique(connection, flake_db, "APPLICATION_TYPE")
purposes = get_unique(connection, flake_db, "PURPOSE")
ownerships = get_unique(connection, flake_db, "HOME_OWNERSHIP")[1:]

# Close connection
connection.close()

# Make some calculations based on value range retrieved
loan_marks = loan_max // 4
loan_min //= loan_marks


# Define app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP],)
app.title = "Interest Rate Modeling"
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


# Build component parts
avp_graph = dcc.Graph(id="avp-graph", style={"height": "500px"})
div_alert = dbc.Spinner(html.Div(id="alert-msg"))
query_card = dbc.Card(
    [
        html.H4("Auto-generated SnowSQL Query", className="card-title"),
        dcc.Markdown(id="sql-query"),
    ],
    body=True,
)

controls = [
    CustomRangeSlider(
        id="loan-amount",
        label="Loan Amount($)",
        values=range(loan_min, loan_max + 1, loan_marks),
    ),
    CustomRangeSlider(
        id="annual-income",
        label="Annual Income ($)",
        values=[0, 20000, 50000, 100000, 200000],
    ),
    OptionMenu(id="app-type", label="Application Type", values=app_types),
    OptionMenu(id="home-ownership", label="Home Ownership", values=ownerships),
    OptionMenu(id="purpose", label="Purpose", values=purposes),
    dbc.Button("Query and Train Ridge", color="primary", id="button-train"),
]


# Define Layout
app.layout = dbc.Container(
    fluid=True,
    children=[
        html.H1("Dash Interest Rate Modeling with Snowflake"),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col([dbc.Card(controls, body=True), div_alert], md=3),
                dbc.Col([avp_graph, query_card], md=4),
                dbc.Col(dcc.Graph(id="coef-graph", style={"height": "800px"}), md=5),
            ]
        ),
    ],
    style={"margin": "auto"},
)


@app.callback(
    [
        Output("alert-msg", "children"),
        Output("avp-graph", "figure"),
        Output("coef-graph", "figure"),
        Output("sql-query", "children"),
    ],
    [Input("button-train", "n_clicks")],
    [
        State("loan-amount", "value"),
        State("annual-income", "value"),
        State("app-type", "value"),
        State("home-ownership", "value"),
        State("purpose", "value"),
    ],
)
def query_and_train(n_clicks, loan_range, inc_range, app_type, owner, purpose):
    t0 = time.time()
    query = dedent(
        f"""
    SELECT *
    FROM {flake_db}.PUBLIC.LOAN_CLEAN
    WHERE
        LOAN_AMNT BETWEEN {loan_range[0]} AND {loan_range[1]} AND
        ANNUAL_INC BETWEEN {inc_range[0]} AND {inc_range[1]} AND
        HOME_OWNERSHIP = '{owner}' AND
        APPLICATION_TYPE = '{app_type}' AND
        PURPOSE = '{purpose}';
    """
    )
    df = connect_read_sql(query=query, engine=engine)

    # Model Training
    df = df.drop(columns=["grade", "policy_code"])
    cat_cols = df.columns[df.dtypes == "object"]
    num_cols = df.columns[df.dtypes != "object"]
    enc_df = pd.get_dummies(df, columns=cat_cols, prefix_sep=": ")
    X = enc_df.drop(columns="int_rate")
    y = enc_df["int_rate"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=0
    )

    model = Ridge(normalize=True, alpha=0.02)
    model.fit(X_train, y_train)

    # Score model
    y_pred = model.predict(X_test)
    test_mae = metrics.mean_absolute_error(y_test, y_pred)

    # Coefficient Figure
    coef_fig = px.bar(
        y=X_test.columns,
        x=model.coef_,
        orientation="h",
        color=X_test.columns.isin(num_cols),
        labels={"color": "Is numerical", "x": "Weight on Prediction", "y": "Features"},
        title="Ridge feature importance for predicting interest rate",
    )

    # Actual vs Predicted Figure
    avp_fig = px.scatter(
        x=y_test,
        y=y_pred,
        labels={"x": "Actual", "y": "Predicted"},
        title=f"Actual vs predicted interest rate (MAE={test_mae:.2f})",
    )

    avp_fig.add_shape(
        type="line", x0=y_test.min(), y0=y_test.min(), x1=y_test.max(), y1=y_test.max()
    )

    t1 = time.time()
    exec_time = t1 - t0
    query_size = df.shape[0]
    alert_msg = f"Queried {query_size} records. Total time: {exec_time:.2f}s."
    alert = dbc.Alert(alert_msg, color="success", dismissable=True)

    return alert, avp_fig, coef_fig, f"```\n{query}\n```"


if __name__ == "__main__":
    app.run_server(debug=True)
