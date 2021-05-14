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
from sqlalchemy import create_engine
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split
from sklearn import metrics

import utils

os.environ["REDIS_URL"] = os.getenv("REDIS_URL", os.getenv("EXTERNAL_REDIS_URL"))
# os.environ['DATABASE_URL'] = os.getenv('DATABASE_URL', os.getenv('EXTERNAL_DATABASE_URL'))

SYNAPSE_HOST = os.environ["SYNAPSE_HOST"]
SYNAPSE_USER = os.environ["SYNAPSE_USER"]
SYNAPSE_PW = os.environ["SYNAPSE_PW"]

driver = "ODBC+Driver+17+for+SQL+Server"
db_name = "loan"
table_name = "cleanLoan"
port = "1433"

engine = create_engine(
    f"mssql+pyodbc://{SYNAPSE_USER}:{SYNAPSE_PW}@{SYNAPSE_HOST}:{port}/{db_name}?driver={driver}"
)

connection = engine.connect()


# Get distincts and range
loan_min, loan_max = utils.get_range(connection, db_name, table_name, "loan_amnt")
inc_min, inc_max = utils.get_range(connection, db_name, table_name, "annual_inc")
app_types = utils.get_unique(connection, db_name, table_name, "application_type")
purposes = utils.get_unique(connection, db_name, table_name, "purpose")
ownerships = utils.get_unique(connection, db_name, table_name, "home_ownership")[1:]

# Close connection
connection.close()

# Make some calculations based on value range retrieved
loan_marks = loan_max // 4
loan_min //= loan_marks


# Define app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUMEN])
app.title = "Synapse Analytics in Dash"
server = app.server
cache = Cache(
    app.server,
    config={"CACHE_TYPE": "redis", "CACHE_REDIS_URL": os.environ["REDIS_URL"]},
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
        html.H4("Auto-generated Synapse SQL Query", className="card-title"),
        dcc.Markdown(id="sql-query"),
    ],
    body=True,
)

navbar = dbc.Navbar(
    [
        html.A(
            # Use row and col to control vertical alignment of logo / brand
            dbc.Row(
                [
                    dbc.Col(
                        html.Img(src=app.get_asset_url("dash-logo.png"), height="40px")
                    ),
                    dbc.Col(dbc.NavbarBrand("Synapse Analytics Integration")),
                ],
                align="center",
                no_gutters=True,
            ),
            href="https://plotly.com/dash",
        ),
    ],
)

controls = [
    utils.OptionMenu(id="app-type", label="Application Type", values=app_types),
    utils.OptionMenu(
        id="home-ownership", label="Home Ownership", values=ownerships, value="OWN"
    ),
    utils.OptionMenu(
        id="purpose", label="Purpose", values=purposes, value="credit_card"
    ),
    utils.CustomRangeSlider(
        id="loan-amount",
        label="Loan Amount($)",
        values=range(loan_min, loan_max + 1, loan_marks),
    ),
    utils.CustomRangeSlider(
        id="annual-income",
        label="Annual Income ($)",
        values=[0, 20000, 50000, 100000, 200000],
    ),
    dbc.Button("Query and Train Ridge", color="primary", id="button-train"),
    html.Br(),
]


# Define Layout
app.layout = dbc.Container(
    fluid=True,
    children=[
        navbar,
        html.Br(),
        dbc.Row(
            [
                dbc.Col([avp_graph, query_card], md=4),
                dbc.Col(dcc.Graph(id="coef-graph", style={"height": "800px"}), md=5),
                dbc.Col([dbc.Card(controls, body=True), div_alert], md=3),
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
    SELECT TOP 2500 *
    FROM {db_name}.dbo.{table_name}
    WHERE
        loan_amnt BETWEEN {loan_range[0]} AND {loan_range[1]} AND
        annual_inc BETWEEN {inc_range[0]} AND {inc_range[1]} AND
        home_ownership = '{owner}' AND
        application_type = '{app_type}' AND
        purpose = '{purpose}';
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
