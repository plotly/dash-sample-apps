import os
import time

import cudf
import cuml
import cupy as cp
import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.express as px


# Load CSV into a cudf
data_dir = os.environ.get("DATA_DIR", "")
path = os.path.join(data_dir, "creditcard.csv")
gdf = cudf.read_csv(path)
gdf.Time = gdf.Time / 3600
gdf.loc[gdf.Amount > 500, "Amount"] = 500

# Define app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server


controls = dbc.Row(
    [
        dbc.Col(
            dbc.FormGroup(
                [
                    dbc.Label("Time since start (h)"),
                    dcc.RangeSlider(
                        id="slider-hours",
                        min=0,
                        max=50,
                        step=1,
                        value=[20, 30],
                        marks={i: str(i) for i in range(0, 51, 10)},
                    ),
                ]
            ),
            md=6,
        ),
        dbc.Col(
            dbc.FormGroup(
                [
                    dbc.Label("Transaction Amount ($)"),
                    dcc.RangeSlider(
                        id="slider-amount",
                        min=0,
                        max=500,
                        step=5,
                        value=[200, 300],
                        marks={i: str(i) for i in range(0, 501, 100)},
                    ),
                ]
            ),
            md=6,
        ),
    ],
)


# Define Layout
app.layout = dbc.Container(
    fluid=True,
    children=[
        html.H1("Dash cuML UMAP Demo"),
        html.Hr(),
        dbc.Card(controls, body=True),
        dcc.Graph(id="graph-umap", style={"height": "70vh", "max-height": "90vw"}),
        html.Div(id="output-info"),
    ],
    style={"max-width": "960px", "margin": "auto"},
)


@app.callback(
    [Output("graph-umap", "figure"), Output("output-info", "children")],
    [Input("slider-amount", "value"), Input("slider-hours", "value"),],
)
def update_graph(amt, hrs):
    t0 = time.time()
    # First, filter based on the slider values
    time_mask = (gdf.Time >= hrs[0]) & (gdf.Time <= hrs[1])
    amount_mask = (gdf.Amount >= amt[0]) & (gdf.Amount <= amt[1])
    filt_df = gdf.loc[time_mask & amount_mask]

    # Then, select the features and train a UMAP model with cuML
    features = filt_df.loc[:, "V1":"V28"].values
    reducer = cuml.UMAP()
    embedding = reducer.fit_transform(features)

    # Convert the embedding back to numpy
    embedding = cp.asnumpy(embedding)
    amount = cp.asnumpy(filt_df.Amount.values.round(2))

    # Create a plotly.express scatter plot
    fig = px.scatter(
        x=embedding[:, 0],
        y=embedding[:, 1],
        color=amount,
        labels={"color": "Amount ($)"},
        title="UMAP projection of credit card transactions",
    )

    t1 = time.time()
    out_msg = f"Projected {embedding.shape[0]} transactions in {t1-t0:.2f}s."
    alert = dbc.Alert(out_msg, color="success", dismissable=True)

    return fig, alert


if __name__ == "__main__":
    app.run_server(debug=True)
