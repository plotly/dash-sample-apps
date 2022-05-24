from dash import Dash, dcc, Input, Output, State
import dash_bootstrap_components as dbc
import numpy as np

from utils.components import Header, controls_first, controls_second, controls_third
from utils.model import svm_prediction


app = Dash(
    __name__,
    title = "Support Vector Machine",
    external_stylesheets=[dbc.themes.CYBORG],
)
server = app.server

app.layout = dbc.Container(
    children=[
        Header(app), 
        dbc.Row([
            dbc.Col([
                dbc.Card(controls_first),
                dbc.Card(controls_second),
                dbc.Card(controls_third)
            ], xl=2, lg=3, sm=4, xs=12, className="control-pannel"),
            dbc.Col(
                dcc.Loading(dcc.Graph(id="graph-sklearn-svm")), 
                xl=7, lg=6, sm=4, xs=12
            ),
            dbc.Col(
                dbc.Row([
                    dcc.Loading(dcc.Graph(id="graph-line-roc-curve", style={'height': '40vh'})),
                    dcc.Loading(dcc.Graph(id="graph-pie-confusion-matrix", style={'height': '40vh'})),
                ]), xl=3, lg=3, sm=4, xs=12
            ),
        ], className="app-body")
    ],
    fluid=True
)


@app.callback(
    Output("slider-svm-parameter-gamma-coef", "marks"),
    Input("slider-svm-parameter-gamma-power", "value"),
)
def update_slider_svm_parameter_gamma_coef(power):
    scale = 10 ** power
    if power < 1:
        return {i: str(round(i * scale, 5)) for i in range(1, 10, 4)}
    else:
        return {i: str(int(i * scale)) for i in range(1, 10, 4)}


@app.callback(
    Output("slider-svm-parameter-C-coef", "marks"),
    Input("slider-svm-parameter-C-power", "value"),
)
def update_slider_svm_parameter_C_coef(power):
    scale = 10 ** power
    if power < 1:
        return {i: str(round(i * scale, 2)) for i in range(1, 10, 4)}
    else:
        return {i: str(int(i * scale)) for i in range(1, 10, 4)}


@app.callback(
    Output("slider-threshold", "value"),
    Input("button-zero-threshold", "n_clicks"),
    State("graph-sklearn-svm", "figure"),
)
def reset_threshold_center(n_clicks, figure):
    if n_clicks:
        Z = np.array(figure["data"][0]["z"])
        value = -Z.min() / (Z.max() - Z.min())
    else:
        value = 0.4
    return value


@app.callback(
    Output("slider-svm-parameter-degree", "disabled"),
    Input("dropdown-svm-parameter-kernel", "value"),
)
def disable_slider_param_degree(kernel):
    return kernel != "poly"


@app.callback(
    Output("slider-svm-parameter-gamma-coef", "disabled"),
    Input("dropdown-svm-parameter-kernel", "value"),
)
def disable_slider_param_gamma_coef(kernel):
    return kernel not in ["rbf", "poly", "sigmoid"]


@app.callback(
    Output("slider-svm-parameter-gamma-power", "disabled"),
    Input("dropdown-svm-parameter-kernel", "value"),
)
def disable_slider_param_gamma_power(kernel):
    return kernel not in ["rbf", "poly", "sigmoid"]


@app.callback(
    Output("graph-sklearn-svm", "figure"),
    Output("graph-line-roc-curve", "figure"),
    Output("graph-pie-confusion-matrix", "figure"),
    Input("dropdown-svm-parameter-kernel", "value"),
    Input("slider-svm-parameter-degree", "value"),
    Input("slider-svm-parameter-C-coef", "value"),
    Input("slider-svm-parameter-C-power", "value"),
    Input("slider-svm-parameter-gamma-coef", "value"),
    Input("slider-svm-parameter-gamma-power", "value"),
    Input("dropdown-select-dataset", "value"),
    Input("slider-dataset-noise-level", "value"),
    Input("radio-svm-parameter-shrinking", "value"),
    Input("slider-threshold", "value"),
    Input("slider-dataset-sample-size", "value"),
)
def update_svm_graph(kernel, degree, C_coef, C_power, gamma_coef, gamma_power, dataset, noise, shrinking, threshold, sample_size):
    return svm_prediction(kernel, degree, C_coef, C_power, gamma_coef, gamma_power, dataset, noise, shrinking, threshold, sample_size)


if __name__ == "__main__":
    app.run_server(debug=True)
