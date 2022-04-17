# -*- utf-8 -*-

import dash_bootstrap_components as dbc
from dash import (
    Input,
    Output,
    State,
    html,
    Dash,
    dcc,
    dash_table,
    get_asset_url,
    ALL,
    MATCH,
    ClientsideFunction,
    callback_context,
    no_update,
)
from dash.exceptions import PreventUpdate

# ==========================================

# ==========================================
import time
import numpy as np
import pandas as pd

# ==========================================
from utils.sampling import sampling, df_split, data_split
from utils.modeling import modeling
from utils.charting import prediction_plot, confusion_matrix_plot, roc_curve_plot

# ==========================================
from utils.handle_func import *
from utils.split_components import *

# ==========================================

# ==========================================

# ==========================================

# ==========================================

# ==========================================

# ==========================================

meta_tags = [  # Please modify before deploying.
    dict(
        name="description",
        content="This is a 2022 replica of dash-svm. I used dash's recent new components and new features, and I've added new page elements.",
    ),
    dict(property="og:url", content="https://plotly.com/"),
    dict(property="og:type", content="website"),
    dict(property="og:title", content="Support Vector Machine Explorer"),
    dict(
        property="og:description",
        content="This is a 2022 replica of dash-svm. I used dash's recent new components and new features, and I've added new page elements.",
    ),
    dict(
        property="og:image",
        content="https://plotly-marketing-website.cdn.prismic.io/plotly-marketing-website/948b6663-9429-4bd6-a4cc-cb33231d4532_logo-plotly.svg",
    ),
    dict(name="twitter:card", content="summary_large_image"),
    dict(property="twitter:domain", content="www.plotly.com/"),
    dict(property="twitter:url", content="https://plotly.com/"),
    dict(name="twitter:title", content="Support Vector Machines"),
    dict(
        name="twitter:description",
        content="This is a 2022 replica of dash-svm. I used dash's recent new components and new features, and I've added new page elements.",
    ),
    dict(
        name="twitter:image",
        content="https://plotly-marketing-website.cdn.prismic.io/plotly-marketing-website/948b6663-9429-4bd6-a4cc-cb33231d4532_logo-plotly.svg",
    ),
]
# ==========================================

app = Dash(
    __name__,
    title="Support Vector Machine Explorer",
    update_title="Eating...",
    # external_scripts=external_scripts,
    meta_tags=meta_tags,
    external_stylesheets=[dbc.themes.FLATLY],
)

server = app.server  # for deployment

app.config.suppress_callback_exceptions = True
# ==========================================
# ==========================================

# =============layout=======================

app.layout = html.Div(
    [
        dbc.Navbar(
            dbc.Container(
                [
                    html.A(
                        dbc.Row(
                            [
                                dbc.Col(
                                    html.Img(
                                        src="https://plotly-marketing-website.cdn.prismic.io/plotly-marketing-website/948b6663-9429-4bd6-a4cc-cb33231d4532_logo-plotly.svg",
                                        height="30px",
                                    )
                                ),
                                dbc.Col(
                                    dbc.NavbarBrand(
                                        "Support Vector Machine Explorer",
                                        className="ms-2",
                                    )
                                ),
                            ],
                            align="center",
                            className="g-0",
                        ),
                        href="/",
                        style={"textDecoration": "none"},
                    ),
                ]
            )
        ),
        dbc.Row(
            [
                dbc.Col(
                    dbc.Container(
                        [
                            html.Br(),
                            dbc.Row(
                                dbc.Col(
                                    fig_0 := dcc.Graph(
                                        config=dict(displayModeBar="hover")
                                    )
                                )
                            ),
                            html.Br(),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        fig_1 := dcc.Graph(
                                            config=dict(displayModeBar=False)
                                        ),
                                        width=6,
                                        align="center",
                                    ),
                                    dbc.Col(
                                        fig_2 := dcc.Graph(
                                            config=dict(displayModeBar=False)
                                        ),
                                        width=6,
                                        align="center",
                                    ),
                                ]
                            ),
                            html.Br(),
                            dbc.Row(
                                alert := dbc.Alert(
                                    is_open=False,
                                    dismissable=True,
                                    duration=2000,
                                    style={
                                        "padding-left": "45px",
                                        "margin-top": "55px",
                                    },
                                )
                            ),
                        ],
                    ),
                    width=8,
                ),
                dbc.Col(
                    dbc.Container(
                        [
                            html.Br(),
                            uploader_btn,
                            html.Br(),
                            threshold,
                            html.Br(),
                            kernel,
                            html.Br(),
                            formula,
                            html.Br(),
                            cost,
                            html.Br(),
                            degree,
                            html.Br(),
                            gamma,
                            html.Br(),
                            html.Br(),
                            shrinking,
                            html.Br(),
                        ]
                    ),
                    width=4,
                ),
            ]
        ),
        dbc.Row(
            [
                svm_params := html.Div(style={"display": "none"}),
                # It seems like the current version of dcc.Store doesn't support assignment expressions.
                dcc.Store(
                    id="datasets_params",
                    storage_type="memory",
                    data=[
                        [
                            {"type": "dataset_parameter", "index": "dataset"},
                            {"type": "dataset_parameter", "index": "sample_size"},
                            {"type": "dataset_parameter", "index": "noise"},
                            {"type": "dataset_parameter", "index": "test_size"},
                        ],
                        ["LS", 100, 0.4, 0.25],
                        [],
                        [],
                        [],
                        [],
                        [],
                        [],
                        "tab-0",
                    ],
                ),
            ]
        ),
    ]
)

# ==========================================

# ================callbacks=================
app.callback(
    [
        Output({"type": "uploader_parameter", "index": "x"}, "options"),
        Output({"type": "uploader_parameter", "index": "y"}, "options"),
        Output({"type": "uploader_parameter", "index": "c"}, "options"),
    ],
    [Input({"type": "uploader_parameter", "index": "uploader"}, "contents")],
    [State({"type": "uploader_parameter", "index": "uploader"}, "filename")],
    prevent_initial_call=True,
)(
    lambda data, filename: 3 * [parse_contents(data, filename, header=True)]
    if data
    else 3 * [no_update]
)


@app.callback(
    Output({"type": "tabs-table", "id": "t2"}, "data"),
    [Input({"type": "uploader_parameter", "index": ALL}, "value")],
    [
        State({"type": "uploader_parameter", "index": ALL}, "id"),
        State({"type": "uploader_parameter", "index": "uploader"}, "contents"),
        State({"type": "uploader_parameter", "index": "uploader"}, "filename"),
    ],
)
def update_output(xyc, idx, uploaded_data, filename):
    data_params = {j["index"]: xyc[i] for i, j in enumerate(idx)}
    if (
        not all([v for k, v in data_params.items() if k not in ["uploader"]])
        or not uploaded_data
    ):
        raise PreventUpdate
    else:

        df0 = parse_contents(
            uploaded_data,
            filename,
            header=False,
            usecols=[
                v for k, v in data_params.items() if k not in ["uploader", "test_size"]
            ],
        )

        df0.rename(
            columns={v: k for k, v in data_params.items() if k != "test_size"},
            inplace=1,
        )

        data_params |= {"df": df0}

        data = df_split(**data_params)

        split_data = data[0]

        df1 = pd.DataFrame(split_data[0], columns=["x", "y"])
        df1["c"] = split_data[2]
        df1["s"] = "TRAIN"
        df2 = pd.DataFrame(split_data[1], columns=["x", "y"])
        df2["c"] = split_data[3]
        df2["s"] = "TEST"
        df = pd.concat([df1, df2])

        return df.to_dict("records")


@app.callback(
    Output({"type": "tabs-table", "id": "t3"}, "data"),
    [
        Input({"type": "canvas_parameter", "index": "test_size"}, "value"),
        Input({"type": "canvas_parameter", "index": "canvas"}, "json_data"),
    ],
)
def canvas_output(test_size, canvas_data):
    if not canvas_data:
        raise PreventUpdate

    X, y = handle_json(canvas_data)

    params = {"test_size": test_size, "X": X, "y": y}

    data = data_split(**params)

    split_data = data[0]

    df1 = pd.DataFrame(split_data[0], columns=["x", "y"])
    df1["c"] = split_data[2]
    df1["s"] = "TRAIN"
    df2 = pd.DataFrame(split_data[1], columns=["x", "y"])
    df2["c"] = split_data[3]
    df2["s"] = "TEST"
    df = pd.concat([df1, df2])

    return df.to_dict("records")


@app.callback(
    Output({"type": "tabs-table", "id": "t1"}, "data"),
    [Input({"type": "dataset_parameter", "index": ALL}, "value")],
    [State({"type": "dataset_parameter", "index": ALL}, "id")],
)
def generate_data(value, idx):

    data_params = {j["index"]: value[i] for i, j in enumerate(idx)}

    data = sampling(**data_params)

    split_data = data[0]

    df1 = pd.DataFrame(split_data[0], columns=["x", "y"])
    df1["c"] = split_data[2]
    df1["s"] = "TRAIN"
    df2 = pd.DataFrame(split_data[1], columns=["x", "y"])
    df2["c"] = split_data[3]
    df2["s"] = "TEST"
    df = pd.concat([df1, df2])
    return df.to_dict("records")


@app.callback(
    Output(svm_params, "children"),
    [Input({"type": "svm_parameter", "index": ALL}, "value")],
    [State({"type": "svm_parameter", "index": ALL}, "id")],
)
def params_update(value, idx):
    df = pd.DataFrame({"index": [i["index"] for i in idx], "value": value})
    return dash_table.DataTable(
        df.to_dict("records"), [{"name": i, "id": i} for i in df.columns]
    )


@app.callback(
    [
        Output(fig_0, "figure"),
        Output(fig_1, "figure"),
        Output(fig_2, "figure"),
        Output(alert, "children"),
        Output(alert, "is_open"),
    ],
    [
        Input(save_btn, "n_clicks"),
        Input({"type": "svm_parameter", "index": ALL}, "value"),
    ],
    [
        State({"type": "svm_parameter", "index": ALL}, "id"),
        State({"type": "dataset_parameter", "index": ALL}, "id"),
        State({"type": "dataset_parameter", "index": ALL}, "value"),
        State({"type": "uploader_parameter", "index": ALL}, "id"),
        State({"type": "uploader_parameter", "index": ALL}, "value"),
        State({"type": "uploader_parameter", "index": ALL}, "contents"),
        State({"type": "uploader_parameter", "index": ALL}, "filename"),
        State({"type": "canvas_parameter", "index": ALL}, "json_data"),
        State({"type": "canvas_parameter", "index": ALL}, "value"),
        State(tabs, "active_tab"),
        State("datasets_params", "data"),
    ],
)
def params_update(
    n_clicks,
    value,
    idx,
    tab_1_idx,
    tab_1_values,
    tab_2_idx,
    tab_2_values,
    uploaded_data,
    filename,
    canvas_data,
    tab_3_params,
    at,
    tabs_cache,
):
    t1 = time.perf_counter()

    if callback_context.triggered[0]["prop_id"].split(".")[0] != save_btn.id:
        [
            tab_1_idx,
            tab_1_values,
            tab_2_idx,
            tab_2_values,
            uploaded_data,
            filename,
            canvas_data,
            tab_3_params,
            at,
        ] = tabs_cache

    if at == "tab-0":
        data_1_params = {j["index"]: tab_1_values[i] for i, j in enumerate(tab_1_idx)}

        data = sampling(**data_1_params)

    elif at == "tab-1":
        data_2_params = {j["index"]: tab_2_values[i] for i, j in enumerate(tab_2_idx)}

        df0 = parse_contents(
            list(filter(None, uploaded_data))[0],
            #  This is not recommended. I did this to bypass a bug in the current version where when I wanted to call back a component placed in dbc.Tab, the process would report an error saying that the component could not be found.
            list(filter(None, filename))[0],
            header=False,
            usecols=[
                v
                for k, v in data_2_params.items()
                if k not in ["uploader", "test_size"]
            ],
        )

        df0.rename(
            columns={v: k for k, v in data_2_params.items() if k != "test_size"},
            inplace=1,
        )

        data_2_params |= {"df": df0}

        data = df_split(**data_2_params)

    elif at == "tab-2":
        X, y = handle_json(
            list(filter(None, canvas_data))[0],
        )

        split_params = {
            "test_size": list(filter(None, tab_3_params))[0],
            "X": X,
            "y": y,
        }

        data = data_split(**split_params)

    if not data:
        raise PreventUpdate

    params = {j["index"]: value[i] for i, j in enumerate(idx)}

    params |= {"data": data}
    params |= {"cost": 10 ** params["cost_power"] * params["cost_coef"]}
    params |= {"gamma": 10 ** params["gamma_power"] * params["gamma_coef"]}

    model = modeling(**params)
    params |= {"model": model}

    fig_0 = prediction_plot(**params)
    fig_1 = roc_curve_plot(**params)
    fig_2 = confusion_matrix_plot(**params)

    t2 = time.perf_counter()

    alert_info = "Takes {:.3} seconds".format(t2 - t1)

    return fig_0, fig_1, fig_2, alert_info, True


@app.callback(Output(offcanvas_content, "children"), [Input(tabs, "active_tab")])
def switch_tab(at):
    if at == "tab-0":
        return tab_1_content
    elif at == "tab-1":
        return tab_2_content
    elif at == "tab-2":
        return tab_3_content
    return html.P("Something is wrong...")


# ==========================================

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="reset_threshold"),
    Output({"type": "svm_parameter", "index": "threshold"}, "value"),
    [Input(threshold_btn, "n_clicks")],
    [State(fig_0, "figure")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="canvas_toggle"),
    Output({"type": "canvas_parameter", "index": "canvas"}, "lineColor"),
    [Input({"type": "canvas_parameter", "index": "toggle"}, "n_clicks")],
    [State({"type": "canvas_parameter", "index": "canvas"}, "lineColor")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="datasets_params_store"),
    Output("datasets_params", "data"),
    [Input(save_btn, "n_clicks")],
    [
        # tab_1_idx, tab_1_values, tab_2_idx,tab_2_values, uploaded_data, filename, canvas_data, tab_3_params
        State({"type": "dataset_parameter", "index": ALL}, "id"),
        State({"type": "dataset_parameter", "index": ALL}, "value"),
        State({"type": "uploader_parameter", "index": ALL}, "id"),
        State({"type": "uploader_parameter", "index": ALL}, "value"),
        State({"type": "uploader_parameter", "index": ALL}, "contents"),
        State({"type": "uploader_parameter", "index": ALL}, "filename"),
        State({"type": "canvas_parameter", "index": ALL}, "json_data"),
        State({"type": "canvas_parameter", "index": ALL}, "value"),
        State(tabs, "active_tab"),
    ],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="open_offcanvas"),
    Output(offcanvas, "is_open"),
    [Input(offcanvas_btn, "n_clicks"), Input(save_btn, "n_clicks")],
    [State(offcanvas, "is_open")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="btn_disabled"),
    Output(save_btn, "disabled"),
    [Input({"type": "tabs-table", "id": ALL}, "is_loading")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="disable_param_degree"),
    Output({"type": "svm_parameter", "index": "degree"}, "disabled"),
    [Input({"type": "svm_parameter", "index": "kernel"}, "value")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="kernel_formula"),
    Output(latex_formula, "children"),
    [Input({"type": "svm_parameter", "index": "kernel"}, "value")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="disable_param_gamma"),
    [
        Output({"type": "svm_parameter", "index": "gamma_power"}, "disabled"),
        Output({"type": "svm_parameter", "index": "gamma_coef"}, "disabled"),
    ],
    [Input({"type": "svm_parameter", "index": "kernel"}, "value")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="scale_param"),
    Output({"type": "svm_parameter", "index": "cost_coef"}, "marks"),
    [Input({"type": "svm_parameter", "index": "cost_power"}, "value")],
)

app.clientside_callback(
    ClientsideFunction(namespace="clientside", function_name="scale_param"),
    Output({"type": "svm_parameter", "index": "gamma_coef"}, "marks"),
    Input({"type": "svm_parameter", "index": "gamma_power"}, "value"),
)

# ==========================================
# ==========================================
# ==========================================
# ==========================================

if __name__ == "__main__":
    app.run_server(debug=True)
