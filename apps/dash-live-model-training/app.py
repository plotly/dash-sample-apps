import dash
import pathlib
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from plotly import tools

from demo_utils import demo_callbacks, demo_explanation

# get relative data folder
DATA_PATH = pathlib.Path(__file__).parent.joinpath("data").resolve()

LOGFILE = "examples/run_log.csv"

# app = dash.Dash(__name__)
app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)

server = app.server
demo_mode = True


def div_graph(name):
    # Generates an html Div containing graph and control options for smoothing and display, given the name
    return html.Div(
        className="row",
        children=[
            html.Div(
                className="two columns",
                style={"padding-bottom": "5%"},
                children=[
                    html.Div(
                        [
                            html.Div(
                                className="graph-checkbox-smoothing",
                                children=["Smoothing:"],
                            ),
                            dcc.Checklist(
                                options=[
                                    {"label": " Training", "value": "train"},
                                    {"label": " Validation", "value": "val"},
                                ],
                                value=[],
                                id=f"checklist-smoothing-options-{name}",
                                className="checklist-smoothing",
                            ),
                        ],
                        style={"margin-top": "10px"},
                    ),
                    html.Div(
                        [
                            dcc.Slider(
                                min=0,
                                max=1,
                                step=0.05,
                                marks={i / 5: str(i / 5) for i in range(0, 6)},
                                value=0.6,
                                updatemode="drag",
                                id=f"slider-smoothing-{name}",
                            )
                        ],
                        style={"margin-bottom": "40px"},
                        className="slider-smoothing",
                    ),
                    html.Div(
                        [
                            html.P(
                                "Plot Display Mode:",
                                style={"font-weight": "bold", "margin-bottom": "0px"},
                                className="plot-display-text",
                            ),
                            html.Div(
                                [
                                    dcc.RadioItems(
                                        options=[
                                            {
                                                "label": " Overlapping",
                                                "value": "overlap",
                                            },
                                            {
                                                "label": " Separate (Vertical)",
                                                "value": "separate_vertical",
                                            },
                                            {
                                                "label": " Separate (Horizontal)",
                                                "value": "separate_horizontal",
                                            },
                                        ],
                                        value="overlap",
                                        id=f"radio-display-mode-{name}",
                                        labelStyle={"verticalAlign": "middle"},
                                        className="plot-display-radio-items",
                                    )
                                ],
                                className="radio-item-div",
                            ),
                            html.Div(id=f"div-current-{name}-value"),
                        ],
                        className="entropy-div",
                    ),
                ],
            ),
            html.Div(id=f"div-{name}-graph", className="ten columns"),
        ],
    )


app.layout = html.Div(
    style={"height": "100%"},
    children=[
        # Banner display
        html.Div(
            [
                html.H2(
                    "Live Model Training Viewer",
                    id="title",
                    className="eight columns",
                    style={"margin-left": "3%"},
                ),
                html.Button(
                    id="learn-more-button",
                    className="two columns",
                    children=["Learn More"],
                ),
                html.Img(
                    src=app.get_asset_url("dash-logo.png"),
                    className="two columns",
                    id="plotly-logo",
                ),
            ],
            className="banner row",
        ),
        html.Div(html.Div(id="demo-explanation", children=[])),
        html.Div(
            className="container",
            style={"padding": "35px 25px"},
            children=[
                dcc.Store(id="storage-simulated-run", storage_type="memory"),
                # Increment the simulation step count at a fixed time interval
                dcc.Interval(
                    id="interval-simulated-step",
                    interval=125,  # Updates every 100 milliseconds, i.e. every step takes 25 ms
                    n_intervals=0,
                ),
                html.Div(
                    className="row",
                    style={"margin": "8px 0px"},
                    children=[
                        html.Div(
                            className="twelve columns",
                            children=[
                                html.Div(
                                    className="eight columns",
                                    children=[
                                        html.Div(
                                            dcc.Dropdown(
                                                id="dropdown-demo-dataset",
                                                options=[
                                                    {
                                                        "label": "CIFAR 10",
                                                        "value": "cifar",
                                                    },
                                                    {
                                                        "label": "MNIST",
                                                        "value": "mnist",
                                                    },
                                                    {
                                                        "label": "Fashion MNIST",
                                                        "value": "fashion",
                                                    },
                                                ],
                                                value="mnist",
                                                placeholder="Select a demo dataset",
                                                searchable=False,
                                            ),
                                            className="six columns dropdown-box-first",
                                        ),
                                        html.Div(
                                            dcc.Dropdown(
                                                id="dropdown-simulation-model",
                                                options=[
                                                    {
                                                        "label": "1-Layer Neural Net",
                                                        "value": "softmax",
                                                    },
                                                    {
                                                        "label": "Simple Conv Net",
                                                        "value": "cnn",
                                                    },
                                                ],
                                                value="cnn",
                                                placeholder="Select Model to Simulate",
                                                searchable=False,
                                            ),
                                            className="six columns dropdown-box-second",
                                        ),
                                        html.Div(
                                            dcc.Dropdown(
                                                id="dropdown-interval-control",
                                                options=[
                                                    {
                                                        "label": "No Updates",
                                                        "value": "no",
                                                    },
                                                    {
                                                        "label": "Slow Updates",
                                                        "value": "slow",
                                                    },
                                                    {
                                                        "label": "Regular Updates",
                                                        "value": "regular",
                                                    },
                                                    {
                                                        "label": "Fast Updates",
                                                        "value": "fast",
                                                    },
                                                ],
                                                value="regular",
                                                className="twelve columns dropdown-box-third",
                                                clearable=False,
                                                searchable=False,
                                            )
                                        ),
                                    ],
                                ),
                                html.Div(
                                    className="four columns",
                                    id="div-interval-control",
                                    children=[
                                        html.Div(
                                            id="div-total-step-count",
                                            className="twelve columns",
                                        ),
                                        html.Div(
                                            id="div-step-display",
                                            className="twelve columns",
                                        ),
                                    ],
                                ),
                            ],
                        )
                    ],
                ),
                dcc.Interval(id="interval-log-update", n_intervals=0),
                dcc.Store(id="run-log-storage", storage_type="memory"),
            ],
        ),
        html.Div(className="container", children=[div_graph("accuracy")]),
        html.Div(
            className="container",
            style={"margin-bottom": "30px"},
            children=[div_graph("cross-entropy")],
        ),
    ],
)


def update_graph(
    graph_id,
    graph_title,
    y_train_index,
    y_val_index,
    run_log_json,
    display_mode,
    checklist_smoothing_options,
    slider_smoothing,
    yaxis_title,
):
    """
    :param graph_id: ID for Dash callbacks
    :param graph_title: Displayed on layout
    :param y_train_index: name of column index for y train we want to retrieve
    :param y_val_index: name of column index for y val we want to retrieve
    :param run_log_json: the json file containing the data
    :param display_mode: 'separate' or 'overlap'
    :param checklist_smoothing_options: 'train' or 'val'
    :param slider_smoothing: value between 0 and 1, at interval of 0.05
    :return: dcc Graph object containing the updated figures
    """

    def smooth(scalars, weight=0.6):
        last = scalars[0]
        smoothed = list()
        for point in scalars:
            smoothed_val = last * weight + (1 - weight) * point
            smoothed.append(smoothed_val)
            last = smoothed_val
        return smoothed

    if run_log_json:
        layout = go.Layout(
            title=graph_title,
            margin=go.layout.Margin(l=50, r=50, b=50, t=50),
            yaxis={"title": yaxis_title},
        )

        run_log_df = pd.read_json(run_log_json, orient="split")

        step = run_log_df["step"]
        y_train = run_log_df[y_train_index]
        y_val = run_log_df[y_val_index]

        # Apply Smoothing if needed
        if "train" in checklist_smoothing_options:
            y_train = smooth(y_train, weight=slider_smoothing)

        if "val" in checklist_smoothing_options:
            y_val = smooth(y_val, weight=slider_smoothing)

        # line charts
        trace_train = go.Scatter(
            x=step,
            y=y_train,
            mode="lines",
            name="Training",
            line=dict(color="rgb(54, 218, 170)"),
            showlegend=False,
        )

        trace_val = go.Scatter(
            x=step,
            y=y_val,
            mode="lines",
            name="Validation",
            line=dict(color="rgb(246, 236, 145)"),
            showlegend=False,
        )

        if display_mode == "separate_vertical":
            figure = tools.make_subplots(rows=2, cols=1, print_grid=False)

            figure.append_trace(trace_train, 1, 1)
            figure.append_trace(trace_val, 2, 1)

            figure["layout"].update(
                title=layout.title,
                margin=layout.margin,
                scene={"domain": {"x": (0.0, 0.5), "y": (0.5, 1)}},
            )

            figure["layout"]["yaxis1"].update(title=yaxis_title)
            figure["layout"]["yaxis2"].update(title=yaxis_title)

        elif display_mode == "separate_horizontal":
            figure = tools.make_subplots(
                rows=1, cols=2, print_grid=False, shared_xaxes=True
            )

            figure.append_trace(trace_train, 1, 1)
            figure.append_trace(trace_val, 1, 2)

            figure["layout"].update(title=layout.title, margin=layout.margin)
            figure["layout"]["yaxis1"].update(title=yaxis_title)
            figure["layout"]["yaxis2"].update(title=yaxis_title)

        elif display_mode == "overlap":
            figure = go.Figure(data=[trace_train, trace_val], layout=layout)

        else:
            figure = None

        return dcc.Graph(figure=figure, id=graph_id)

    return dcc.Graph(id=graph_id)


demo_callbacks(app, demo_mode)


@app.callback(
    [Output("demo-explanation", "children"), Output("learn-more-button", "children")],
    [Input("learn-more-button", "n_clicks")],
)
def learn_more(n_clicks):
    if n_clicks is None:
        n_clicks = 0
    if (n_clicks % 2) == 1:
        n_clicks += 1
        return (
            html.Div(
                className="container",
                style={"margin-bottom": "30px"},
                children=[demo_explanation(demo_mode)],
            ),
            "Close",
        )

    n_clicks += 1
    return (html.Div(), "Learn More")


@app.callback(
    Output("interval-log-update", "interval"),
    [Input("dropdown-interval-control", "value")],
)
def update_interval_log_update(interval_rate):
    if interval_rate == "fast":
        return 500

    elif interval_rate == "regular":
        return 1000

    elif interval_rate == "slow":
        return 5 * 1000

    # Refreshes every 24 hours
    elif interval_rate == "no":
        return 24 * 60 * 60 * 1000


if not demo_mode:

    @app.callback(
        Output("run-log-storage", "data"), [Input("interval-log-update", "n_intervals")]
    )
    def get_run_log(_):
        names = [
            "step",
            "train accuracy",
            "val accuracy",
            "train cross entropy",
            "val cross entropy",
        ]

        try:
            run_log_df = pd.read_csv(DATA_PATH.joinpath(LOGFILE), names=names)
            json = run_log_df.to_json(orient="split")
        except FileNotFoundError as error:
            print(error)
            print(
                "Please verify if the csv file generated by your model is placed in the correct directory."
            )
            return None

        return json


@app.callback(
    Output("div-step-display", "children"), [Input("run-log-storage", "data")]
)
def update_div_step_display(run_log_json):
    if run_log_json:
        run_log_df = pd.read_json(run_log_json, orient="split")
        return html.H6(
            f"Step: {run_log_df['step'].iloc[-1]}",
            style={"margin-top": "3px", "float": "right"},
        )


@app.callback(
    Output("div-accuracy-graph", "children"),
    [
        Input("run-log-storage", "data"),
        Input("radio-display-mode-accuracy", "value"),
        Input("checklist-smoothing-options-accuracy", "value"),
        Input("slider-smoothing-accuracy", "value"),
    ],
)
def update_accuracy_graph(
    run_log_json, display_mode, checklist_smoothing_options, slider_smoothing
):
    graph = update_graph(
        "accuracy-graph",
        "Prediction Accuracy",
        "train accuracy",
        "val accuracy",
        run_log_json,
        display_mode,
        checklist_smoothing_options,
        slider_smoothing,
        "Accuracy",
    )

    try:
        if display_mode in ["separate_horizontal", "overlap"]:
            graph.figure.layout.yaxis["range"] = [0, 1.3]
            graph.figure.layout.yaxis2["range"] = [0, 1.3]
        else:
            graph.figure.layout.yaxis1["range"] = [0, 1.3]
            graph.figure.layout.yaxis2["range"] = [0, 1.3]

    except AttributeError:
        pass

    return [graph]


@app.callback(
    Output("div-cross-entropy-graph", "children"),
    [
        Input("run-log-storage", "data"),
        Input("radio-display-mode-cross-entropy", "value"),
        Input("checklist-smoothing-options-cross-entropy", "value"),
        Input("slider-smoothing-cross-entropy", "value"),
    ],
)
def update_cross_entropy_graph(
    run_log_json, display_mode, checklist_smoothing_options, slider_smoothing
):
    graph = update_graph(
        "cross-entropy-graph",
        "Cross Entropy Loss",
        "train cross entropy",
        "val cross entropy",
        run_log_json,
        display_mode,
        checklist_smoothing_options,
        slider_smoothing,
        "Loss",
    )
    return [graph]


@app.callback(
    Output("div-current-accuracy-value", "children"), [Input("run-log-storage", "data")]
)
def update_div_current_accuracy_value(run_log_json):
    if run_log_json:
        run_log_df = pd.read_json(run_log_json, orient="split")
        return [
            html.P(
                "Current Accuracy:",
                style={
                    "font-weight": "bold",
                    "margin-top": "15px",
                    "margin-bottom": "0px",
                },
            ),
            html.Div(f"Training: {run_log_df['train accuracy'].iloc[-1]:.4f}"),
            html.Div(f"Validation: {run_log_df['val accuracy'].iloc[-1]:.4f}"),
        ]


@app.callback(
    Output("div-current-cross-entropy-value", "children"),
    [Input("run-log-storage", "data")],
)
def update_div_current_cross_entropy_value(run_log_json):
    if run_log_json:
        run_log_df = pd.read_json(run_log_json, orient="split")
        return [
            html.P(
                "Current Loss:",
                style={
                    "font-weight": "bold",
                    "margin-top": "15px",
                    "margin-bottom": "0px",
                },
            ),
            html.Div(f"Training: {run_log_df['train cross entropy'].iloc[-1]:.4f}"),
            html.Div(f"Validation: {run_log_df['val cross entropy'].iloc[-1]:.4f}"),
        ]


# Running the server
if __name__ == "__main__":
    app.run_server(debug=True)
