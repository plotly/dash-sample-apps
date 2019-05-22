import dash
import pathlib
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
from plotly import tools

from demo_utils import demo_components, demo_callbacks, demo_explanation

# get relative data folder
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()

LOGFILE = "examples/run_log.csv"

app = dash.Dash(__name__)
server = app.server

demo_mode = True


def div_graph(name):
    """Generates an html Div containing graph and control options for smoothing and display, given the name"""
    return html.Div(
        className="row",
        children=[
            html.Div(
                className="two columns",
                children=[
                    html.Div(
                        [
                            html.P(
                                "Smoothing:",
                                style={"font-weight": "bold", "margin-bottom": "0px"},
                            ),
                            dcc.Checklist(
                                options=[
                                    {"label": " Training", "value": "train"},
                                    {"label": " Validation", "value": "val"},
                                ],
                                values=[],
                                id=f"checklist-smoothing-options-{name}",
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
                    ),
                    html.Div(
                        [
                            html.P(
                                "Plot Display mode:",
                                style={"font-weight": "bold", "margin-bottom": "0px"},
                            ),
                            dcc.RadioItems(
                                options=[
                                    {"label": " Overlapping", "value": "overlap"},
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
                            ),
                            html.Div(id=f"div-current-{name}-value"),
                        ]
                    ),
                ],
            ),
            html.Div(id=f"div-{name}-graph", className="ten columns"),
        ],
    )


app.layout = html.Div(
    [
        # Banner display
        html.Div(
            [
                html.H2("Live Model Training Viewer", id="title"),
                html.Button("Learn More", id="learn-more-button"),
                html.Img(src=app.get_asset_url("dash-by-plotly-logo.png")),
            ],
            className="banner",
        ),
        # Body
        html.Div(
            className="container",
            children=[
                # Extract the demo components if we are in demo mode
                *demo_components(demo_mode),
                html.Div(
                    className="row",
                    id="div-interval-control",
                    children=[
                        html.Div(
                            id="div-step-display",
                            className="two columns",
                            style={"float": "right"},
                        ),
                        dcc.Dropdown(
                            id="dropdown-interval-control",
                            options=[
                                {"label": "No Updates", "value": "no"},
                                {"label": "Slow Updates", "value": "slow"},
                                {"label": "Regular Updates", "value": "regular"},
                                {"label": "Fast Updates", "value": "fast"},
                            ],
                            value="regular",
                            className="ten columns",
                            clearable=False,
                            searchable=False,
                        ),
                    ],
                ),
                dcc.Interval(id="interval-log-update", n_intervals=0),
                # Hidden Div Storing JSON-serialized dataframe of run log
                html.Div(id="run-log-storage", style={"display": "none"}),
                # The html divs storing the graphs and display parameters
                div_graph("accuracy"),
            ],
        ),
        html.Div(
            className="container",
            children=[
                div_graph("cross-entropy"),
                # Explanation for the demo version of the app
            ],
        ),
        html.Div(
            id="demo-explanation",
            children=[
                # demo_explanation(demo_mode),
            ],
        ),
    ]
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

    if run_log_json:  # exists
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

        trace_train = go.Scatter(
            x=step,
            y=y_train,
            mode="lines",
            name="Training",
            line=dict(color="rgb(54, 218, 170)"),
        )

        trace_val = go.Scatter(
            x=step,
            y=y_val,
            mode="lines",
            name="Validation",
            line=dict(color="rgb(246, 236, 145)"),
        )

        if display_mode == "separate_vertical":
            figure = tools.make_subplots(
                rows=2, cols=1, print_grid=False, shared_yaxes=True
            )

            figure.append_trace(trace_train, 1, 1)
            figure.append_trace(trace_val, 2, 1)

            figure["layout"].update(
                title=layout.title,
                margin=layout.margin,
                scene={"domain": {"x": (0.0, 0.5), "y": (0.5, 1)}},
            )

        elif display_mode == "separate_horizontal":
            figure = tools.make_subplots(
                rows=1, cols=2, shared_yaxes=True, print_grid=False
            )

            figure.append_trace(trace_train, 1, 1)
            figure.append_trace(trace_val, 1, 2)

            figure["layout"].update(title=layout.title, margin=layout.margin)

        elif display_mode == "overlap":
            figure = go.Figure(data=[trace_train, trace_val], layout=layout)

        else:
            figure = None

        return dcc.Graph(figure=figure, id=graph_id)

    return dcc.Graph(id=graph_id)


demo_callbacks(app, demo_mode)


@app.callback(
    Output("demo-explanation", "children"), [Input("learn-more-button", "n_clicks")]
)
def learn_more(n_clicks):
    if n_clicks == None:
        n_clicks = 0
        return
    else:
        if (n_clicks % 2) == 1:
            n_clicks += 1
            return (
                html.Div(className="container", children=[demo_explanation(demo_mode)]),
            )

        else:
            n_clicks += 1
            return (
                html.Div(
                    children=[
                        # demo_explanation(demo_mode),
                    ]
                ),
            )


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
        Output("run-log-storage", "children"),
        [Input("interval-log-update", "n_intervals")],
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
    Output("div-step-display", "children"), [Input("run-log-storage", "children")]
)
def update_div_step_display(run_log_json):
    if run_log_json:
        run_log_df = pd.read_json(run_log_json, orient="split")
        return html.H6(
            f"Step: {run_log_df['step'].iloc[-1]}", style={"margin-top": "3px"}
        )


@app.callback(
    Output("div-accuracy-graph", "children"),
    [
        Input("run-log-storage", "children"),
        Input("radio-display-mode-accuracy", "value"),
        Input("checklist-smoothing-options-accuracy", "values"),
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
            graph.figure.layout.yaxis["range"] = [0, 1]
        else:
            graph.figure.layout.yaxis1["range"] = [0, 1]
            graph.figure.layout.yaxis2["range"] = [0, 1]

    except AttributeError:
        pass

    return [graph]


@app.callback(
    Output("div-cross-entropy-graph", "children"),
    [
        Input("run-log-storage", "children"),
        Input("radio-display-mode-cross-entropy", "value"),
        Input("checklist-smoothing-options-cross-entropy", "values"),
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
    Output("div-current-accuracy-value", "children"),
    [Input("run-log-storage", "children")],
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
    [Input("run-log-storage", "children")],
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

