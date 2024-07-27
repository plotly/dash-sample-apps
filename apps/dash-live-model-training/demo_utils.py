from dash import Dash, html, dcc, Input, Output, State, callback, callback_context
import pandas as pd
import pathlib

# get relative data folder
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()


def demo_explanation(demo_mode):
    if demo_mode:
        # Markdown files
        with open(PATH.joinpath("demo.md"), "r") as file:
            demo_md = file.read()

        return html.Div(
            html.Div([dcc.Markdown(demo_md, className="markdown")]),
            style={"margin": "10px"},
        )


def demo_callbacks(app, demo_mode):
    if demo_mode:

        @app.server.before_first_request
        def load_demo_run_logs():
            global data_dict, demo_md

            names = [
                "step",
                "train accuracy",
                "val accuracy",
                "train cross entropy",
                "val cross entropy",
            ]

            data_dict = {
                "softmax": {
                    "cifar": pd.read_csv(
                        DATA_PATH.joinpath("cifar_softmax_run_log.csv"), names=names
                    ),
                    "mnist": pd.read_csv(
                        DATA_PATH.joinpath("mnist_softmax_run_log.csv"), names=names
                    ),
                    "fashion": pd.read_csv(
                        DATA_PATH.joinpath("fashion_softmax_run_log.csv"), names=names
                    ),
                },
                "cnn": {
                    "cifar": pd.read_csv(
                        DATA_PATH.joinpath("cifar_cnn_run_log.csv"), names=names
                    ),
                    "mnist": pd.read_csv(
                        DATA_PATH.joinpath("mnist_cnn_run_log.csv"), names=names
                    ),
                    "fashion": pd.read_csv(
                        DATA_PATH.joinpath("fashion_cnn_run_log.csv"), names=names
                    ),
                },
            }

        @app.callback(
            Output("storage-simulated-run", "data"),
            [Input("interval-simulated-step", "n_intervals")],
            [
                State("dropdown-demo-dataset", "value"),
                State("dropdown-simulation-model", "value"),
            ],
        )
        def simulate_run(n_intervals, demo_dataset, simulation_model):
            if simulation_model and demo_dataset and n_intervals > 0:
                step = n_intervals * 5
                run_logs = data_dict[simulation_model][demo_dataset]

                run_below_steps = run_logs[run_logs["step"] <= step]
                json = run_below_steps.to_json(orient="split")

                return json

        @app.callback(
            Output("interval-simulated-step", "n_intervals"),
            [
                Input("dropdown-demo-dataset", "value"),
                Input("dropdown-simulation-model", "value"),
            ],
        )
        def reset_interval_simulated_step(*_):
            return 0

        @app.callback(
            Output("run-log-storage", "data"),
            [Input("interval-log-update", "n_intervals")],
            [State("storage-simulated-run", "data")],
        )
        def get_run_log(_, simulated_run):
            if simulate_run:
                return simulated_run

        @app.callback(
            Output("div-total-step-count", "children"),
            [Input("dropdown-demo-dataset", "value")],
        )
        def total_step_count(dataset_name):
            if dataset_name is not None:
                dataset = data_dict["softmax"][dataset_name]
                return html.H6(
                    f"Total Steps: {dataset['step'].iloc[-1]}",
                    style={"margin-top": "3px", "float": "right"},
                )
