import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
from dash.dependencies import Input, Output, State
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


