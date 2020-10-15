from base64 import b64encode

import dash
import dash_table
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output
import numpy as np
import pandas as pd

app = dash.Dash(__name__)
server = app.server

my_dfs = {
    "small": pd.DataFrame(np.random.randint(50, size=(10, 5))),
    "large": pd.DataFrame(np.random.randint(50, size=(50000, 5))),
    "do not click": pd.DataFrame(data=np.random.randint(50, size=(500000, 10))),
    "do not click": pd.DataFrame(data=np.random.randint(50, size=(5000000, 10))),
}

app.layout = html.Div(
    [
        dcc.RadioItems(
            id="choose-file",
            options=[{"label": x, "value": x} for x in my_dfs],
            value="small",
        ),
        dcc.Loading(
            html.A(html.Button("Download"), id="download", download="data.csv")
        ),
    ]
)


@app.callback(Output("download", "href"), [Input("choose-file", "value")])
def table_to_csv(val):
    df = my_dfs[val]
    return "data:text/csv," + df.to_csv(index=False)


if __name__ == "__main__":
    app.run_server(debug=True)