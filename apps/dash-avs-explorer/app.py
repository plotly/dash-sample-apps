import os
from textwrap import dedent

import dash_avs_ui
import dash
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html


def Header(name, app):
    title = html.H1(name, style={"margin-top": 5})
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 60}
    )
    link = html.A(logo, href="https://plotly.com/dash/")

    return dbc.Row([dbc.Col(title, md=8), dbc.Col(link, md=4)])


mapbox_token = os.getenv("MAPBOX_ACCESS_TOKEN")

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "Autonomous Vehicle Visualization"
server = app.server

STYLES = ["light", "dark", "satellite"]

DATASETS = {
    "kitti": "https://raw.githubusercontent.com/plotly/xviz-data-kitti/master/2011_09_26/",
    "nuScenes": "https://raw.githubusercontent.com/uber/xviz-data/master/nutonomy/",
}

SCENES = {
    "kitti": {
        x: f"2011_09_26_drive_{x}_sync/" for x in ["0001", "0002", "0005", "0011"]
    },
    "nuScenes": {x: f"scene-{x}/" for x in ["0006"]},
}

controls = [
    dbc.FormGroup(
        [
            dbc.Label("Map Style"),
            dbc.Select(
                id="select-style",
                options=[{"label": s.capitalize(), "value": s} for s in STYLES],
                value=STYLES[0],
            ),
        ]
    ),
    dbc.FormGroup(
        [
            dbc.Label("Dataset"),
            dbc.RadioItems(
                id="select-dataset",
                options=[{"label": x, "value": x} for x in DATASETS],
                value="kitti",
            ),
        ]
    ),
    dbc.FormGroup([dbc.Label("Scene"), dbc.RadioItems(id="select-scene", inline=True)]),
    dbc.FormGroup(
        [
            dbc.Label("Mode"),
            dbc.Checklist(
                id="select-mode",
                options=[{"label": "advanced", "value": "advanced"}],
                value=[],
                switch=True,
            ),
        ]
    ),
    dbc.FormGroup(
        [dbc.Label("About this app"), html.Br(), dbc.Button("More Info", id="btn-info")]
    ),
]

compatibility_alert = dbc.Alert(
    "This app works best in Chrome and Firefox.",
    color="info",
    is_open=True,
    dismissable=True,
)

info_collapse = dbc.Collapse(
    [
        html.Br(),
        dbc.Card(
            dcc.Markdown(
                dedent(
                    """
            This app is based on [avs.auto](http://avs.auto/). This app
            loads scenes collected from self-driving car trips, and display
            both sensor (e.g. lidar, path) and human annotated data (bounding
            boxes). You can choose between two different UIs, a basic and an advanced one.
            Both UIs were created using [streetscape.gl](https://avs.auto/#/streetscape.gl/overview/introduction)
            and are available as Dash components inside [this repo](https://github.com/plotly/dash-avs-ui).
            The source code for this demo is available [here](https://github.com/plotly/dash-avs-demo).
            """
                )
            ),
            body=True,
        ),
    ],
    id="info-collapse",
)


app.layout = dbc.Container(
    [
        Header("Dash Autonomous Driving Demo", app),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Card(controls, body=True),
                        info_collapse,
                        html.Br(),
                        compatibility_alert,
                    ],
                    md=3,
                ),
                dbc.Col(
                    dbc.Spinner(
                        html.Div(id="div-ui"), spinner_style={"margin": "auto"}
                    ),
                    md=9,
                ),
            ]
        ),
        html.Div(id="temp"),
    ],
    fluid=True,
)


@app.callback(
    Output("info-collapse", "is_open"),
    [Input("btn-info", "n_clicks")],
    [State("info-collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    [Output("select-scene", "options"), Output("select-scene", "value")],
    [Input("select-dataset", "value")],
)
def update_scene_selection(dataset):
    scenes = SCENES[dataset]
    options = [{"label": k, "value": v} for k, v in scenes.items()]
    value = list(scenes.values())[0]

    return options, value


@app.callback(
    Output("div-ui", "children"),
    [
        Input("select-scene", "value"),
        Input("select-style", "value"),
        Input("select-mode", "value"),
    ],
    [State("select-dataset", "value")],
)
def update_log(scene, style, mode, dataset):
    map_style = f"mapbox://styles/mapbox/{style}-v9"
    base_url = DATASETS[dataset]

    log = {
        "timingsFilePath": base_url + scene + "0-frame.json",
        "getFilePath": base_url + scene + "${index}-frame.glb",
        "worker": True,
        "maxConcurrency": 4,
    }

    print(log)

    if "advanced" in mode:
        UI = dash_avs_ui.AdvancedUI
    else:
        UI = dash_avs_ui.BasicUI

    return UI(
        id={"name": "avs-ui"},
        log=log,
        mapStyle=map_style,
        containerStyle={"height": "calc(100vh - 100px)", "width": "73%"},
    )


if __name__ == "__main__":
    app.run_server(debug=True)
