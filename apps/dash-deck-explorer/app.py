from importlib import import_module
import inspect
from textwrap import dedent
import os

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from tqdm import tqdm


def Header(name, app):
    title = html.H1(name, style={"margin-top": 5})
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 60}
    )
    link = html.A(logo, href="https://plotly.com/dash/")

    return dbc.Row([dbc.Col(title, md=8), dbc.Col(link, md=4)])


def format_demo_name(demo):
    return demo.replace("usage-", "").replace("-", " ").title()


ignored_demos = ["usage-events.py", "usage-style-prop.py"]

deck_demos = [
    n.replace(".py", "").replace("usage-", "")
    for n in sorted(os.listdir("./demos"))
    if ".py" in n and n not in ignored_demos
]


deck_modules = {demo: import_module(f"demos.usage-{demo}") for demo in tqdm(deck_demos)}

print("Demos Loaded:", deck_demos)

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])
server = app.server

app_selection = dbc.FormGroup(
    [
        dbc.Label("Select Demo", width=3),
        dbc.Col(
            dbc.Select(
                id="demo-selection",
                options=[
                    {"label": demo.replace("-", " ").title(), "value": demo}
                    for demo in deck_demos
                ],
                className="form-control-plaintext",
                value=deck_demos[0],
            ),
            width=9,
        ),
    ],
    row=True,
)

tab_style = {"height": "calc(100vh - 230px)", "padding": "15px"}
# tab_style = {'max-height': 'calc(100vh - 210px)'}
tabs = dbc.Tabs(
    [
        dbc.Tab(dcc.Markdown(id="description", style=tab_style), label="Description"),
        dbc.Tab(dcc.Markdown(id="source-code", style=tab_style), label="Source Code"),
    ]
)

layout = [
    Header("Dash Deck Explorer", app),
    html.Br(),
    dcc.Location(id="url", refresh=False),
    dbc.Row(
        [
            dbc.Col(
                dbc.Card(
                    id="deck-card", style={"height": "calc(100vh - 110px)"}, body=True
                ),
                md=6,
            ),
            dbc.Col([app_selection, tabs], md=6),
        ]
    ),
]

app.layout = dbc.Container(layout, fluid=True)


@app.callback(Output("url", "pathname"), Input("demo-selection", "value"))
def update_url(name):
    return "/deck-explorer/" + name


@app.callback(
    [
        Output("deck-card", "children"),
        Output("description", "children"),
        Output("source-code", "children"),
    ],
    Input("url", "pathname"),
)
def update_demo(pathname):
    if pathname in ["/dash-deck-explorer/", None, "/"]:
        return dash.no_update

    name = pathname.split("/")[-1]

    module = deck_modules[name]
    deck_component = module.app.layout
    desc = module.__doc__
    code = f"```\n{inspect.getsource(module)}\n```"

    end = dedent(
        f"""

    -----
    * Source Code on GitHub: [Link to demo](https://github.com/plotly/dash-deck/blob/master/demos/usage-{name}.py)
    * Dash Deck for enterprises: [Contact us](https://plotly.com/contact-us)
    * Download it now: [PyPi](https://pypi.org/project/dash-deck)
    * About Dash Deck: [Readme](https://github.com/plotly/dash-deck/blob/master/README.md) | [Announcement](https://community.plotly.com/)
    """
    )

    return deck_component, desc + end, code


if __name__ == "__main__":
    app.run_server(debug=True)
