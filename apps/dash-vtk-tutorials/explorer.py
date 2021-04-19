from importlib import import_module
from inspect import getsource
from copy import deepcopy
import json
import os

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc


def prepend_recursive(component, prefix: str) -> None:
    """in-place modifications"""
    if hasattr(component, "id"):
        if type(component.id) == str:
            component.id = prefix + component.id
        elif type(component.id) == dict:
            key = "type"
            if key in component.id:
                component.id[key] = prefix + component.id[key]

    if hasattr(component, "children") and component.children is not None:
        for child in component.children:
            prepend_recursive(child, prefix)


def prepend_list_of_dict(ls: list, prefix: str) -> list:
    new_ls = []

    for di in ls:
        di = deepcopy(di)
        try:  # is a dictionary
            di_id = json.loads(di["id"])
            key = "type"
            if key in di_id:
                di_id[key] = prefix + di_id[key]

            di["id"] = json.dumps(di_id).replace(" ", "")

        except ValueError:  # is a string
            di["id"] = prefix + di["id"]

        new_ls.append(di)
    return new_ls


def prepend_callback_map(di: dict, prefix: str) -> dict:
    new_di = {}
    for k, v in di.items():
        v = deepcopy(v)
        v["inputs"] = prepend_list_of_dict(v["inputs"], prefix)
        v["state"] = prepend_list_of_dict(v["state"], prefix)
        new_di[prefix + k] = v

    return new_di


def prepend_callback_list(ls: list, prefix: str) -> list:
    new_ls = []
    for di in ls:
        di = deepcopy(di)
        if type(di["output"]) == list:
            di["output"] = prepend_list_of_dict(di["output"], prefix)
        else:
            di["output"] = prefix + di["output"]
        di["inputs"] = prepend_list_of_dict(di["inputs"], prefix)
        di["state"] = prepend_list_of_dict(di["state"], prefix)

        new_ls.append(di)

    return new_ls


def name_to_label(x):
    return (
        x.replace("_", " ")
        .replace("t0", "Tutorial #")
        .replace("s0", "Webinar Demo #")
        .title()
        .replace("Vtk", "VTK")
    )


def Header(name, app):
    title = html.H2(name, style={"display": "inline-flex"})
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"),
        style={
            "height": 60,
            "display": "inline-flex",
            "margin-top": "-10px",
            "margin-right": "5px",
        },
    )
    link = html.A(logo, href="https://plotly.com/dash/", target="_blank")

    return html.Div([link, title])


def display_demo(name, layout, code):
    download_btn = html.A(
        html.Button(
            "Download",
            style={
                "width": "90px",
                "margin": "auto",
                "padding": "0px",
                "font-size": "10px",
                "height": "35px",
                "border-radius": "2px",
            },
        ),
        href=app.get_asset_url(name + ".py"),
        download="app.py",
        style={"position": "absolute", "top": "1.5em", "right": "1.5em"},
    )
    return html.Div(
        [
            html.Div(
                [download_btn, dcc.Markdown(f"```\n{code}\n```"),],
                style={
                    "float": "left",
                    "width": "49%",
                    "height": "85vh",
                    "overflow-y": "auto",
                    "position": "relative",
                    "background-color": "#F7FAFC",
                    "border": "1px solid #A1ACC3",
                    "border-right": "none",
                },
            ),
            html.Div(
                layout,
                style={
                    "float": "left",
                    "width": "48%",
                    "padding": "5px 1% 5px 1%",
                    "height": "calc(85vh - 10px)",
                    "overflow-y": "auto",
                    "border": "1px solid #A1ACC3",
                },
            ),
        ]
    )


prefix_ignored = []

ignored_pages = ["data", "requirements.txt"]


app = dash.Dash(
    __name__,
    suppress_callback_exceptions=True,
    external_stylesheets=[dbc.themes.COSMO],
)
server = app.server

app_subdomain = os.getenv("APP_SUBDOMAIN", "dash-vtk-tutorials")

pages = [
    p.replace(".py", "")
    for p in sorted(os.listdir("demos"))
    if p not in ignored_pages and p.endswith(".py")
]
print(pages)
modules = {p: import_module(f"demos.{p}") for p in pages}
apps = {p: m.app for p, m in modules.items()}
source_codes = {p: getsource(m) for p, m in modules.items()}
notfound_404 = html.Div(
    [
        html.H1("404"),
        "Webpage not found. Please contact us if a page is supposed to be here.",
    ]
)


app.layout = dbc.Container(
    children=[
        dbc.Row(
            style={"height": "10%", "align-items": "center"},
            children=[
                dbc.Col([Header("VTK Tutorials", app),], width=8,),
                dbc.Col(
                    dbc.Spinner(
                        dbc.Select(
                            id="app-choice",
                            placeholder="Please select an app...",
                            style={"width": "100%"},
                            options=[
                                {"label": name_to_label(x), "value": x} for x in pages
                            ],
                        ),
                    ),
                    width=4,
                ),
            ],
        ),
        html.Div(id="display", style={"height": "90%"}),
        dcc.Location(id="url", refresh=True),
    ],
    style={"height": "calc(100vh - 15px)"},
    fluid=True,
)

for k in apps:
    new_callback_map = apps[k].callback_map
    new_callback_list = apps[k]._callback_list

    # Prepend to layout IDs recursively in-place
    # if k in prefix_ignored:
    #     new_callback_map = apps[k].callback_map
    #     new_callback_list = apps[k]._callback_list
    # else:
    #     prepend_recursive(apps[k].layout, prefix=k + "-")
    #     new_callback_map = prepend_callback_map(apps[k].callback_map, prefix=k + "-")
    #     new_callback_list = prepend_callback_list(apps[k]._callback_list, prefix=k + "-")

    app.callback_map.update(new_callback_map)
    app._callback_list.extend(new_callback_list)


@app.callback(
    [Output("url", "pathname"), Output("url", "refresh")], Input("app-choice", "value")
)
def update_url(name):
    if name is None:
        return dash.no_update, dash.no_update

    return f"/{app_subdomain}/{name}", True


@app.callback(
    [Output("display", "children"), Output("app-choice", "options")],
    [Input("url", "pathname")],
)
def display_content(pathname):
    if app_subdomain in pathname:
        name = pathname.split("/")[-1]

        if name == "":
            return html.P("Please select an app from the dropdown"), dash.no_update

        elif name in pages:
            # return display_demo(
            #     name=name, layout=apps[name].layout, code=source_codes[name]
            # )
            return apps[name].layout.children, dash.no_update

        else:
            return notfound_404, dash.no_update

    return dash.no_update, dash.no_update


if __name__ == "__main__":
    app.run_server(debug=True)
