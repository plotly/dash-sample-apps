import json

import dash
import dash_bootstrap_components as dbc
import dash_table
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go


def Header(name, app):
    title = html.H1(name, style={"margin-top": 5})
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 60}
    )
    link = html.A(logo, href="https://plotly.com/dash/")

    return dbc.Row([dbc.Col(title, md=8), dbc.Col(link, md=4)])


styles = {"pre": {"border": "thin lightgrey solid", "overflowX": "scroll"}}


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

fig = go.Figure()
text = "Click and drag here <br> to draw a rectangle <br><br> or select another shape <br>in the modebar"
fig.add_annotation(
    x=0.5, y=0.5, text=text, xref="paper", yref="paper", showarrow=False, font_size=20
)
# shape defined programatically
fig.add_shape(editable=True, x0=-1, x1=0, y0=2, y1=3, name="one", xref="x1", yref="y1")

fig.add_shape(editable=True, x0=1, x1=2, y0=3, y1=4, name="two", xref="x1", yref="y1")

# define dragmode and add modebar buttons
fig.update_layout(dragmode="drawrect")
fig_config = {"modeBarButtonsToAdd": ["drawrect", "eraseshape"]}

controls = []

app.layout = dbc.Container(
    [
        Header("App Title", app),
        html.Hr(),
        dcc.Graph(id="graph", figure=fig, config=fig_config),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H3("relayoutData's shapes"),
                        html.Pre(id="relayout-data", style=styles["pre"]),
                    ]
                ),
                dbc.Col(
                    [
                        html.H3("layoutData's shapes"),
                        html.Pre(id="layout-data", style=styles["pre"]),
                    ]
                ),
                dbc.Col(
                    [
                        dash_table.DataTable(
                            id="table",
                            data=[],
                            columns=[
                                {"id": c, "name": c}
                                for c in [
                                    "scene",
                                    "time",
                                    "order",
                                    "object",
                                    "xmin",
                                    "xmax",
                                    "ymin",
                                    "ymax",
                                ]
                            ],
                        )
                    ]
                ),
            ]
        ),
        dcc.Store(id="shapes"),
    ],
    fluid=False,
)


@app.callback(
    [Output("relayout-data", "children"), Output("layout-data", "children")],
    [Input("graph", "relayoutData")],
    [State("graph", "figure")],
)
def display_relayout_data(relayoutData, figure):
    return (
        json.dumps(relayoutData, indent=2),
        json.dumps(figure["layout"]["shapes"], indent=2),
    )


@app.callback(
    Output("table", "data"),
    [Input("graph", "relayoutData")],
    [State("graph", "figure"), State("table", "data")],
)
def update_table(relayout_data, figure, table_data):
    curr_time = 213.31231923
    scene = 1

    if relayout_data is None:
        return dash.no_update

    shapes = figure["layout"]["shapes"]

    keys = list(relayout_data.keys())
    if "shapes[" not in keys[0]:
        return dash.no_update

    i = int(keys[0].replace("shapes[", "").split("].")[0])

    print(relayout_data)
    relayout = {k.split(".")[-1]: v for k, v in relayout_data.items()}

    if i >= len(shapes):
        return dash.no_update

    filtered_table_data = [
        row
        for row in table_data
        if not (
            row["order"] == i
            and row["time"] == round(curr_time, 6)
            and row["scene"] == scene
        )
    ]

    new_shape = shapes[i]
    new = {
        "time": round(curr_time, 6),
        "scene": scene,
        "object": new_shape["name"],
        "order": i,
        "xmin": round(relayout["x0"], 1),
        "xmax": round(relayout["x1"], 1),
        "ymin": round(relayout["y0"], 1),
        "ymax": round(relayout["y1"], 1),
    }

    filtered_table_data.append(new)

    return filtered_table_data


if __name__ == "__main__":
    app.run_server(debug=True)
