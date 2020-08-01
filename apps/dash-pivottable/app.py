import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_pivottable

from data import data


def Header(name, app):
    img_style = {"float": "right", "height": 40, "margin-right": 10}
    dash_logo = html.Img(src=app.get_asset_url("dash.png"), style=img_style)
    ghub_logo = html.Img(src=app.get_asset_url("github.png"), style=img_style)

    return html.Div(
        [
            html.H1(name, style={"margin": 10, "display": "inline"}),
            html.A(dash_logo, href="https://plotly.com/dash/"),
            html.A(ghub_logo, href="https://github.com/plotly/dash-pivottable"),
            html.Hr(),
        ]
    )


app = dash.Dash(__name__)
app.title = "Dash Pivottable"
server = app.server

app.layout = html.Div(
    [
        Header("Dash Pivottable", app),
        dash_pivottable.PivotTable(
            id="table",
            data=data,
            cols=["Day of Week"],
            colOrder="key_a_to_z",
            rows=["Party Size"],
            rowOrder="key_a_to_z",
            rendererName="Grouped Column Chart",
            aggregatorName="Average",
            vals=["Total Bill"],
            valueFilter={"Day of Week": {"Thursday": False}},
        ),
        html.Div(id="output"),
    ]
)


@app.callback(
    Output("output", "children"),
    [
        Input("table", "cols"),
        Input("table", "rows"),
        Input("table", "rowOrder"),
        Input("table", "colOrder"),
        Input("table", "aggregatorName"),
        Input("table", "rendererName"),
    ],
)
def display_props(cols, rows, row_order, col_order, aggregator, renderer):
    return [
        html.P(str(cols), id="columns"),
        html.P(str(rows), id="rows"),
        html.P(str(row_order), id="row_order"),
        html.P(str(col_order), id="col_order"),
        html.P(str(aggregator), id="aggregator"),
        html.P(str(renderer), id="renderer"),
    ]


if __name__ == "__main__":
    app.run_server(debug=True)
