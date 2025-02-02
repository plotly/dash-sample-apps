from dash import Dash, html, Input, Output
import dash_pivottable

from data.data import data
from utils.components import header

app = Dash(__name__, title="Dash Pivottable")
server = app.server

app.layout = html.Div([
    header(app, "black", "Dash PivotTable"),

    html.Div([
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
    ],
    className="app-body")
])


@app.callback(
    Output("output", "children"),
    Input("table", "cols"),
    Input("table", "rows"),
    Input("table", "rowOrder"),
    Input("table", "colOrder"),
    Input("table", "aggregatorName"),
    Input("table", "rendererName"),
)
def display_props(cols, rows, row_order, col_order, aggregator, renderer):
    """
        This callback demonstrates how to access the properties of the PivotTable, which can then be used as needed.
    """
    return [
        html.H3("PivotTable parameters:"),
        html.P("cols: " + str(cols)),
        html.P("rows: " + str(rows)),
        html.P("row_order: " + str(row_order)),
        html.P("col_order: " + str(col_order)),
        html.P("aggregator: " + str(aggregator)),
        html.P("renderer: " + str(renderer)),
    ]


if __name__ == "__main__":
    app.run_server(debug=True)
