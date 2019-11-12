import lasio
import argparse
import os
import re
import pandas
import pathlib
from plotly import tools
import plotly.graph_objs as go
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt


app = dash.Dash(__name__)

app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

server = app.server

DATA_PATH = pathlib.Path(__file__).parent.resolve()


def parse_args():
    parser = argparse.ArgumentParser(description="Launch a Dash app to view a LAS log.")

    parser.add_argument(
        "lasfile",
        type=argparse.FileType(mode="r"),
        help="Log ASCII Standard (LAS) file",
    )

    parser.add_argument("--debug", "-d", action="store_true", help="enable debug mode")

    args = parser.parse_args()

    return args.lasfile, args.debug


if "DASH_APP_NAME" in os.environ:
    lasfile = open(DATA_PATH.joinpath("alcor1.las"))
    debug = True
else:
    lasfile, debug = parse_args()

lf = lasio.read(lasfile)


def generate_frontpage():
    filename = os.path.basename(lasfile.name)

    return html.Div(
        id="las-header",
        children=[
            html.Img(id="las-logo", src=app.get_asset_url("logo.png")),
            html.Div(
                id="las-header-text",
                children=[
                    html.H1("LAS Report"),
                    html.Div(
                        id="las-file-info",
                        children=[
                            html.Span(id="las-filename", children=filename),
                            html.Span(
                                " ({0})".format(
                                    lf.version["VERS"].descr
                                    if "VERS" in lf.version
                                    else "Unknown version"
                                )
                            ),
                        ],
                    ),
                ],
            ),
            html.Img(id="dash-logo", src=app.get_asset_url("dash-logo.png")),
        ],
    )


def generate_axis_title(descr, unit):
    title_words = descr.split(" ")

    current_line = ""
    lines = []
    for word in title_words:
        if len(current_line) + len(word) > 15:
            lines.append(current_line[:-1])
            current_line = ""
        current_line += "{} ".format(word)
    lines.append(current_line)

    title = "<br>".join(lines)
    title += "<br>({})".format(unit)

    return title


def generate_curves(
    height=950,
    width=800,
    bg_color="white",
    font_size=10,
    tick_font_size=8,
    line_width=0.5,
):
    # include one graph for all curves, since they have the same x axis
    yvals = "DEPT"

    cols = list(lf.curves.keys())

    plots = []

    plots.append(["BTVPVS", "DGRC"])
    plots.append(
        list(filter(lambda x: x == "EWXT" or re.search(r"R[0-9][0-9]P", x), cols))
    )
    plots.append(["ALCDLC", "ALDCLC"])
    plots.append(["TNPS"])
    plots.append(["BTCSS", "BTCS"])

    fig = tools.make_subplots(
        rows=1, cols=len(plots), shared_yaxes=True, horizontal_spacing=0
    )

    for i in range(len(plots)):
        for column in plots[i]:
            fig.append_trace(
                go.Scatter(
                    x=lf.curves[column].data,
                    y=lf.curves[yvals].data,
                    name=column,
                    line={
                        "width": line_width,
                        "dash": "dashdot" if column in plots[1] else "solid",
                    },
                ),
                row=1,
                col=i + 1,
            )
            fig["layout"]["xaxis{}".format(i + 1)].update(
                title=generate_axis_title(
                    lf.curves[plots[i][0]]["descr"], lf.curves[plots[i][0]]["unit"]
                ),
                type="log" if column in plots[1] else "linear",
            )

    fig["data"][1]["xaxis"] = "x6"
    fig["data"][6]["xaxis"] = "x7"
    fig["data"][8]["xaxis"] = "x8"
    fig["data"][11]["xaxis"] = "x9"

    # DGRC on graph 1
    fig["layout"]["xaxis6"] = dict(
        overlaying="x1",
        anchor="y",
        side="top",
        title=generate_axis_title(
            lf.curves["DGRC"]["descr"], lf.curves["DGRC"]["unit"]
        ),
    )

    # EWXT on graph 2
    fig["layout"]["xaxis7"] = dict(
        overlaying="x2",
        anchor="y",
        side="top",
        title=generate_axis_title(
            lf.curves["EWXT"]["descr"], lf.curves["EWXT"]["unit"]
        ),
    )

    # ALDCLC on graph 3
    fig["layout"]["xaxis8"] = dict(
        overlaying="x3",
        anchor="y",
        side="top",
        title=generate_axis_title(
            lf.curves["ALDCLC"]["descr"], lf.curves["ALDCLC"]["unit"]
        ),
    )

    # BTCS on graph 5
    fig["layout"]["xaxis9"] = dict(
        overlaying="x5",
        anchor="y",
        side="top",
        title=generate_axis_title(
            lf.curves["BTCS"]["descr"], lf.curves["BTCS"]["unit"]
        ),
    )

    # y axis title
    fig["layout"]["yaxis"].update(
        title=generate_axis_title(lf.curves[yvals]["descr"], lf.curves[yvals]["unit"]),
        autorange="reversed",
    )

    for axis in fig["layout"]:
        if re.search(r"[xy]axis[0-9]*", axis):
            fig["layout"][axis].update(
                mirror="all",
                automargin=True,
                showline=True,
                title=dict(font=dict(family="Arial, sans-serif", size=font_size)),
                tickfont=dict(family="Arial, sans-serif", size=tick_font_size),
            )

    fig["layout"].update(
        height=height,
        width=width,
        plot_bgcolor=bg_color,
        paper_bgcolor=bg_color,
        hovermode="y",
        legend={"font": {"size": tick_font_size}},
        margin=go.layout.Margin(r=100, t=100, b=50, l=80, autoexpand=False),
    )

    return dcc.Graph(figure=fig)


def generate_table():
    cols = ["mnemonic", "descr", "unit", "value"]
    data = {
        lf.well[i]["mnemonic"]: {col: lf.well[i][col] for col in cols}
        for i in range(len(lf.well))
    }

    df = pandas.DataFrame(data=data)

    df = df.transpose()
    df = df[cols]
    return dt.DataTable(
        id="table",
        sort_action="native",
        filter_action="native",
        row_deletable=True,
        css={
            "rule": "display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;"
        },
        style_data={"whiteSpace": "normal"},
        style_cell={
            "padding": "15px",
            "midWidth": "0px",
            "width": "25%",
            "textAlign": "center",
            "border": "white",
        },
        style_cell_conditional=[
            {"if": {"row_index": "even"}, "backgroundColor": "#f9f9f9"}
        ],
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict("rows"),
    )


app.layout = html.Div(
    [
        html.Div(id="controls", children=[html.Button("Print", id="las-print")]),
        html.Div(id="frontpage", className="page", children=generate_frontpage()),
        html.Div(
            className="section",
            children=[
                html.Div(className="section-title", children="LAS well"),
                html.Div(
                    className="page",
                    children=[
                        html.Div(id="las-table", children=generate_table()),
                        html.Div(id="las-table-print"),
                    ],
                ),
            ],
        ),
        html.Div(
            className="section",
            children=[
                html.Div(className="section-title", children="LAS curves"),
                html.Div(
                    className="page",
                    children=[html.Div(id="las-curves", children=generate_curves())],
                ),
            ],
        ),
    ]
)


@app.callback(Output("las-table-print", "children"), [Input("table", "data")])
def update_table_print(data):
    colwidths = {
        "mnemonic": "100px",
        "descr": "300px",
        "unit": "25px",
        "value": "300px",
    }
    tables_list = []
    num_tables = int(len(data) / 34) + 1  # 34 rows max per page
    for i in range(num_tables):
        table_rows = []
        for j in range(34):
            if i * 34 + j >= len(data):
                break
            table_rows.append(
                html.Tr([html.Td(data[i * 34 + j][key]) for key in data[0].keys()])
            )
        table_rows.insert(
            0,
            html.Tr(
                [
                    html.Th(key.title(), style={"width": colwidths[key]})
                    for key in data[0].keys()
                ]
            ),
        )
        tables_list.append(
            html.Div(className="tablepage", children=html.Table(table_rows))
        )
    return tables_list


if __name__ == "__main__":
    app.run_server(debug=debug)
