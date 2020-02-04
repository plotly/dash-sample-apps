import datetime as dt
import re
from textwrap import dedent

import dash
import dash_core_components as dcc
import dash_html_components as html

import dash_table
import pandas as pd

import sqlalchemy as db
from dash.dependencies import Input, Output, State

# SQL Engine
disk_engine = db.create_engine(
    "sqlite:///data_entry.db", connect_args={"check_same_thread": False}
)
connection = disk_engine.connect()
metadata = db.MetaData()
SQL_table = db.Table(
    "data_entries",
    metadata,
    db.Column("operator_id", db.String(255)),
    db.Column("reagent", db.String(255)),
    db.Column("time_elapsed", db.String(255)),
    db.Column("amount_pipetted", db.Float()),
    db.Column("time_stamp", db.DATETIME, primary_key=True),
)

app = dash.Dash(__name__)
server = app.server

app.config.suppress_callback_exceptions = False

app.layout = html.Div(
    [
        html.Div(
            [
                html.Img(src="assets/dash-logo.png", className="app__logo"),
                html.H4("MES FORM SUBMISSION", className="header__text"),
            ],
            className="app__header",
        ),
        html.Div(
            [
                dcc.Tabs(
                    id="tabs",
                    value="data-entry",
                    children=[
                        dcc.Tab(
                            label="DATA ENTRY",
                            value="data-entry",
                            children=[
                                html.Div(
                                    [
                                        html.Div(
                                            [
                                                html.P(
                                                    "Operator ID:",
                                                    className="input__heading",
                                                ),
                                                dcc.Input(
                                                    id="enter-operator-id",
                                                    value="User",
                                                    placeholder="Enter operator ID.",
                                                    className="oper__input",
                                                ),
                                            ],
                                            className="input__container",
                                        ),
                                        html.Div(
                                            [
                                                html.P(
                                                    "Reagent:",
                                                    className="input__heading",
                                                ),
                                                dcc.Dropdown(
                                                    id="select-reagent",
                                                    options=[
                                                        {"label": i, "value": i}
                                                        for i in [
                                                            "Acetone",
                                                            "Isopropanol",
                                                            "Toluene",
                                                        ]
                                                    ],
                                                    value="Acetone",
                                                    placeholder="Select reagent.",
                                                    className="reag__select",
                                                ),
                                            ],
                                            className="dropdown__container",
                                        ),
                                        html.Div(
                                            [
                                                html.P(
                                                    "Time (HH:MM:SS)",
                                                    className="input__heading",
                                                ),
                                                dcc.Input(
                                                    id="enter-time",
                                                    value="00:00:00",
                                                    placeholder="Enter Time (HH:MM:SS)",
                                                    className="time__input",
                                                ),
                                            ],
                                            className="input__container",
                                        ),
                                        html.Div(
                                            [
                                                html.P(
                                                    "Amount Pipetted (mL):",
                                                    className="input__heading",
                                                ),
                                                dcc.Input(
                                                    id="enter-pipetted",
                                                    type="number",
                                                    value="0",
                                                    placeholder="Enter amount in (mL).",
                                                    className="vol__input",
                                                ),
                                            ],
                                            className="input__container",
                                        ),
                                        html.Div(
                                            [
                                                html.Button(
                                                    "SUBMIT ENTRY",
                                                    id="submit-entry",
                                                    className="submit__button",
                                                )
                                            ]
                                        ),
                                    ],
                                    className="container__1",
                                )
                            ],
                        ),
                        dcc.Tab(
                            label="VIEW DATA ENTRY",
                            value="view-entry",
                            children=[
                                html.Div(
                                    [
                                        dcc.Graph(
                                            id="entry-graph",
                                            figure={
                                                "data": [],
                                                "layout": {
                                                    "title": "Timestamp vs. Amount Pipetted (mL)",
                                                    "xaxis": {
                                                        "title": "Timestamp (YYYY-MM-DD HH-MM-SS)"
                                                    },
                                                    "yaxis": {
                                                        "title": "Amount Pipetted (mL)"
                                                    },
                                                    "annotations": [],
                                                    "margin": {"l": 50},
                                                },
                                            },
                                            config={"editable": True},
                                            className="graph__1",
                                        ),
                                        html.Br(),
                                        html.Div(
                                            children=[
                                                html.Button(
                                                    "ADD ANNOTATION",
                                                    id="add-button",
                                                    style={"margin-right": "2.5%"},
                                                    className="add__button",
                                                ),
                                                html.Button(
                                                    "DELETE ANNOTATION",
                                                    id="delete-button",
                                                    style={"margin-right": "2.5%"},
                                                    className="del__button",
                                                ),
                                                html.Button(
                                                    "CLEAR DATA",
                                                    id="clear-button",
                                                    style={"margin-right": "2.5%"},
                                                    className="clear__button",
                                                ),
                                            ]
                                        ),
                                    ],
                                    className="graph_container",
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            [
                                                dash_table.DataTable(
                                                    id="entry-table",
                                                    style_cell={
                                                        "minWidth": "0px",
                                                        "maxWidth": "180px",
                                                        "whiteSpace": "normal",
                                                    },
                                                )
                                            ],
                                            className="table__1",
                                        )
                                    ],
                                    className="table__container",
                                ),
                            ],
                        ),
                    ],
                )
            ],
            className="tabs__container",
        ),
    ],
    className="app__container",
)


@app.callback(
    Output("tabs", "value"),
    [Input("submit-entry", "n_clicks")],
    [
        State("enter-operator-id", "value"),
        State("select-reagent", "value"),
        State("enter-time", "value"),
        State("enter-pipetted", "value"),
    ],
)
def entry_to_db(submit_entry, operator_id, reagent, time_elapsed, amount_pipetd):
    if submit_entry:
        sample_entry = [
            {
                "operator_id": operator_id,
                "reagent": reagent,
                "time_elapsed": time_elapsed,
                "amount_pipetted": amount_pipetd,
                "time_stamp": dt.datetime.now(),
            }
        ]
        insert_entry = connection.execute(db.insert(SQL_table), sample_entry)
        return "view-entry"
    raise dash.exceptions.PreventUpdate


@app.callback(
    [
        Output("entry-graph", "figure"),
        Output("entry-table", "columns"),
        Output("entry-table", "data"),
    ],
    [
        Input("tabs", "value"),
        Input("add-button", "n_clicks"),
        Input("delete-button", "n_clicks"),
        Input("clear-button", "n_clicks"),
        Input("entry-graph", "relayoutData"),
    ],
    [
        State("entry-graph", "figure"),
        State("entry-table", "columns"),
        State("entry-table", "data"),
    ],
)
def entry_table(
    tab,
    add_button,
    delete_button,
    clear_button,
    relayData,
    figure,
    entry_columns,
    entry_data,
):
    if tab == "view-entry":
        try:
            callback = dash.callback_context.triggered[0]["prop_id"]
        except:
            pass

        df = pd.read_sql_query(
            dedent(
                """
        SELECT * from data_entries
        """
            ),
            disk_engine,
        )

        if "add-button" in callback:
            if len(df) == 0:
                return figure, entry_columns, entry_data
            else:
                annotations = figure["layout"]["annotations"]
                add_annote = dict(
                    x=df.iloc[0]["time_stamp"],
                    y=df.iloc[0]["amount_pipetted"],
                    xref="x",
                    yref="y",
                    text="Add Text",
                    showarrow=True,
                    arrowhead=7,
                    ax=0,
                    ay=-40,
                )
                annotations.append(add_annote)
                figure["layout"]["annotations"] = annotations
                return figure, entry_columns, entry_data

        elif "delete-button" in callback:
            annotations = figure["layout"]["annotations"]
            try:
                del annotations[-1]
            except:
                raise dash.exceptions.PreventUpdate
            figure["layout"]["annotations"] = annotations
            return figure, entry_columns, entry_data
        elif "clear-button" in callback:
            figure["layout"] = {
                "title": "Timestamp vs. Amount Pipetted (mL)",
                "xaxis": {"title": "Timestamp (YYYY-MM-DD HH-MM-SS)"},
                "yaxis": {"title": "Amount Pipetted (mL)"},
                "annotations": [],
                "margin": {"l": 50},
            }
            figure["data"] = []
            query = db.delete(SQL_table)
            results = connection.execute(query)
            return figure, entry_columns, []
        elif "entry-graph" in callback:
            if any("annotations" in key for key in relayData):
                relay_dict = {
                    key.split("annotations")[1]: value
                    for key, value in relayData.items()
                    if "annotations" in key.lower()
                }
                for k, v in relay_dict.items():
                    annot_idx = int(re.findall("\[(.*?)\]", k)[0])
                    annot_param = k.split(".")[1]
                    figure["layout"]["annotations"][annot_idx][annot_param] = v
                return figure, entry_columns, entry_data
            raise dash.exceptions.PreventUpdate
        figure["data"] = [
            {
                "x": df.time_stamp,
                "y": df.amount_pipetted,
                "mode": "lines+markers",
                "fill": "tozeroy",
            }
        ]
        columns = [{"name": i, "id": i} for i in df.columns]
        data = df.to_dict("records")
        return figure, columns, data
    raise dash.exceptions.PreventUpdate


if __name__ == "__main__":
    app.run_server(debug=True)
