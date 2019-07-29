import os
import dash
import dash_core_components as dcc
import dash_html_components as html
from datetime import datetime
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from sqlalchemy import create_engine
from components import Header, Row, Column


# connect to database
con_string = f"postgresql+pg8000{os.getenv('DATABASE_URL').lstrip('postgres')}"
engine = create_engine(con_string)
connection = engine.connect()

app = dash.Dash(__name__)

# expose the server variable for deployments
server = app.server

# standard Dash app code below
app.layout = html.Div(
    [
        Header("PostgreSQL DEMO APP", app),
        Row(
            [
                Column(
                    [
                        html.Div(
                            [
                                html.H5("SQL Editor"),
                                dcc.Textarea(
                                    id="sql__textarea__input",
                                    placeholder="> Enter a SQL command...",
                                    value="",
                                    className="sql__textarea",
                                ),
                                html.Div(
                                    [
                                        html.Button(
                                            "Execute query",
                                            id="sql__execute__handler",
                                            className="sql__button",
                                            n_clicks_timestamp=0,
                                            title="Attempt to execute SQL query in the editor.",
                                        ),
                                        html.Button(
                                            "Reset DB",
                                            id="sql__reset__handler",
                                            className="sql__button--red",
                                            n_clicks_timestamp=0,
                                            title="Reset the database schema and drop all tables.",
                                        ),
                                    ],
                                    className="btn--group",
                                ),
                            ]
                        )
                    ],
                    width=7,
                ),
                Column(
                    [
                        html.Div(
                            [
                                dcc.Markdown(
                                    """
                    ##### Sample Queries

                    ```
                    CREATE TABLE Account (
                        Email varchar(255) PRIMARY KEY,
                        Username varchar(255)
                    );
                    ```

                    ```
                    INSERT INTO Account (Email, Username)
                    VALUES ('chriswoo@gmail.com', 'Chris');
                    ```

                    ```
                    SELECT * FROM Account;
                    ```

                    ```
                    DROP TABLE Account;
                    ```

                """
                                )
                            ]
                        ),
                        html.Div(
                            [
                                html.H5("Log History"),
                                html.Pre(
                                    ["No history."],
                                    id="sql__history__output",
                                    className="sql__history__pre",
                                ),
                            ],
                            className="mt-2 log__history",
                        ),
                    ],
                    width=5,
                ),
            ]
        ),
        Row(
            [
                html.H5("Query Output"),
                html.Table(id="table__output"),
                html.Div(
                    [html.Pre(["No output."], id="sql__output", className="sql__pre")]
                ),
            ]
        ),
        dcc.Store(id="log_storage"),
    ],
    className="app__container",
)


@app.callback(
    [
        Output("sql__output", "children"),
        Output("table__output", "children"),
        Output("sql__history__output", "children"),
        Output("log_storage", "data"),
    ],
    [
        Input("sql__execute__handler", "n_clicks_timestamp"),
        Input("sql__reset__handler", "n_clicks_timestamp"),
    ],
    [State("sql__textarea__input", "value"), State("log_storage", "data")],
)
def execute_query(btn_execute, btn_reset, value, history):

    # reset button clicked
    if btn_execute < btn_reset:
        history_storage = add_to_history(int(btn_reset), "DROP TABLES", history)
        payload = drop_all_tables()
        return payload["status"], None, output_history(history_storage), history_storage

    # execute button clicked
    if value == "" or value is None:
        raise PreventUpdate

    history_storage = add_to_history(int(btn_execute), value, history)
    payload = execute_query(value)

    # no table to show
    if not isinstance(payload["status"], list):
        return payload["status"], None, output_history(history_storage), history_storage

    # if the data is a list/table
    return None, payload["status"], output_history(history_storage), history_storage


def execute_query(statement):
    """ Execute PostgreSQL statement. """

    try:
        results = connection.execute(statement).fetchall()
        headers = connection.execute(statement).keys()

        table = [html.Tr([html.Th([header]) for header in headers])]
        rows = [
            html.Tr([html.Td(item[header]) for header in headers]) for item in results
        ]
        table.extend(rows)

        return {"status": table}

    except Exception as error:
        print(error)
        return {"status": f"{error}"}


def drop_all_tables():
    """ Drops all tables and schema. """

    try:
        first_stmt = connection.execute("DROP SCHEMA public CASCADE;")
        second_stmt = connection.execute("CREATE SCHEMA public;")
        return {"status": f"\n {first_stmt} \n {second_stmt}"}
    except Exception as error:
        return {"status": f"\n {error}"}


def add_to_history(timestamp, statement, history):
    """ Add dictionnary to dcc store. """

    ts = datetime.fromtimestamp(timestamp // 1000.0).strftime("%Y-%m-%d %H:%M:%S")
    if history is None:
        history = []
    history.append(
        {"ts": ts, "statement": statement.replace("\r", "").replace("\n", "")}
    )
    return history


def output_history(history):
    """ Format the history logs to output. """

    output = [f"{item['ts']}: {item['statement']}" for item in reversed(history)]
    return "\n".join(output)


if __name__ == "__main__":
    app.run_server(debug=True)
