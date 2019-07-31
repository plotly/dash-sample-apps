# -*- coding: utf-8 -*-
import pandas as pd
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
from plotly import graph_objs as go

from app import app, indicator, sf_manager

colors = {"background": "#F3F6FA", "background_div": "white"}

accounts = sf_manager.get_accounts()
contacts = sf_manager.get_contacts()
users = sf_manager.get_users()

# returns pie chart based on filters values
# column makes the function reusable


def pie_chart(df, column, priority, origin):
    df = df.dropna(subset=["Type", "Reason", "Origin"])
    nb_cases = len(df.index)
    types = []
    values = []

    # filter priority and origin
    if priority == "all_p":
        if origin == "all":
            types = df[column].unique().tolist()
        else:
            types = df[df["Origin"] == origin][column].unique().tolist()
    else:
        if origin == "all":
            types = df[df["Priority"] == priority][column].unique().tolist()
        else:
            types = (
                df[(df["Priority"] == priority) & (df["Origin"] == origin)][column]
                .unique()
                .tolist()
            )

    # if no results were found
    if types == []:
        layout = dict(
            autosize=True, annotations=[dict(text="No results found", showarrow=False)]
        )
        return {"data": [], "layout": layout}

    for case_type in types:
        nb_type = df.loc[df[column] == case_type].shape[0]
        values.append(nb_type / nb_cases * 100)

    layout = go.Layout(
        autosize=True,
        margin=dict(l=0, r=0, b=0, t=4, pad=8),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    trace = go.Pie(
        labels=types,
        values=values,
        marker={"colors": ["#264e86", "#0074e4", "#74dbef", "#eff0f4"]},
    )

    return {"data": [trace], "layout": layout}


def cases_by_period(df, period, priority, origin):
    df = df.dropna(subset=["Type", "Reason", "Origin"])
    stages = df["Type"].unique()

    # priority filtering
    if priority != "all_p":
        df = df[df["Priority"] == priority]

    # period filtering
    df["CreatedDate"] = pd.to_datetime(df["CreatedDate"], format="%Y-%m-%d")
    if period == "W-MON":
        df["CreatedDate"] = pd.to_datetime(df["CreatedDate"]) - pd.to_timedelta(
            7, unit="d"
        )
    df = df.groupby([pd.Grouper(key="CreatedDate", freq=period), "Type"]).count()

    dates = df.index.get_level_values("CreatedDate").unique()
    dates = [str(i) for i in dates]

    co = {  # colors for stages
        "Electrical": "#264e86",
        "Other": "#0074e4",
        "Structural": "#74dbef",
        "Mechanical": "#eff0f4",
        "Electronic": "rgb(255, 127, 14)",
    }

    data = []
    for stage in stages:
        stage_rows = []
        for date in dates:
            try:
                row = df.loc[(date, stage)]
                stage_rows.append(row["IsDeleted"])
            except Exception as e:
                stage_rows.append(0)

        data_trace = go.Bar(
            x=dates, y=stage_rows, name=stage, marker=dict(color=co[stage])
        )
        data.append(data_trace)

    layout = go.Layout(
        autosize=True,
        barmode="stack",
        margin=dict(l=40, r=25, b=40, t=0, pad=4),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return {"data": data, "layout": layout}


def cases_by_account(cases):
    cases = cases.dropna(subset=["AccountId"])
    cases = pd.merge(cases, accounts, left_on="AccountId", right_on="Id")
    cases = cases.groupby(["AccountId", "Name"]).count()
    cases = cases.sort_values("IsDeleted")
    data = [
        go.Bar(
            y=cases.index.get_level_values("Name"),
            x=cases["IsDeleted"],
            orientation="h",
            marker=dict(color="#0073e4"),
        )
    ]  # x could be any column value since its a count

    layout = go.Layout(
        autosize=True,
        barmode="stack",
        margin=dict(l=210, r=25, b=20, t=0, pad=4),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return {"data": data, "layout": layout}


# returns modal (hidden by default)
def modal():
    contacts["Name"] = (
        contacts["Salutation"]
        + " "
        + contacts["FirstName"]
        + " "
        + contacts["LastName"]
    )
    return html.Div(
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.Span(
                                    "New Case",
                                    style={
                                        "color": "#506784",
                                        "fontWeight": "bold",
                                        "fontSize": "20",
                                    },
                                ),
                                html.Span(
                                    "×",
                                    id="cases_modal_close",
                                    n_clicks=0,
                                    style={
                                        "float": "right",
                                        "cursor": "pointer",
                                        "marginTop": "0",
                                        "marginBottom": "17",
                                    },
                                ),
                            ],
                            className="row",
                            style={"borderBottom": "1px solid #C8D4E3"},
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.P(
                                            "Account name",
                                            style={
                                                "textAlign": "left",
                                                "marginBottom": "2",
                                                "marginTop": "4",
                                            },
                                        ),
                                        html.Div(
                                            dcc.Dropdown(
                                                id="new_case_account",
                                                options=[
                                                    {
                                                        "label": row["Name"],
                                                        "value": row["Id"],
                                                    }
                                                    for index, row in accounts.iterrows()
                                                ],
                                                clearable=False,
                                                value=accounts.iloc[0].Id,
                                            )
                                        ),
                                        html.P(
                                            "Priority",
                                            style={
                                                "textAlign": "left",
                                                "marginBottom": "2",
                                                "marginTop": "4",
                                            },
                                        ),
                                        dcc.Dropdown(
                                            id="new_case_priority",
                                            options=[
                                                {"label": "High", "value": "High"},
                                                {"label": "Medium", "value": "Medium"},
                                                {"label": "Low", "value": "Low"},
                                            ],
                                            value="Medium",
                                            clearable=False,
                                        ),
                                        html.P(
                                            "Origin",
                                            style={
                                                "textAlign": "left",
                                                "marginBottom": "2",
                                                "marginTop": "4",
                                            },
                                        ),
                                        dcc.Dropdown(
                                            id="new_case_origin",
                                            options=[
                                                {"label": "Phone", "value": "Phone"},
                                                {"label": "Web", "value": "Web"},
                                                {"label": "Email", "value": "Email"},
                                            ],
                                            value="Phone",
                                            clearable=False,
                                        ),
                                        html.P(
                                            "Reason",
                                            style={
                                                "textAlign": "left",
                                                "marginBottom": "2",
                                                "marginTop": "4",
                                            },
                                        ),
                                        dcc.Dropdown(
                                            id="new_case_reason",
                                            options=[
                                                {
                                                    "label": "Installation",
                                                    "value": "Installation",
                                                },
                                                {
                                                    "label": "Equipment Complexity",
                                                    "value": "Equipment Complexity",
                                                },
                                                {
                                                    "label": "Performance",
                                                    "value": "Performance",
                                                },
                                                {
                                                    "label": "Breakdown",
                                                    "value": "Breakdown",
                                                },
                                                {
                                                    "label": "Equipment Design",
                                                    "value": "Equipment Design",
                                                },
                                                {
                                                    "label": "Feedback",
                                                    "value": "Feedback",
                                                },
                                                {"label": "Other", "value": "Other"},
                                            ],
                                            value="Installation",
                                            clearable=False,
                                        ),
                                        html.P(
                                            "Subject",
                                            style={
                                                "float": "left",
                                                "marginTop": "4",
                                                "marginBottom": "2",
                                            },
                                            className="row",
                                        ),
                                        dcc.Input(
                                            id="new_case_subject",
                                            placeholder="The Subject of the case",
                                            type="text",
                                            value="",
                                            style={"width": "100%"},
                                        ),
                                    ],
                                    className="six columns",
                                    style={"paddingRight": "15"},
                                ),
                                html.Div(
                                    [
                                        html.P(
                                            "Contact name",
                                            style={
                                                "textAlign": "left",
                                                "marginBottom": "2",
                                                "marginTop": "4",
                                            },
                                        ),
                                        html.Div(
                                            dcc.Dropdown(
                                                id="new_case_contact",
                                                options=[
                                                    {
                                                        "label": row["Name"],
                                                        "value": row["Id"],
                                                    }
                                                    for index, row in contacts.iterrows()
                                                ],
                                                clearable=False,
                                                value=contacts.iloc[0].Id,
                                            )
                                        ),
                                        html.P(
                                            "Type",
                                            style={
                                                "textAlign": "left",
                                                "marginBottom": "2",
                                                "marginTop": "4",
                                            },
                                        ),
                                        dcc.Dropdown(
                                            id="new_case_type",
                                            options=[
                                                {
                                                    "label": "Electrical",
                                                    "value": "Electrical",
                                                },
                                                {
                                                    "label": "Mechanical",
                                                    "value": "Mechanical",
                                                },
                                                {
                                                    "label": "Electronic",
                                                    "value": "Electronic",
                                                },
                                                {
                                                    "label": "Structural",
                                                    "value": "Structural",
                                                },
                                                {"label": "Other", "value": "Other"},
                                            ],
                                            value="Electrical",
                                        ),
                                        html.P(
                                            "Status",
                                            style={
                                                "textAlign": "left",
                                                "marginBottom": "2",
                                                "marginTop": "4",
                                            },
                                        ),
                                        dcc.Dropdown(
                                            id="new_case_status",
                                            options=[
                                                {"label": "New", "value": "New"},
                                                {
                                                    "label": "Working",
                                                    "value": "Working",
                                                },
                                                {
                                                    "label": "Escalated",
                                                    "value": "Escalated",
                                                },
                                                {"label": "Closed", "value": "Closed"},
                                            ],
                                            value="New",
                                        ),
                                        html.P(
                                            "Supplied Email",
                                            style={
                                                "textAlign": "left",
                                                "marginBottom": "2",
                                                "marginTop": "4",
                                            },
                                        ),
                                        dcc.Input(
                                            id="new_case_email",
                                            placeholder="email",
                                            type="email",
                                            value="",
                                            style={"width": "100%"},
                                        ),
                                        html.P(
                                            "Description",
                                            style={
                                                "float": "left",
                                                "marginTop": "4",
                                                "marginBottom": "2",
                                            },
                                            className="row",
                                        ),
                                        dcc.Textarea(
                                            id="new_case_description",
                                            placeholder="Description of the case",
                                            value="",
                                            style={"width": "100%"},
                                        ),
                                    ],
                                    className="six columns",
                                    style={"paddingLeft": "15"},
                                ),
                            ],
                            style={"marginTop": "10", "textAlign": "center"},
                            className="row",
                        ),
                        html.Span(
                            "Submit",
                            id="submit_new_case",
                            n_clicks=0,
                            className="button button--primary add pretty_container",
                        ),
                    ],
                    className="modal-content",
                    style={"textAlign": "center", "border": "1px solid #C8D4E3"},
                )
            ],
            className="modal",
        ),
        id="cases_modal",
        style={"display": "none"},
    )


layout = [
    html.Div(
        id="cases_grid",
        children=[
            html.Div(
                className="control dropdown-styles",
                children=dcc.Dropdown(
                    id="cases_period_dropdown",
                    options=[
                        {"label": "By day", "value": "D"},
                        {"label": "By week", "value": "W-MON"},
                        {"label": "By month", "value": "M"},
                    ],
                    value="D",
                    clearable=False,
                ),
            ),
            html.Div(
                className="control dropdown-styles",
                children=dcc.Dropdown(
                    id="priority_dropdown",
                    options=[
                        {"label": "All priority", "value": "all_p"},
                        {"label": "High priority", "value": "High"},
                        {"label": "Medium priority", "value": "Medium"},
                        {"label": "Low priority", "value": "Low"},
                    ],
                    value="all_p",
                    clearable=False,
                ),
            ),
            html.Div(
                className="control dropdown-styles",
                children=dcc.Dropdown(
                    id="origin_dropdown",
                    options=[
                        {"label": "All origins", "value": "all"},
                        {"label": "Phone", "value": "Phone"},
                        {"label": "Web", "value": "Web"},
                        {"label": "Email", "value": "Email"},
                    ],
                    value="all",
                    clearable=False,
                ),
            ),
            html.Span(
                "Add new",
                id="new_case",
                n_clicks=0,
                className="button button--primary add pretty_container",
            ),
            html.Div(
                id="cases_indicators",
                className="row indicators",
                children=[
                    indicator("#00cc96", "Low priority cases", "left_cases_indicator"),
                    indicator(
                        "#119DFF", "Medium priority cases", "middle_cases_indicator"
                    ),
                    indicator(
                        "#EF553B", "High priority cases", "right_cases_indicator"
                    ),
                ],
            ),
            html.Div(
                id="cases_types_container",
                className="pretty_container chart_div",
                children=[
                    html.P("Cases Type"),
                    dcc.Graph(
                        id="cases_types",
                        config=dict(displayModeBar=False),
                        style={"height": "89%", "width": "98%"},
                    ),
                ],
            ),
            html.Div(
                id="cases_reasons_container",
                className="chart_div pretty_container",
                children=[
                    html.P("Cases Reasons"),
                    dcc.Graph(id="cases_reasons", config=dict(displayModeBar=False)),
                ],
            ),
            html.Div(
                id="cases_by_period_container",
                className="pretty_container chart_div",
                children=[
                    html.P("Cases over Time"),
                    dcc.Graph(id="cases_by_period", config=dict(displayModeBar=False)),
                ],
            ),
            html.Div(
                id="cases_by_account_container",
                className="pretty_container chart_div",
                children=[
                    html.P("Cases by Company"),
                    dcc.Graph(id="cases_by_account", config=dict(displayModeBar=False)),
                ],
            ),
        ],
    ),
    modal(),
]


@app.callback(Output("left_cases_indicator", "children"), [Input("cases_df", "data")])
def left_cases_indicator_callback(df):
    df = pd.read_json(df, orient="split")
    low = len(df[(df["Priority"] == "Low") & (df["Status"] == "New")]["Priority"].index)
    return dcc.Markdown("**{}**".format(low))


@app.callback(Output("middle_cases_indicator", "children"), [Input("cases_df", "data")])
def middle_cases_indicator_callback(df):
    df = pd.read_json(df, orient="split")
    medium = len(
        df[(df["Priority"] == "Medium") & (df["Status"] == "New")]["Priority"].index
    )
    return dcc.Markdown("**{}**".format(medium))


@app.callback(Output("right_cases_indicator", "children"), [Input("cases_df", "data")])
def right_cases_indicator_callback(df):
    df = pd.read_json(df, orient="split")
    high = len(
        df[(df["Priority"] == "High") & (df["Status"] == "New")]["Priority"].index
    )
    return dcc.Markdown("**{}**".format(high))


@app.callback(
    Output("cases_reasons", "figure"),
    [
        Input("priority_dropdown", "value"),
        Input("origin_dropdown", "value"),
        Input("cases_df", "data"),
    ],
)
def cases_reasons_callback(priority, origin, df):
    df = pd.read_json(df, orient="split")
    chart = pie_chart(df, "Reason", priority, origin)
    return chart


@app.callback(
    Output("cases_types", "figure"),
    [
        Input("priority_dropdown", "value"),
        Input("origin_dropdown", "value"),
        Input("cases_df", "data"),
    ],
)
def cases_types_callback(priority, origin, df):
    df = pd.read_json(df, orient="split")
    chart = pie_chart(df, "Type", priority, origin)
    chart["layout"]["legend"]["orientation"] = "h"
    return chart


@app.callback(
    Output("cases_by_period", "figure"),
    [
        Input("cases_period_dropdown", "value"),
        Input("origin_dropdown", "value"),
        Input("priority_dropdown", "value"),
        Input("cases_df", "data"),
    ],
)
def cases_period_callback(period, origin, priority, df):
    df = pd.read_json(df, orient="split")
    return cases_by_period(df, period, priority, origin)


@app.callback(Output("cases_by_account", "figure"), [Input("cases_df", "data")])
def cases_account_callback(df):
    df = pd.read_json(df, orient="split")
    return cases_by_account(df)


@app.callback(Output("cases_modal", "style"), [Input("new_case", "n_clicks")])
def display_cases_modal_callback(n):
    if n > 0:
        return {"display": "block"}
    return {"display": "none"}


@app.callback(
    Output("new_case", "n_clicks"),
    [Input("cases_modal_close", "n_clicks"), Input("submit_new_case", "n_clicks")],
)
def close_modal_callback(n, n2):
    return 0


@app.callback(
    Output("cases_df", "data"),
    [Input("submit_new_case", "n_clicks")],
    [
        State("new_case_account", "value"),
        State("new_case_origin", "value"),
        State("new_case_reason", "value"),
        State("new_case_subject", "value"),
        State("new_case_contact", "value"),
        State("new_case_type", "value"),
        State("new_case_status", "value"),
        State("new_case_description", "value"),
        State("new_case_priority", "value"),
        State("cases_df", "data"),
    ],
)
def add_case_callback(
    n_clicks,
    account_id,
    origin,
    reason,
    subject,
    contact_id,
    case_type,
    status,
    description,
    priority,
    current_df,
):
    if n_clicks > 0:
        query = {
            "AccountId": account_id,
            "Origin": origin,
            "Reason": reason,
            "Subject": subject,
            "ContactId": contact_id,
            "Type": case_type,
            "Status": status,
            "Description": description,
            "Priority": priority,
        }

        sf_manager.add_case(query)
        df = sf_manager.get_cases()
        return df.to_json(orient="split")

    return current_df
