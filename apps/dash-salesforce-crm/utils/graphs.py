
from plotly import graph_objs as go
import pandas as pd





### Opportunity page graphs ### ### ### ###




def leads_choropleth_map(status, df):
    if status == "open":
        df = df[
            (df["Status"] == "Open - Not Contacted")
            | (df["Status"] == "Working - Contacted")
        ]

    elif status == "converted":
        df = df[df["Status"] == "Closed - Converted"]

    elif status == "lost":
        df = df[df["Status"] == "Closed - Not Converted"]

    df = df.groupby("State").count()

    scl = [[0.0, "rgb(38, 78, 134)"], [1.0, "#0091D5"]]  # colors scale

    data = [
        dict(
            type="choropleth",
            colorscale=scl,
            locations=df.index,
            z=df["Id"],
            locationmode="USA-states",
            marker=dict(line=dict(color="rgb(255,255,255)", width=2)),
            colorbar=dict(len=0.8),
        )
    ]

    layout = dict(
        autosize=True,
        geo=dict(
            scope="usa",
            projection=dict(type="albers usa"),
            lakecolor="rgb(255, 255, 255)",
        ),
        margin=dict(l=10, r=10, t=0, b=0),
    )
    return dict(data=data, layout=layout)


def lead_source(status, df):
    if status == "open":
        df = df[
            (df["Status"] == "Open - Not Contacted")
            | (df["Status"] == "Working - Contacted")
        ]

    elif status == "converted":
        df = df[df["Status"] == "Closed - Converted"]

    elif status == "lost":
        df = df[df["Status"] == "Closed - Not Converted"]

    nb_leads = len(df.index)
    types = df["LeadSource"].unique().tolist()
    values = []

    # compute % for each leadsource type
    for case_type in types:
        nb_type = df[df["LeadSource"] == case_type].shape[0]
        values.append(nb_type / nb_leads * 100)

    trace = go.Pie(
        labels=types,
        values=values,
        marker={"colors": ["#264e86", "#0074e4", "#74dbef", "#eff0f4"]},
    )

    layout = dict(autosize=True, margin=dict(l=15, r=10, t=0, b=65))
    return dict(data=[trace], layout=layout)


def converted_leads_count(period, df):
    df["CreatedDate"] = pd.to_datetime(df["CreatedDate"], format="%Y-%m-%d")
    df = df[df["Status"] == "Closed - Converted"]

    df = (
        df.groupby([pd.Grouper(key="CreatedDate", freq=period)])
        .count()
        .reset_index()
        .sort_values("CreatedDate")
    )

    trace = go.Scatter(
        x=df["CreatedDate"],
        y=df["Id"],
        name="converted leads",
        fill="tozeroy",
        fillcolor="#e6f2ff",
    )

    data = [trace]

    layout = go.Layout(
        autosize=True,
        xaxis=dict(showgrid=False),
        margin=dict(l=33, r=25, b=37, t=5, pad=4),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return {"data": data, "layout": layout}





### Opportunity page graphs ### ### ### ###




def converted_opportunities(period, source, df):
    df["CreatedDate"] = pd.to_datetime(df["CreatedDate"], format="%Y-%m-%d")

    # source filtering
    if source == "all_s":
        df = df[df["IsWon"] == 1]
    else:
        df = df[(df["LeadSource"] == source) & (df["IsWon"] == 1)]

    # period filtering
    if period == "W-MON":
        df["CreatedDate"] = pd.to_datetime(df["CreatedDate"]) - pd.to_timedelta(
            7, unit="d"
        )
    df = (
        df.groupby([pd.Grouper(key="CreatedDate", freq=period)])
        .count()
        .reset_index()
        .sort_values("CreatedDate")
    )

    # if no results were found
    if df.empty:
        layout = dict(
            autosize=True, annotations=[dict(text="No results found", showarrow=False)]
        )
        return {"data": [], "layout": layout}

    trace = go.Scatter(
        x=df["CreatedDate"],
        y=df["IsWon"],
        name="converted opportunities",
        fill="tozeroy",
        fillcolor="#e6f2ff",
    )

    data = [trace]

    layout = go.Layout(
        autosize=True,
        xaxis=dict(showgrid=False),
        margin=dict(l=35, r=25, b=23, t=5, pad=4),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return {"data": data, "layout": layout}


def opportunities_heat_map_fig(df, stage):
    df = df[pd.notnull(df["Type"])]
    x, y, z = [], df["Type"].unique(), []

    if stage == "all_s":
        x = df["StageName"].unique()
    elif stage == "cold":
        x = ["Needs Analysis", "Prospecting", "Qualification"]
    elif stage == "warm":
        x = ["Value Proposition", "Id. Decision Makers", "Perception Analysis"]
    else:
        x = ["Proposal/Price Quote", "Negotiation/Review", "Closed Won"]

    for lead_type in y:
        z_row = []
        for stage in x:
            probability = df[(df["StageName"] == stage) & (df["Type"] == lead_type)][
                "Probability"
            ].mean()
            z_row.append(probability)
        z.append(z_row)

    trace = dict(
        type="heatmap", z=z, x=x, y=y, name="mean probability", colorscale="Blues"
    )
    layout = dict(
        autosize=True,
        margin=dict(t=25, l=210, b=85, pad=4),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return go.Figure(data=[trace], layout=layout)





### Cases page graphs ### ### ### ###




def cases_pie_chart(df, column, priority, origin, h_orientation=None):
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

    if h_orientation:
        layout["legend"]["orientation"] = "h"

    return {"data": [trace], "layout": layout}


def cases_by_period(df, period, priority):
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


def cases_by_account(cases, salesforce_manager):
    accounts = salesforce_manager.get_accounts()
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
