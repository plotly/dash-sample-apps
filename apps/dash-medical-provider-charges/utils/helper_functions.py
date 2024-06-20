from constants import state_map, data_dict
import dash
from dash.exceptions import PreventUpdate
from dash import dash_table
import pandas as pd


def generate_aggregation(df, metric):
    aggregation = {
        metric[0]: ["min", "mean", "max"],
        metric[1]: ["min", "mean", "max"],
        metric[2]: ["min", "mean", "max"],
    }
    grouped = (
        df.groupby(["Hospital Referral Region (HRR) Description", "Provider Name"])
        .agg(aggregation)
        .reset_index()
    )

    grouped["lat"] = grouped["lon"] = grouped["Provider Street Address"] = grouped[
        "Provider Name"
    ]
    grouped["lat"] = grouped["lat"].apply(lambda x: get_lat_lon_add(df, x)[0])
    grouped["lon"] = grouped["lon"].apply(lambda x: get_lat_lon_add(df, x)[1])
    grouped["Provider Street Address"] = grouped["Provider Street Address"].apply(
        lambda x: get_lat_lon_add(df, x)[2]
    )

    return grouped


def get_lat_lon_add(df, name):
    return [
        df.groupby(["Provider Name"]).get_group(name)["lat"].tolist()[0],
        df.groupby(["Provider Name"]).get_group(name)["lon"].tolist()[0],
        df.groupby(["Provider Name"])
        .get_group(name)["Provider Street Address"]
        .tolist()[0],
    ]


def region_dropdown(select_all, state_select):
    state_raw_data = data_dict[state_select]
    regions = state_raw_data["Hospital Referral Region (HRR) Description"].unique()
    options = [{"label": i, "value": i} for i in regions]

    ctx = dash.callback_context
    if ctx.triggered[0]["prop_id"].split(".")[0] == "region-select-all":
        if select_all == ["All"]:
            value = [i["value"] for i in options]
        else:
            value = dash.no_update
    else:
        value = regions[:4]
    return (
        value,
        options,
        "Medicare Provider Charges in the State of {}".format(state_map[state_select]),
    )


def checklist(selected, select_options, checked):
    if len(selected) < len(select_options) and len(checked) == 0:
        raise PreventUpdate()

    elif len(selected) < len(select_options) and len(checked) == 1:
        return dcc.Checklist(
            id="region-select-all",
            options=[{"label": "Select All Regions", "value": "All"}],
            value=[],
        )

    elif len(selected) == len(select_options) and len(checked) == 1:
        raise PreventUpdate()

    return dcc.Checklist(
        id="region-select-all",
        options=[{"label": "Select All Regions", "value": "All"}],
        value=["All"],
    )


def procedure_stats(procedure_select, geo_select, cost_select, state_select):
    procedure_dict = {
        "DRG": [],
        "Procedure": [],
        "Provider Name": [],
        "Cost Summary": [],
    }

    ctx = dash.callback_context
    prop_id = ""
    if ctx.triggered:
        prop_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if prop_id == "procedure-plot" and procedure_select is not None:
        for point in procedure_select["points"]:
            procedure_dict["DRG"].append(point["y"].split(" - ")[0])
            procedure_dict["Procedure"].append(point["y"].split(" - ")[1])

            procedure_dict["Provider Name"].append(point["customdata"])
            procedure_dict["Cost Summary"].append(("${:,.2f}".format(point["x"])))

    # Display all procedures at selected hospital
    provider_select = []

    if prop_id == "geo-map" and geo_select is not None:
        for point in geo_select["points"]:
            provider = point["customdata"][0]
            provider_select.append(provider)

        state_raw_data = data_dict[state_select]
        provider_filtered = state_raw_data[
            state_raw_data["Provider Name"].isin(provider_select)
        ]

        for i in range(len(provider_filtered)):
            procedure_dict["DRG"].append(
                provider_filtered.iloc[i]["DRG Definition"].split(" - ")[0]
            )
            procedure_dict["Procedure"].append(
                provider_filtered.iloc[i]["DRG Definition"].split(" - ")[1]
            )
            procedure_dict["Provider Name"].append(
                provider_filtered.iloc[i]["Provider Name"]
            )
            procedure_dict["Cost Summary"].append(
                "${:,.2f}".format(provider_filtered.iloc[0][cost_select])
            )

    procedure_data_df = pd.DataFrame(data=procedure_dict)

    return dash_table.DataTable(
        id="procedure-stats-table",
        columns=[{"name": i, "id": i} for i in procedure_dict.keys()],
        data=procedure_data_df.to_dict(orient="records"),
        filter_action="native",
        sort_action="native",
        style_cell={
            "textOverflow": "ellipsis",
            "background-color": "#242a3b",
            "color": "#7b7d8d",
        },
        sort_mode="multi",
        page_size=5,
        style_as_list_view=False,
        style_header={"background-color": "#1f2536", "padding": "2px 12px 0px 12px"},
    )
