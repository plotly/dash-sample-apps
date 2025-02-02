import plotly.graph_objs as go
from constants import mapbox_access_token
import dash
from dash import dash_table
import pandas as pd
from utils.helper_functions import generate_aggregation
from constants import data_dict, cost_metric


def generate_geo_map(geo_data, selected_metric, region_select, procedure_select):
    filtered_data = geo_data[
        geo_data["Hospital Referral Region (HRR) Description"].isin(region_select)
    ]

    colors = ["#21c7ef", "#76f2ff", "#ff6969", "#ff1717"]

    hospitals = []

    lat = filtered_data["lat"].tolist()
    lon = filtered_data["lon"].tolist()
    average_covered_charges_mean = filtered_data[selected_metric]["mean"].tolist()
    regions = filtered_data["Hospital Referral Region (HRR) Description"].tolist()
    provider_name = filtered_data["Provider Name"].tolist()

    # Cost metric mapping from aggregated data

    cost_metric_data = {}
    cost_metric_data["min"] = filtered_data[selected_metric]["mean"].min()
    cost_metric_data["max"] = filtered_data[selected_metric]["mean"].max()
    cost_metric_data["mid"] = (cost_metric_data["min"] + cost_metric_data["max"]) / 2
    cost_metric_data["low_mid"] = (
        cost_metric_data["min"] + cost_metric_data["mid"]
    ) / 2
    cost_metric_data["high_mid"] = (
        cost_metric_data["mid"] + cost_metric_data["max"]
    ) / 2

    for i in range(len(lat)):
        val = average_covered_charges_mean[i]
        region = regions[i]
        provider = provider_name[i]

        if val <= cost_metric_data["low_mid"]:
            color = colors[0]
        elif cost_metric_data["low_mid"] < val <= cost_metric_data["mid"]:
            color = colors[1]
        elif cost_metric_data["mid"] < val <= cost_metric_data["high_mid"]:
            color = colors[2]
        else:
            color = colors[3]

        selected_index = []
        if provider in procedure_select["hospital"]:
            selected_index = [0]

        hospital = go.Scattermapbox(
            lat=[lat[i]],
            lon=[lon[i]],
            mode="markers",
            marker=dict(
                color=color,
                showscale=True,
                colorscale=[
                    [0, "#21c7ef"],
                    [0.33, "#76f2ff"],
                    [0.66, "#ff6969"],
                    [1, "#ff1717"],
                ],
                cmin=cost_metric_data["min"],
                cmax=cost_metric_data["max"],
                size=10
                * (1 + (val + cost_metric_data["min"]) / cost_metric_data["mid"]),
                colorbar=dict(
                    x=0.9,
                    len=0.7,
                    title=dict(
                        text="Average Cost",
                        font={"color": "#737a8d", "family": "Open Sans"},
                    ),
                    titleside="top",
                    tickmode="array",
                    tickvals=[cost_metric_data["min"], cost_metric_data["max"]],
                    ticktext=[
                        "${:,.2f}".format(cost_metric_data["min"]),
                        "${:,.2f}".format(cost_metric_data["max"]),
                    ],
                    ticks="outside",
                    thickness=15,
                    tickfont={"family": "Open Sans", "color": "#737a8d"},
                ),
            ),
            opacity=0.8,
            selectedpoints=selected_index,
            selected=dict(marker={"color": "#ffff00"}),
            customdata=[(provider, region)],
            hoverinfo="text",
            text=provider
            + "<br>"
            + region
            + "<br>Average Procedure Cost:"
            + " ${:,.2f}".format(val),
        )
        hospitals.append(hospital)

    layout = go.Layout(
        margin=dict(l=10, r=10, t=20, b=10, pad=5),
        plot_bgcolor="#171b26",
        paper_bgcolor="#171b26",
        clickmode="event+select",
        hovermode="closest",
        showlegend=False,
        mapbox=go.layout.Mapbox(
            accesstoken=mapbox_access_token,
            bearing=10,
            center=go.layout.mapbox.Center(
                lat=filtered_data.lat.mean(), lon=filtered_data.lon.mean()
            ),
            pitch=5,
            zoom=5,
            style="mapbox://styles/plotlymapbox/cjvppq1jl1ips1co3j12b9hex",
        ),
    )

    return {"data": hospitals, "layout": layout}


def generate_procedure_plot(raw_data, cost_select, region_select, provider_select):
    procedure_data = raw_data[
        raw_data["Hospital Referral Region (HRR) Description"].isin(region_select)
    ].reset_index()

    traces = []
    selected_index = procedure_data[
        procedure_data["Provider Name"].isin(provider_select)
    ].index

    text = (
        procedure_data["Provider Name"]
        + "<br>"
        + "<b>"
        + procedure_data["DRG Definition"].map(str)
        + "/<b> <br>"
        + "Average Procedure Cost: $ "
        + procedure_data[cost_select].map(str)
    )

    provider_trace = go.Box(
        y=procedure_data["DRG Definition"],
        x=procedure_data[cost_select],
        name="",
        customdata=procedure_data["Provider Name"],
        boxpoints="all",
        jitter=0,
        pointpos=0,
        hoveron="points",
        fillcolor="rgba(0,0,0,0)",
        line=dict(color="rgba(0,0,0,0)"),
        hoverinfo="text",
        hovertext=text,
        selectedpoints=selected_index,
        selected=dict(marker={"color": "#FFFF00", "size": 13}),
        unselected=dict(marker={"opacity": 0.2}),
        marker=dict(
            line=dict(width=1, color="#000000"),
            color="#21c7ef",
            opacity=0.7,
            symbol="square",
            size=12,
        ),
    )

    traces.append(provider_trace)

    layout = go.Layout(
        showlegend=False,
        hovermode="closest",
        dragmode="select",
        clickmode="event+select",
        xaxis=dict(
            zeroline=False,
            automargin=True,
            showticklabels=True,
            title=dict(text="Procedure Cost", font=dict(color="#737a8d")),
            linecolor="#737a8d",
            tickfont=dict(color="#737a8d"),
            type="log",
        ),
        yaxis=dict(
            automargin=True,
            showticklabels=True,
            tickfont=dict(color="#737a8d"),
            gridcolor="#171b26",
        ),
        plot_bgcolor="#171b26",
        paper_bgcolor="#171b26",
    )
    # x : procedure, y: cost,
    return {"data": traces, "layout": layout}


def hospital_datatable(geo_select, procedure_select, cost_select, state_select):
    state_agg = generate_aggregation(data_dict[state_select], cost_metric)
    # make table from geo-select
    geo_data_dict = {
        "Provider Name": [],
        "City": [],
        "Street Address": [],
        "Maximum Cost ($)": [],
        "Minimum Cost ($)": [],
    }

    ctx = dash.callback_context
    if ctx.triggered:
        prop_id = ctx.triggered[0]["prop_id"].split(".")[0]

        # make table from procedure-select
        if prop_id == "procedure-plot" and procedure_select is not None:

            for point in procedure_select["points"]:
                provider = point["customdata"]

                dff = state_agg[state_agg["Provider Name"] == provider]

                geo_data_dict["Provider Name"].append(point["customdata"])
                city = dff["Hospital Referral Region (HRR) Description"].tolist()[0]
                geo_data_dict["City"].append(city)

                address = dff["Provider Street Address"].tolist()[0]
                geo_data_dict["Street Address"].append(address)

                geo_data_dict["Maximum Cost ($)"].append(
                    dff[cost_select]["max"].tolist()[0]
                )
                geo_data_dict["Minimum Cost ($)"].append(
                    dff[cost_select]["min"].tolist()[0]
                )

        if prop_id == "geo-map" and geo_select is not None:

            for point in geo_select["points"]:
                provider = point["customdata"][0]
                dff = state_agg[state_agg["Provider Name"] == provider]

                geo_data_dict["Provider Name"].append(point["customdata"][0])
                geo_data_dict["City"].append(point["customdata"][1].split("- ")[1])

                address = dff["Provider Street Address"].tolist()[0]
                geo_data_dict["Street Address"].append(address)

                geo_data_dict["Maximum Cost ($)"].append(
                    dff[cost_select]["max"].tolist()[0]
                )
                geo_data_dict["Minimum Cost ($)"].append(
                    dff[cost_select]["min"].tolist()[0]
                )

        geo_data_df = pd.DataFrame(data=geo_data_dict)
        data = geo_data_df.to_dict(orient="records")

    else:
        data = [{}]

    return dash_table.DataTable(
        id="cost-stats-table",
        columns=[{"name": i, "id": i} for i in geo_data_dict.keys()],
        data=data,
        filter_action="native",
        page_size=5,
        style_cell={"background-color": "#242a3b", "color": "#7b7d8d"},
        style_as_list_view=False,
        style_header={"background-color": "#1f2536", "padding": "0px 5px"},
    )


def update_geo_map(cost_select, region_select, procedure_select, state_select):
    # generate geo map from state-select, procedure-select
    state_agg_data = generate_aggregation(data_dict[state_select], cost_metric)

    provider_data = {"procedure": [], "hospital": []}
    if procedure_select is not None:
        for point in procedure_select["points"]:
            provider_data["procedure"].append(point["y"])
            provider_data["hospital"].append(point["customdata"])

    return generate_geo_map(state_agg_data, cost_select, region_select, provider_data)


def update_procedure_plot(cost_select, region_select, geo_select, state_select):
    # generate procedure plot from selected provider
    state_raw_data = data_dict[state_select]

    provider_select = []
    if geo_select is not None:
        for point in geo_select["points"]:
            provider_select.append(point["customdata"][0])
    return generate_procedure_plot(
        state_raw_data, cost_select, region_select, provider_select
    )
