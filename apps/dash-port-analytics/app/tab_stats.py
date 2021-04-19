import numpy as np
import pandas as pd
import plotly.graph_objects as go
from app import helpers
from config import strings, styles


def get_stats_card1_data(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
):
    """
    Gets values for the first card in the Stats tab.

    :param df: Pandas DataFrame, input data
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: list - [pct, direction] if no errors, [0, '-'] if there are errors
    """
    data = helpers.filter_by_vessel_and_time(
        df=df, vessel_type=vessel_type, year=year, month=month
    )
    df_port = data[data["port"] == port]
    df_other = data[data["port"] != port]

    try:
        pct = -np.round(
            (
                (df_other["n"].mean() - df_port["n"].mean())
                / np.abs(df_other["n"].mean())
            )
            * 100
        )
        pct = int(pct)
        direction = "lower" if pct < 0 else "higher"
        return [np.abs(pct), direction]
    except Exception as _:
        return [0, "-"]


def get_stats_card2_data(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
):
    """
    Gets values for the second card in the Stats tab.

    :param df: Pandas DataFrame, input data
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: list - [pct, direction] if no errors, [0, '-'] if there are errors
    """
    data = helpers.filter_by_vessel_and_time(
        df=df, vessel_type=vessel_type, year=year, month=month
    )
    data = data[data["len_stop"] > 20]
    data = data.groupby(by=["port", "ship_type"]).mean()
    data = data.reset_index()
    port_stop_mean = data[data["port"] == port]["len_stop"].mean()
    port_other_mean = data[data["port"] != port]["len_stop"].mean()

    try:
        pct = -int(
            np.round(((port_other_mean - port_stop_mean) / port_other_mean) * 100)
        )
        direction = "shorter" if pct < 0 else "longer"
        return [np.abs(pct), direction]
    except Exception as _:
        return [0, "-"]


def get_stats_card3_data(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
):
    """
    Gets values for the third card in the Stats tab.

    :param df: Pandas DataFrame, input data
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: list - [pct, direction] if no errors, [0, '-'] if there are errors
    """
    data = helpers.filter_by_vessel_and_time(
        df=df, vessel_type=vessel_type, year=year, month=month
    )
    df_port = data[data["port"] == port]
    df_other = data[data["port"] != port]

    try:
        pct = -np.round(
            (
                (df_other["sum_dwt"].mean() - df_port["sum_dwt"].mean())
                / np.abs(df_other["sum_dwt"].mean())
            )
            * 100
        )
        pct = int(pct)
        direction = "lower" if pct < 0 else "higher"
        return [np.abs(pct), direction]
    except Exception as _:
        return [0, 0]


def plot_stats_total_num_vessels(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
) -> go.Figure:
    """
    Returns a figure for the first chart on the Stats tab. It shows the total number of vessels in port
    by applied conditions.

    :param df: Pandas DataFrame, input data
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: Plotly figure
    """
    data = helpers.filter_by_port_vessel_and_time(
        df=df, port=port, vessel_type=vessel_type, year=year, month=month
    )
    if len(data) > 0:
        plot_data = []
        for dt in data["date"].unique():
            for stype in data["ship_type"].unique():
                curr = data[(data["date"] == dt) & (data["ship_type"] == stype)]
                if len(curr) > 0:
                    plot_data.append(
                        {"date": dt, "ship_type": stype, "num": curr["n"].values[0]}
                    )

        plot_data = pd.DataFrame(plot_data)
        plot_data["color"] = plot_data["ship_type"].apply(helpers.generate_color)
        fig_data = []
        for stype in plot_data["ship_type"].unique():
            ss = plot_data[plot_data["ship_type"] == stype]
            fig_data.append(
                go.Bar(
                    name=stype,
                    x=ss["date"].tolist(),
                    y=ss["num"].tolist(),
                    marker_color=ss.iloc[0]["color"],
                )
            )
    else:
        fig_data = go.Bar(x=np.arange(1, 9), y=[0] * 8)

    return go.Figure(
        data=fig_data,
        layout=styles.generate_plot_layout(
            x_title=strings.CHART_STATS_TOTAL_VESSELS_X,
            y_title=strings.CHART_STATS_TOTAL_VESSELS_Y,
            bar_mode="stack",
        ),
    )


def plot_avg_vessel_stop_duration(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
) -> go.Figure:
    """
    Returns a figure for the second chart on the Stats tab. It shows the average stop duration
    by applied conditions.

    :param df: Pandas DataFrame, input data
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: Plotly figure
    """
    data = helpers.filter_by_port_vessel_and_time(
        df=df, port=port, vessel_type=vessel_type, year=year, month=month
    )
    if len(data) > 0:
        data = data.groupby(by="ship_type").mean().reset_index()
        data = data[["ship_type", "len_stop"]]
        data["len_stop"] = data["len_stop"].apply(lambda x: np.round(x, 2))
        data["color"] = data["ship_type"].apply(helpers.generate_color)
        fig_data = go.Bar(
            x=data["ship_type"], y=data["len_stop"], marker_color=data["color"]
        )
    else:
        fig_data = go.Bar(x=np.arange(1, 9), y=[0] * 8)

    return go.Figure(
        data=fig_data,
        layout=styles.generate_plot_layout(
            x_title=strings.CHART_STATS_STOP_DUR_X,
            y_title=strings.CHART_STATS_STOP_DUR_Y,
            bar_mode="stack",
        ),
    )


def plot_total_capacity_of_vessels(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
) -> go.Figure:
    """
    Returns a figure for the third chart on the Stats tab. It shows the total capacity of vessels
    by applied conditions.

    :param df: Pandas DataFrame, input data
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: Plotly figure
    """
    data = helpers.filter_by_port_vessel_and_time(
        df=df, port=port, vessel_type=vessel_type, year=year, month=month
    )
    if len(data) > 0:
        fig_data = []
        data = data.groupby(by=["date", "ship_type"]).sum().reset_index()
        data = data[["date", "ship_type", "sum_dwt"]]
        data["color"] = data["ship_type"].apply(helpers.generate_color)

        for stype in data["ship_type"].unique():
            ss = data[data["ship_type"] == stype]
            fig_data.append(
                go.Bar(
                    name=stype,
                    x=ss["date"].tolist(),
                    y=ss["sum_dwt"].tolist(),
                    marker_color=ss.iloc[0]["color"],
                )
            )
    else:
        fig_data = go.Bar(x=np.arange(1, 9), y=[0] * 8)

    return go.Figure(
        data=fig_data,
        layout=styles.generate_plot_layout(
            x_title=strings.CHART_STATS_TOTAL_CAP_VESSELS_X,
            y_title=strings.CHART_STATS_TOTAL_CAP_VESSELS_Y,
            bar_mode="stack",
        ),
    )
