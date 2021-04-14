import pandas as pd
import plotly.graph_objects as go
from scipy import signal
from app import helpers
from config import styles, strings


def plot_num_vessels_comparison(
    df: pd.DataFrame, port1: str, port2: str, vessel_type: str
) -> go.Figure:
    """
    Returns a figure for the first chart on the Compare tab. It shows per day comparison between
    number of vessels in port by applied conditions.

    :param df: Pandas DataFrame, input data
    :param port1: str, a port to compare
    :param port2: str, a port to compare
    :param vessel_type: str, vessel type of interest
    :return: Plotly figure
    """
    data = helpers.filter_by_vessel_and_port(
        df=df, port1=port1, port2=port2, vessel_type=vessel_type
    )
    data_port_1 = (
        data[data["port"] == port1]
        .groupby(by="date")
        .sum()
        .reset_index()[["date", "n"]]
    )
    data_port_2 = (
        data[data["port"] == port2]
        .groupby(by="date")
        .sum()
        .reset_index()[["date", "n"]]
    )
    fig_data = [
        go.Bar(
            x=data_port_1["date"].tolist(),
            y=data_port_1["n"].tolist(),
            name=port1,
            marker_color=styles.COLOR_APPSILON_1,
        ),
        go.Bar(
            x=data_port_2["date"].tolist(),
            y=data_port_2["n"].tolist(),
            name=port2,
            marker_color=styles.COLOR_APPSILON_8,
        ),
    ]
    return go.Figure(
        data=fig_data,
        layout=styles.generate_plot_layout(
            x_title=strings.CHART_COMPARE_NUM_VESSELS_X,
            y_title=strings.CHART_COMPARE_NUM_VESSELS_Y,
            bar_mode="group",
        ),
    )


def plot_avg_stop_duration_comparison(
    df: pd.DataFrame, port1: str, port2: str, vessel_type: str
) -> go.Figure:
    """
    Returns a figure for the second chart on the Compare tab. It shows per day comparison between
    average stop duration by applied conditions.

    :param df: Pandas DataFrame, input data
    :param port1: str, a port to compare
    :param port2: str, a port to compare
    :param vessel_type: str, vessel type of interest
    :return: Plotly figure
    """
    data = helpers.filter_by_vessel_and_port(
        df=df, port1=port1, port2=port2, vessel_type=vessel_type
    )
    xs = [port1, port2]
    ys = [
        data[data["port"] == port1]["len_stop"].mean(),
        data[data["port"] == port2]["len_stop"].mean(),
    ]

    return go.Figure(
        data=[
            go.Bar(
                x=xs,
                y=ys,
                marker_color=[styles.COLOR_APPSILON_1, styles.COLOR_APPSILON_8],
            )
        ],
        layout=styles.generate_plot_layout(
            x_title=strings.CHART_COMPARE_AVG_DURATION_X,
            y_title=strings.CHART_COMPARE_AVG_DURATION_Y,
            bar_mode="",
        ),
    )


def plot_avg_sum_capacity_comparison(
    df: pd.DataFrame, port1: str, port2: str, vessel_type: str
) -> go.Figure:
    """
    Returns a figure for the first chart on the Compare tab. It shows per day comparison between
    average sum of capacity by applied conditions.

    :param df: Pandas DataFrame, input data
    :param port1: str, a port to compare
    :param port2: str, a port to compare
    :param vessel_type: str, vessel type of interest
    :return: Plotly figure
    """
    data = helpers.filter_by_vessel_and_port(
        df=df, port1=port1, port2=port2, vessel_type=vessel_type
    )
    data_port_1 = (
        data[data["port"] == port1]
        .groupby(by="date")
        .sum()
        .reset_index()[["date", "sum_dwt"]]
    )
    data_port_2 = (
        data[data["port"] == port2]
        .groupby(by="date")
        .sum()
        .reset_index()[["date", "sum_dwt"]]
    )

    fig_data = [
        go.Bar(
            x=data_port_1["date"].tolist(),
            y=data_port_1["sum_dwt"].tolist(),
            name=port1,
            marker_color=styles.COLOR_APPSILON_1,
        ),
        go.Bar(
            x=data_port_2["date"].tolist(),
            y=data_port_2["sum_dwt"].tolist(),
            name=port2,
            marker_color=styles.COLOR_APPSILON_8,
        ),
    ]
    return go.Figure(
        data=fig_data,
        layout=styles.generate_plot_layout(
            x_title=strings.CHART_COMPARE_CAPACITY_X,
            y_title=strings.CHART_COMPARE_CAPACITY_Y,
            bar_mode="group",
        ),
    )


def get_compare_tab_insights(
    df: pd.DataFrame, port1: str, port2: str, vessel_type: str
) -> dict:
    """
    Returns data for comparison of increases/decreases in both ports and overall trend.

    :param df: Pandas DataFrame, input data
    :param port1: str, a port to compare
    :param port2: str, a port to compare
    :param vessel_type: str, vessel type of interest
    :return: insights dictionary
    """
    data = helpers.filter_by_vessel_and_port(
        df=df, port1=port1, port2=port2, vessel_type=vessel_type
    )
    data_port_1 = (
        data[data["port"] == port1]
        .groupby(by="date")
        .sum()
        .reset_index()[["date", "n"]]
    )
    data_port_2 = (
        data[data["port"] == port2]
        .groupby(by="date")
        .sum()
        .reset_index()[["date", "n"]]
    )

    def stats_checkup(x: int):
        if x == 0:
            return "-"
        return x

    def trend_direction(x: list) -> str:
        if len(x) == 0:
            return "-"
        detrend = signal.detrend(x)
        diff = x - detrend
        if diff[-1] > diff[0]:
            return "Positive"
        return "Negative"

    n_increase_port1 = stats_checkup((data_port_1["n"].diff() > 0).sum())
    n_decrease_port1 = stats_checkup((data_port_1["n"].diff() < 0).sum())
    n_increase_port2 = stats_checkup((data_port_2["n"].diff() > 0).sum())
    n_decrease_port2 = stats_checkup((data_port_2["n"].diff() < 0).sum())

    return {
        "Increase_P1": n_increase_port1,
        "Decrease_P1": n_decrease_port1,
        "Trend_P1": trend_direction(data_port_1["n"].tolist()),
        "Increase_P2": n_increase_port2,
        "Decrease_P2": n_decrease_port2,
        "Trend_P2": trend_direction(data_port_2["n"].tolist()),
    }
