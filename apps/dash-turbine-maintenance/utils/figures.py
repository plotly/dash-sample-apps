from datetime import datetime, date
import plotly.graph_objs as go
from utils.helper_functions import df, df_button, x_test, fig_update_layout
import pickle
import pandas as pd


def update_graph(selected_column, start_date, end_date, n_get_new_info, n_pred):
    if n_pred is None:  # here is my work before prediction button is activated.
        value_rul = 0.0
        information_update = (
            "This field is used to display information about a feature displayed "
            "on the graph and estimated RUL. In order to estimate the RUL, use "
            "the button 'Get New Data' and then, 'Predict'. The estimated RUL will be "
            "printed."
        )
        if n_get_new_info is None:
            if selected_column in list(df):
                if start_date and end_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    end_date_object = datetime.strptime(end_date, "%Y-%m-%d")
                    mask = (df.index > start_date_object) & (
                        df.index <= end_date_object
                    )
                    df_within_dates = df.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                elif start_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    mask = df.index > start_date_object
                    df_within_dates = df.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                else:
                    fig = go.Figure(
                        data=[go.Scatter(x=df.index, y=df[selected_column])]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
            else:
                fig = go.Figure()
                fig = fig_update_layout(fig)
                return fig, value_rul, information_update
        else:
            if selected_column in list(df_button):
                if start_date and end_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    end_date_object = datetime.strptime(end_date, "%Y-%m-%d")
                    mask = (df_button.index > start_date_object) & (
                        df_button.index <= end_date_object
                    )
                    df_within_dates = df_button.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    _information_update = (
                        "New information is received for the last week and covers periods from "
                        + str(df_button.index[0])
                        + " to "
                        + str(df_button.index[-1])
                        + ". To predict"
                        " RUL, use 'Predict' button. To view data for the aforementioned period, choose"
                        " appropriate dates on the calendar."
                    )
                    return fig, value_rul, _information_update
                elif start_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    mask = df_button.index > start_date_object
                    df_within_dates = df_button.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    _information_update = (
                        "New information is received for the last week and covers periods from "
                        + str(df_button.index[0])
                        + " to "
                        + str(df_button.index[-1])
                        + ". To predict"
                        " RUL, use 'Predict' button. To view data for the aforementioned period, choose"
                        " appropriate dates on the calendar."
                    )
                    return fig, value_rul, _information_update
                else:
                    fig = go.Figure(
                        data=[
                            go.Scatter(x=df_button.index, y=df_button[selected_column])
                        ]
                    )
                    fig = fig_update_layout(fig)
                    _information_update = (
                        "New information is received for the last week and covers periods from "
                        + str(df_button.index[0])
                        + " to "
                        + str(df_button.index[-1])
                        + ". To predict"
                        " RUL, use 'Predict' button. To view data for the aforementioned period, choose"
                        " appropriate dates on the calendar."
                    )
                    return fig, value_rul, _information_update
            else:
                fig = go.Figure()
                fig = fig_update_layout(fig)
                _information_update = (
                    "New information is received for the last week and covers periods from "
                    + str(df_button.index[0])
                    + " to "
                    + str(df_button.index[-1])
                    + ". To predict"
                    " RUL, use 'Predict' button. To view data for the aforementioned period, choose"
                    " appropriate dates on the calendar."
                )
                return fig, value_rul, _information_update
    else:
        if n_get_new_info is None:
            value_rul = 0.0
            information_update = "To predict RUL, please use 'Get New Data' button."
            if selected_column in list(df):
                if start_date and end_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    end_date_object = datetime.strptime(end_date, "%Y-%m-%d")
                    mask = (df.index > start_date_object) & (
                        df.index <= end_date_object
                    )
                    df_within_dates = df.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                elif start_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    mask = df.index > start_date_object
                    df_within_dates = df.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                else:
                    fig = go.Figure(
                        data=[go.Scatter(x=df.index, y=df[selected_column])]
                    )
                    return fig, value_rul, information_update
            else:
                fig = go.Figure()
                fig = fig_update_layout(fig)
                return fig, value_rul, information_update
        else:
            model = pickle.load(open("assets/xgb_reg.pkl", "rb"))
            y_pred = model.predict(x_test)
            df_out = pd.DataFrame()
            df_out["pred"] = y_pred
            value_rul = round(max(df_out["pred"]))
            information_update = (
                "RUL is estimated based on the readings from the last week: "
                "from " + str(x_test.index[0]) + " to " + str(x_test.index[-1])
            )
            if selected_column in list(df_button):
                if start_date and end_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    end_date_object = datetime.strptime(end_date, "%Y-%m-%d")
                    mask = (df_button.index > start_date_object) & (
                        df_button.index <= end_date_object
                    )
                    df_within_dates = df_button.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                elif start_date:
                    start_date_object = datetime.strptime(start_date, "%Y-%m-%d")
                    mask = df_button.index > start_date_object
                    df_within_dates = df_button.loc[mask]
                    fig = go.Figure(
                        data=[
                            go.Scatter(
                                x=df_within_dates.index,
                                y=df_within_dates[selected_column],
                            )
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
                else:
                    fig = go.Figure(
                        data=[
                            go.Scatter(x=df_button.index, y=df_button[selected_column])
                        ]
                    )
                    fig = fig_update_layout(fig)
                    return fig, value_rul, information_update
            else:
                fig = go.Figure()
                fig = fig_update_layout(fig)
                return fig, value_rul, information_update


def display_click_data(clickData):
    if clickData:
        data_time = clickData["points"][0]["x"]
        value_active_power = df["WEC: ava. Power"].loc[df.index == data_time].values[0]
        value_active_power_wind = (
            df["WEC: ava. available P from wind"].loc[df.index == data_time].values[0]
        )
        value_reactive_power = (
            df["WEC: ava. reactive Power"].loc[df.index == data_time].values[0]
        )
        value_wind_speed = (
            df["WEC: ava. windspeed"].loc[df.index == data_time].values[0]
        )
        value_blade_angle = (
            df["WEC: ava. blade angle A"].loc[df.index == data_time].values[0]
        )
        return (
            value_active_power,
            value_active_power_wind,
            value_wind_speed,
            value_reactive_power,
            value_blade_angle,
        )
    else:
        value_active_power = 0
        value_active_power_wind = 0
        value_reactive_power = 0
        value_wind_speed = 0
        value_blade_angle = 0
        return (
            value_active_power,
            value_active_power_wind,
            value_wind_speed,
            value_reactive_power,
            value_blade_angle,
        )
