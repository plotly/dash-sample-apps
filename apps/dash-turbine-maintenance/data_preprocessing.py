import pandas as pd
from sklearn.model_selection import train_test_split

pd.options.mode.chained_assignment = None


def data_preprocessing():
    df = pd.read_csv("data/SCADA_data.csv")
    status_data_wec = pd.read_csv("data/status_data_wec.csv")

    df["Inverter avg. temp"] = df[
        [
            "CS101 : Sys 1 inverter 1 cabinet temp.",
            "CS101 : Sys 1 inverter 2 cabinet temp.",
            "CS101 : Sys 1 inverter 3 cabinet temp.",
            "CS101 : Sys 1 inverter 4 cabinet temp.",
            "CS101 : Sys 1 inverter 5 cabinet temp.",
            "CS101 : Sys 1 inverter 6 cabinet temp.",
            "CS101 : Sys 1 inverter 7 cabinet temp.",
            "CS101 : Sys 2 inverter 1 cabinet temp.",
            "CS101 : Sys 2 inverter 2 cabinet temp.",
            "CS101 : Sys 2 inverter 3 cabinet temp.",
            "CS101 : Sys 2 inverter 4 cabinet temp.",
        ]
    ].mean(axis=1)
    df["Inverter std. temp"] = df[
        [
            "CS101 : Sys 1 inverter 1 cabinet temp.",
            "CS101 : Sys 1 inverter 2 cabinet temp.",
            "CS101 : Sys 1 inverter 3 cabinet temp.",
            "CS101 : Sys 1 inverter 4 cabinet temp.",
            "CS101 : Sys 1 inverter 5 cabinet temp.",
            "CS101 : Sys 1 inverter 6 cabinet temp.",
            "CS101 : Sys 1 inverter 7 cabinet temp.",
            "CS101 : Sys 2 inverter 1 cabinet temp.",
            "CS101 : Sys 2 inverter 2 cabinet temp.",
            "CS101 : Sys 2 inverter 3 cabinet temp.",
            "CS101 : Sys 2 inverter 4 cabinet temp.",
        ]
    ].std(axis=1)
    df["Time"] = pd.to_datetime(df["Time"], dayfirst=True, errors="coerce")
    df.sort_values(by="Time", axis=0, inplace=True)
    df.reset_index(drop=True, inplace=True)

    af_corr_time_wec_s = status_data_wec.loc[
        (status_data_wec["Main Status"] == 62)
        | (status_data_wec["Main Status"] == 80)
        | (status_data_wec["Main Status"] == 228)
        | (status_data_wec["Main Status"] == 60)
        | (status_data_wec["Main Status"] == 9),
        "Time",
    ]

    af_corr_time_wec_s = pd.to_datetime(af_corr_time_wec_s)
    af_corr_time_wes = af_corr_time_wec_s.round("10min")

    df.rename(columns={"Error": "Fault"}, inplace=True)
    df["Fault"] = [0] * len(df)
    for i, d in enumerate(df["Time"]):
        if d in af_corr_time_wes.values:
            df["Fault"][i] = 1

    nf_times = []
    rul = []
    for i, d in enumerate(df["Fault"]):
        nf_times.append(df["Time"][i])
        if d == 1:
            for j in nf_times:
                rul.append(df["Time"][i] - j)
            nf_times = []

    df_trimmed = df.head(len(rul))

    rul_days = [x.days for x in rul]
    df_trimmed["RUL"] = rul_days

    cols_to_drop = [
        "Fault",
        "CS101 : Sys 1 inverter 1 cabinet temp.",
        "CS101 : Sys 1 inverter 2 cabinet temp.",
        "CS101 : Sys 1 inverter 3 cabinet temp.",
        "CS101 : Sys 1 inverter 4 cabinet temp.",
        "CS101 : Sys 1 inverter 5 cabinet temp.",
        "CS101 : Sys 1 inverter 6 cabinet temp.",
        "CS101 : Sys 1 inverter 7 cabinet temp.",
        "CS101 : Sys 2 inverter 1 cabinet temp.",
        "CS101 : Sys 2 inverter 2 cabinet temp.",
        "CS101 : Sys 2 inverter 3 cabinet temp.",
        "CS101 : Sys 2 inverter 4 cabinet temp.",
    ]

    for i in cols_to_drop:
        if i in list(df):
            df_trimmed.drop(i, axis=1, inplace=True)

    df_trimmed = df_trimmed.head(39298)
    df_trimmed.set_index("Time", inplace=True)

    df = df_trimmed.copy()
    features = df.columns.tolist()
    timesteps = 5

    df_list = [
        df[features].shift(shift_val)
        if (shift_val == 0)
        else df[features].shift(-shift_val).add_suffix(f"_{shift_val}")
        for shift_val in range(0, timesteps)
    ]
    df_concat = pd.concat(df_list, axis=1, sort=False)
    df_concat = df_concat.iloc[:-timesteps, :]

    target = "RUL"
    x_train, x_test, y_train, y_test = train_test_split(
        df_concat,
        df[target].iloc[:-timesteps],
        test_size=0.02642894598,
        random_state=10,
        shuffle=False,
    )

    df_test = x_test.iloc[:, : df.shape[1]]

    return df, df_test, x_test, y_test
