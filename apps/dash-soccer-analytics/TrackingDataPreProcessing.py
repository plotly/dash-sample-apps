import pandas as pd
import numpy as np

# raw_df = pd.read_csv(cv_file, names=header_list, error_bad_lines=False, dtype=str)
raw_df1 = pd.read_csv("data/source/Sample_Game_2_RawTrackingData_Away_Team.csv")
# sample every 5th row
raw_df1 = raw_df1.iloc[::7, :]

raw_df2 = pd.read_csv("data/source/Sample_Game_2_RawTrackingData_Home_Team.csv")
raw_df2 = raw_df2.iloc[::7, :]

column = 1
df = pd.DataFrame(columns=["half", "frame", "time", "x", "y"])
# x range needs adjusted depending on how many columns there are. Should really calculate this not eyeball items
# and do it manually.
for x in range(0, 13):
    column = column + 2
    df_temp = raw_df1.iloc[:, [0, 1, 2, column, column + 1]].copy()
    df_temp.columns = ["half", "frame", "time", "x", "y"]
    df_temp["jersey_number"] = raw_df1.columns[column]
    df = pd.concat([df, df_temp]).reset_index(drop=True)
df["team"] = "Away"
df.loc[df["jersey_number"] == "0", "team"] = "Ball"
df.loc[df["x"].isna(), "x"] = None
df.loc[df["y"].isna(), "y"] = None
df = df[df["x"].notna()]
df.drop(df.loc[df["half"] == "Period"].index, inplace=True)

column = 1
df2 = pd.DataFrame(columns=["half", "frame", "time", "x", "y"])
for x in range(0, 12):
    column = column + 2
    df_temp2 = raw_df2.iloc[:, [0, 1, 2, column, column + 1]].copy()
    df_temp2.columns = ["half", "frame", "time", "x", "y"]
    df_temp2["jersey_number"] = raw_df2.columns[column]
    df2 = pd.concat([df2, df_temp2]).reset_index(drop=True)
df2["team"] = "Home"
df2.loc[df2["x"].isna(), "x"] = None
df2.loc[df2["y"].isna(), "y"] = None
df2 = df2[df2["x"].notna()]
df2.drop(df2.loc[df2["half"] == "Period"].index, inplace=True)

df = df.iloc[1:]
df["frame"] = df["frame"].apply(pd.to_numeric, errors="coerce")
df = df.sort_values(by=["frame"])

df2 = df2.iloc[1:]
df2["frame"] = df2["frame"].apply(pd.to_numeric, errors="coerce")
df2 = df2.sort_values(by=["frame"])

df_export = pd.concat([df, df2]).reset_index(drop=True)
df_export = df_export.sort_values(by=["frame"])
df_export["time"] = df_export["time"].apply(pd.to_numeric, errors="coerce")
df_export["time"] = df_export["time"].div(60).round(4)
export_file_name = input(
    "Please enter a name for the file to be exported (ending with .csv): "
)
export_file_name = "data/" + export_file_name
df_export.to_csv(export_file_name, index=False)
