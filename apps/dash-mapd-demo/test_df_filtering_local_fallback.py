import pandas as pd
import pathlib

DATA_PATH = pathlib.Path(__file__).parent.joinpath("dataset")
fallback_df = pd.read_csv(DATA_PATH.joinpath("init.csv"), index_col=0).dropna()

fallback_df['arr_timestamp'] = pd.to_datetime(fallback_df['arr_timestamp'], errors='coerce')
fallback_df['dep_timestamp'] = pd.to_datetime(fallback_df['dep_timestamp'], errors='coerce')

# choro-select
# choro_query = f"SELECT AVG(depdelay) AS avg_delay, {state_col} AS state FROM " \
 #   f"{table} WHERE dep_timestamp BETWEEN '{start_f}' AND '{end_f}' GROUP BY {state_col}"
# SELECT AVG(depdelay) AS avg_delay, origin_state AS state FROM flights_2008_7M WHERE dep_timestamp BETWEEN
# '2008-01-02 00:00:00' AND '2008-01-02 00:00:00' GROUP BY origin_state
start_f = '2008-01-02 00:00:00'
end_f = '2008-01-09 00:00:00'

choro_df = fallback_df[(fallback_df['dep_timestamp']> start_f) & (fallback_df['dep_timestamp']<end_f)]

choro_df.groupby('origin_state')
print(choro_df)