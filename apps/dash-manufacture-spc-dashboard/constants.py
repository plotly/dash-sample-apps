import pandas as pd
from utils.helper_functions import init_df

## Load data used across the app
df = pd.read_csv("data/spc_data.csv")
state_dict = init_df(df)

## Define constants used across the app
params = list(df)
max_length = len(df)

suffix_row = "_row"
suffix_button_id = "_button"
suffix_sparkline_graph = "_sparkline_graph"
suffix_count = "_count"
suffix_ooc_n = "_OOC_number"
suffix_ooc_g = "_OOC_graph"
suffix_indicator = "_indicator"