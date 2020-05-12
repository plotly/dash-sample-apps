# ========== (c) JP Hwang 2020-04-10  ==========

import logging

# ===== START LOGGER =====
logger = logging.getLogger(__name__)
root_logger = logging.getLogger()
root_logger.setLevel(logging.INFO)
sh = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
sh.setFormatter(formatter)
root_logger.addHandler(sh)


import pandas as pd
import numpy as np

network_df = pd.read_csv("outputs/network_df.csv", index_col=0)
network_df["citations"] = network_df["citations"].fillna("")
network_df["cited_by"] = network_df["cited_by"].fillna("")

network_df_sm = network_df[
    (network_df["cited_by"] != "") | (network_df["citations"] != "")
]
old_indices = [str(i) for i in network_df_sm.index]
for i, row in network_df_sm.iterrows():
    rowloc = row.name
    network_df_sm.loc[rowloc, "cited_by"] = ",".join(
        [
            str(old_indices.index(j))
            for j in network_df_sm.loc[rowloc, "cited_by"].split(",")
            if len(j) > 0
        ]
    )
    network_df_sm.loc[rowloc, "citations"] = ",".join(
        [
            str(old_indices.index(j))
            for j in network_df_sm.loc[rowloc, "citations"].split(",")
            if len(j) > 0
        ]
    )
network_df_sm.reset_index(drop=True, inplace=True)
network_df_sm.to_csv("outputs/network_df_sm.csv")
