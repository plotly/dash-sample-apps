# ========== (c) JP Hwang 2020-04-02  ==========

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
import os
import json

desired_width = 320
pd.set_option("display.max_columns", 20)
pd.set_option("display.width", desired_width)

if "cord-nlp-cytoscape" in os.listdir("."):
    os.chdir("cord-nlp-cytoscape")

lda_df = pd.read_csv("outputs/lda_df.csv", index_col=0)
meta_df = pd.read_csv("metadata.csv")
meta_df["sha"] = meta_df["sha"].fillna("")
datadir = "comm_use_subset"

# Add citation data to the DataFrame
cite_links = list()
journals = list()
pub_dates = list()
authors = list()
counter = 0
cited_by = {str(i): [] for i in range(len(lda_df))}
logger.info("Starting analysis...")
for i, row in lda_df.iterrows():
    datafile = row["filename"]
    data_sha = row["filename"][:-5]

    # Add citation data to the DataFrame
    conn_nodes = list()
    with open(os.path.join(datadir, datafile), "r") as f:
        data_dict = json.load(f)
    for k, v in data_dict["bib_entries"].items():
        ref_title = v["title"].lower()
        if ref_title in lda_df["title"].values:
            tgt_id = str(lda_df[lda_df["title"] == ref_title].index[0])
            if str(i) != tgt_id:
                conn_nodes.append(tgt_id)
                cited_by[tgt_id].append(str(i))
                counter += 1
    cite_links.append(",".join(conn_nodes))

    # add journal & publication (meta) data to the lda dataframe
    meta_row = meta_df[meta_df["sha"].str.contains(data_sha)]
    if len(meta_row) > 0:
        journals.append(meta_row["journal"].values[0])
        pub_dates.append(meta_row["publish_time"].values[0])
        authors.append(meta_row["authors"].values[0])
    else:
        logger.error(f"Warning - {datafile} not found in Metadata")
        journals.append("No data")
        pub_dates.append("No data")
        authors.append("No data")

    if (i + 1) % 500 == 0:
        logger.info(
            f"Processed {i+1} files out of {len(lda_df)}, found {counter} references"
        )

logger.info(f"Setting new columns for citations, journal & publication date")
network_df = lda_df.assign(citations=cite_links)
network_df["citations"] = network_df["citations"].fillna("")
network_df = network_df.assign(journal=journals)
network_df = network_df.assign(pub_date=pub_dates)
network_df = network_df.assign(authors=authors)

cited_by_list = [",".join(cited_by[str(i)]) for i in range(len(network_df))]
network_df = network_df.assign(cited_by=cited_by_list)
network_df = network_df.assign(
    n_cites=[len(i.split(",")) if len(i) > 0 else 0 for i in cited_by_list]
)

for col in ["title", "journal", "pub_date", "authors"]:
    network_df[col] = network_df[col].fillna("No data")

# TODO - filter out data points with missing data (e.g. no title)
# network_df = network_df[-network_df.title.isna()]

network_df_file = "outputs/network_df.csv"
logger.info(f"Saving file {network_df_file}")
network_df.to_csv(network_df_file)
