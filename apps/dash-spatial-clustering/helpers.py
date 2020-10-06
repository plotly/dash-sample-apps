import json
import os
import shutil
import time

import pandas as pd
import numpy as np
import pathlib
import geopandas as gpd
import pysal as ps

from sklearn import cluster
from sklearn.preprocessing import scale

# Data reading & Processing
app_path = pathlib.Path(__file__).parent.resolve()
data_path = pathlib.Path(__file__).parent.joinpath("data")
geo_json_path = data_path.joinpath("Zipcodes.geojson")
austin_listings = pd.read_csv(
    "https://raw.githubusercontent.com/plotly/datasets/master/dash-sample-apps/dash-spatial-clustering/data/listings.csv",
    low_memory=False,
)

# Refractor zipcode outlier, modify in place
zip_outlier = austin_listings[austin_listings["zipcode"] == "TX 78702"].index
austin_listings.loc[zip_outlier, "zipcode"] = "78702"
austin_listings = austin_listings.dropna(axis=0, subset=["zipcode"])


review_columns = [c for c in austin_listings.columns if "review_" in c]

# Geojson loading
with open(geo_json_path) as response:
    zc_link = json.load(response)
    # Add id for choropleth layer
    for feature in zc_link["features"]:
        feature["id"] = feature["properties"]["zipcode"]

listing_zipcode = austin_listings["zipcode"].unique()


def apply_clustering():
    """
     # Apply KMeans clustering to group zipcodes into categories based on type of houses listed(i.e. property type)

    :return: Dataframe.
    db: scaled proportions of house types by zipcode, use for plotting Choropleth map layer.
    barh_df : scaled proportion of house type grouped by cluster, use for prop type chart and review chart.
    """
    variables = ["bedrooms", "bathrooms", "beds"]
    aves = austin_listings.groupby("zipcode")[variables].mean()
    review_aves = austin_listings.groupby("zipcode")[review_columns].mean()

    types = pd.get_dummies(austin_listings["property_type"])
    prop_types = types.join(austin_listings["zipcode"]).groupby("zipcode").sum()
    prop_types_pct = (prop_types * 100.0).div(prop_types.sum(axis=1), axis=0)

    aves_props = aves.join(prop_types_pct)

    # Standardize a dataset along any axis, Center to the mean and component wise scale to unit variance.
    db = pd.DataFrame(
        scale(aves_props), index=aves_props.index, columns=aves_props.columns
    ).rename(lambda x: str(x))

    # Apply clustering on scaled df
    km5 = cluster.KMeans(n_clusters=5)
    km5cls = km5.fit(db.reset_index().values)
    # print(len(km5cls.labels_))

    db["cl"] = km5cls.labels_

    # sort by labels since every time cluster is running, label 0-4 is randomly assigned
    db["count"] = db.groupby("cl")["cl"].transform("count")

    db.sort_values("count", inplace=True, ascending=True)

    barh_df = prop_types_pct.assign(cl=km5cls.labels_).groupby("cl").mean()

    # Join avg review columns for updating review plot
    db = db.join(review_aves)
    grouped = db.groupby("cl")[review_columns].mean()
    barh_df = barh_df.join(grouped)

    return db.reset_index(), barh_df


def rating_clustering(threshold):
    start = time.time()
    # Explore boundaries/ areas where customers are have similar ratings. Different from
    # predefined number of output regions, it takes target variable(num of reviews, and
    # apply a minimum threshold (5% per region) on it.

    # Bring review columns at zipcode level
    rt_av = austin_listings.groupby("zipcode")[review_columns].mean().dropna()

    # Regionalization requires building of spatial weights
    zc = gpd.read_file(geo_json_path)
    zrt = zc[["geometry", "zipcode"]].join(rt_av, on="zipcode").dropna()
    zrt.to_file("tmp")
    w = ps.queen_from_shapefile("tmp/tmp.shp", idVariable="zipcode")

    # Remove temp tmp/* we created for spatial weights
    if os.path.isdir(os.path.join(app_path, "tmp")):
        print("removing tmp folder")
        shutil.rmtree(os.path.join(app_path, "tmp"))

    # Impose that every resulting region has at least 5% of the total number of reviews
    n_review = (
        austin_listings.groupby("zipcode")
        .sum()["number_of_reviews"]
        .rename(lambda x: str(int(x)))
        .reindex(zrt["zipcode"])
    )
    thr = np.round(int(threshold) / 100 * n_review.sum())

    # Set the seed for reproducibility
    np.random.seed(1234)

    z = zrt.drop(["geometry", "zipcode"], axis=1).values
    # Create max-p algorithm, note that this API is upgraded in pysal>1.11.1
    maxp = ps.region.Maxp(w, z, thr, n_review.values[:, None], initial=100)
    maxp.cinference(nperm=99)
    # p value compared with randomly assigned region
    p_value = maxp.cpvalue
    print("p_value:", p_value)
    lbls = pd.Series(maxp.area2region).reindex(zrt["zipcode"])
    regionalization_df = (
        pd.DataFrame(lbls).reset_index().rename(columns={"zipcode": "zipcode", 0: "cl"})
    )

    end = time.time()
    # The larger threshold, the longer time it takes for computing
    print(
        "Computing threshold {}%".format(threshold),
        "time cost for clustering: ",
        end - start,
    )

    types = pd.get_dummies(austin_listings["property_type"])
    prop_types = types.join(austin_listings["zipcode"]).groupby("zipcode").sum()
    merged = pd.merge(
        prop_types.reset_index(), regionalization_df, on="zipcode", how="inner"
    )
    d_merged = merged.drop(["zipcode", "cl"], axis=1)
    prop_types_pct = (d_merged * 100.0).div(d_merged.sum(axis=1), axis=0)

    pct_d = (
        prop_types_pct.assign(cl=merged["cl"], zipcode=merged["zipcode"])
        .groupby("cl")
        .mean()
    )

    zrt = zrt[review_columns].groupby(lbls.values).mean()
    joined_prop = pct_d.join(zrt)
    return regionalization_df, p_value, joined_prop


# #
# rating = rating_clustering(5)
