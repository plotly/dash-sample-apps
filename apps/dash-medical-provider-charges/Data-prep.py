import pandas as pd
import numpy as np
import pathlib

from uszipcode import SearchEngine

DATA_PATH = pathlib.Path(__file__, "/data").resolve()
all_data = DATA_PATH.joinpath(
    "CMS_inpatient_2016_data/Inpatient_summary_2016_all/Medicare_Provider_Charge_Inpatient_DRGALL_FY2016.csv.gz"
)

# Read all
df = pd.read_csv(str(all_data)[1:], index_col=0, low_memory=False)

state = [
    "AL",
    "AZ",
    "AR",
    "CA",
    "CO",
    "CT",
    "DC",
    "FL",
    "GA",
    "IL",
    "IN",
    "IA",
    "KY",
    "MD",
    "MA",
    "MI",
    "MN",
    "MO",
    "NE",
    "NJ",
    "NY",
    "NC",
    "OH",
    "OK",
    "OR",
    "PA",
    "SC",
    "TN",
    "TX",
    "UT",
    "VA",
    "WA",
    "WI",
    "AK",
    "DE",
    "HI",
    "ID",
    "KS",
    "LA",
    "ME",
    "MS",
    "MT",
    "NV",
    "NH",
    "NM",
    "ND",
    "RI",
    "SD",
    "WV",
    "VT",
    "WY",
]

# Generate lat, lon from zip code
search = SearchEngine(
    simple_zipcode=False
)  # set simple_zipcode=False to use rich info database


def generate_lat_lon(df):
    df["lat"] = df["lon"] = np.nan
    zip_code = df["Provider Zip Code"]
    for ind, item in enumerate(zip_code):
        zipcode = search.by_zipcode(str(item))
        zip = zipcode.to_dict()
        df.loc[ind, "lat"] = zip["lat"]
        df.loc[ind, "lon"] = zip["lng"]
    return df


# Replace charges dtype as numeric
metric_list = [
    "Average Covered Charges",
    "Average Total Payments",
    "Average Medicare Payments",
]

for state_name in state:
    state_data = df[df["Provider State"] == state_name].reset_index()
    state_lat_lon = generate_lat_lon(state_data).dropna()
    for metric in metric_list:
        state_lat_lon[metric] = state_lat_lon[metric].apply(
            lambda x: x.replace("$", "")
        )
        state_lat_lon[metric] = state_lat_lon[metric].apply(
            lambda x: float(x.split()[0].replace(",", ""))
        )
    state_lat_lon.to_csv("df_{}_lat_lon.csv".format(state_name), index=False)
