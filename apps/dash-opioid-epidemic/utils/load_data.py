import os
import pandas as pd


# Load data
df_lat_lon = pd.read_csv(
    os.path.join(os.path.dirname(__file__), "../data/lat_lon_counties.csv")
)

df_lat_lon["FIPS "] = df_lat_lon["FIPS "].apply(lambda x: str(x).zfill(5))

df_full_data = pd.read_csv(
    os.path.join(
        os.path.dirname(__file__), "../data/age_adjusted_death_rate_no_quotes.csv"
    )
)

df_full_data["County Code"] = df_full_data["County Code"].apply(
    lambda x: str(x).zfill(5)
)

df_full_data["County"] = (
    df_full_data["Unnamed: 0"] + ", " + df_full_data.County.map(str)
)
