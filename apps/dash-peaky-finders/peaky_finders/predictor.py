from typing import List
import datetime as dt
from datetime import timedelta
import requests
import os
import pickle
from typing import Dict, Tuple

import geopandas as gpd
import pandas as pd
import pytz, datetime
from shapely import wkt
from timezonefinderL import TimezoneFinder

from peaky_finders.data_acquisition.train_model import (
    LoadCollector,
    GEO_COORDS,
    CATEGORICAL_FEATURES,
    MONTH_TO_SEASON,
)
from peaky_finders.training_pipeline import MODEL_OUTPUT_DIR, MODEL_INPUT_DIR
from peaky_finders.data_acquisition.train_model import GEO_COORDS


ISO_MAP_IDS = {
    56669: "MISO",
    14725: "PJM",
    2775: "CAISO",
    13434: "ISONE",
    13501: "NYISO",
}

ISO_LIST = ["NYISO", "ISONE", "PJM", "MISO", "CAISO"]

PEAK_DATA_PATH = os.path.join(os.path.dirname(__file__), "historical_peaks")


tz_finder = TimezoneFinder()


def get_iso_map():
    iso_df = pd.read_csv("iso_map_final.csv")
    iso_df["geometry"] = iso_df["geometry"].apply(wkt.loads)
    iso_gdf = gpd.GeoDataFrame(iso_df, crs="EPSG:4326", geometry="geometry")
    return iso_gdf


class Predictor:
    def __init__(self, iso_name: str, start: str, end: str) -> None:
        self.start = start
        self.end = end
        self.iso_name = iso_name
        self.load_collector: LoadCollector = None

    def get_load(self, start: str, end: str):
        self.load_collector = LoadCollector(self.iso_name, start, end)

    def featurize(self):
        self.load_collector.engineer_features()

    def add_future(self, load: pd.Series) -> pd.Series:
        future = pd.date_range(
            start=load.index[-1], end=(load.index[-1] + timedelta(days=1)), freq="H"
        ).to_frame(name="load_MW")
        tz_finder = TimezoneFinder()
        lon = float(GEO_COORDS[self.iso_name]["lon"])
        lat = float(GEO_COORDS[self.iso_name]["lat"])
        tz_name = tz_finder.timezone_at(lng=lon, lat=lat)
        future["load_MW"] = None
        future.index = future.index.tz_convert(tz_name)
        return future

    def prepare_predictions(self):
        self.get_load(self.start, self.end)
        load = self.load_collector.load
        self.load_collector.engineer_features()
        model_input = self.load_collector.load.copy()
        for feature in CATEGORICAL_FEATURES:
            dummies = pd.get_dummies(
                model_input[feature], prefix=feature, drop_first=True
            )
            model_input = model_input.drop(feature, axis=1)
            model_input = pd.concat([model_input, dummies], axis=1)
        return model_input, load

    def predict_load(self, model_input: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
        model_path = os.path.join(
            MODEL_OUTPUT_DIR, (f"xg_boost_{self.iso_name}_load_model.pkl")
        )
        xgb = pickle.load(open(model_path, "rb"))
        if "holiday_True" not in model_input.columns:
            model_input["holiday_True"] = 0
        X = model_input.drop("load_MW", axis=1).astype(float).dropna()
        weekday_cols = [f"weekday_{i + 1}" for i in range(0, 6)]
        if len(set(weekday_cols) - set(X.columns)) > 0:
            for col in list(set(weekday_cols) - set(X.columns)):
                X[col] = 0
        predictions = xgb.predict(X[xgb.get_booster().feature_names])
        X["predicted_load"] = predictions
        return X["predicted_load"]


def predict_load(self,):
    for iso in ISO_LIST:
        model_input_path = os.path.join(MODEL_INPUT_DIR, MODEL_INPUT_DATA[iso])
        model_path = os.path.join(
            MODEL_OUTPUT_DIR, (f"xg_boost_{self.iso_name}_load_model.pkl")
        )


def predict_all(iso_list: list, start: str, end: str) -> Tuple[Dict[str, pd.DataFrame]]:
    historical_vs_predicted = {}
    for iso in iso_list:
        predictor = Predictor(iso, start, end)
        model_input, historical_load = predictor.prepare_predictions()
        predictions = predictor.predict_load(model_input)
        comparison_df = pd.concat([model_input, predictions], axis=1)[
            ["load_MW", "predicted_load"]
        ]
        historical_vs_predicted[iso] = comparison_df
    return historical_vs_predicted


def get_peak_data(iso_list: list) -> Tuple[Dict[str, pd.DataFrame]]:
    peak_data = {}
    for iso in iso_list:
        path = "https://raw.githubusercontent.com/kbaranko/peaky-finders/master/peaky_finders/historical_peaks"
        iso_data = pd.read_csv(
            f"{path}/{iso}_historical_peaks.csv", parse_dates=["timestamp"]
        )
        iso_data["timestamp"] = iso_data["timestamp"].apply(
            lambda x: x.astimezone(pytz.utc)
        )
        tz_name = tz_finder.timezone_at(
            lng=float(GEO_COORDS[iso]["lon"]), lat=float(GEO_COORDS[iso]["lat"])
        )
        iso_data.index = pd.DatetimeIndex(iso_data["timestamp"])
        iso_data.index = iso_data.index.tz_convert(tz_name)
        peak_data[iso] = iso_data
    return peak_data


def get_temperature_forecast(iso: str) -> dict:
    lon = GEO_COORDS[iso]["lon"]
    lat = GEO_COORDS[iso]["lat"]
    API_KEY = os.environ["DARKSKY_KEY"]
    url = f"https://api.darksky.net/forecast/{API_KEY}/{lat},{lon}"
    response = requests.get(url)
    if response.status_code == 200:
        print(response.status_code)
    else:
        raise ValueError(
            f"Error getting data from DarkSky API: "
            f"Response Code {response.status_code}"
        )
    info = response.json()
    hourly_data = info["hourly"]["data"]
    hourly_temp = {}
    for info in hourly_data:
        timestamp = datetime.datetime.fromtimestamp(info["time"])
        tz = tz_finder.timezone_at(lng=float(lon), lat=float(lat))
        timestamp = timestamp.astimezone(pytz.timezone(tz))
        hourly_temp[timestamp] = info["temperature"]
    return hourly_temp


def create_load_duration(peak_data: Dict[str, pd.DataFrame]) -> Dict[str, pd.Series]:
    load_duration_curves = {}
    for iso in ISO_LIST:
        load_duration_curves[iso] = pd.Series(
            peak_data[iso]["load_MW"].values
        ).sort_values(ascending=False)
    return load_duration_curves


def get_forecasts(iso_list: List[str]):
    predictions = {}
    historical_load = {}
    temperature = {}
    for iso in iso_list:
        path = f"https://raw.githubusercontent.com/kbaranko/peaky-finders/master/peaky_finders/forecasts/{iso}_forecasts.csv"
        iso_data = pd.read_csv(path, parse_dates=["timestamp"])
        iso_data["timestamp"] = iso_data["timestamp"].apply(
            lambda x: x.astimezone(pytz.utc)
        )
        tz_name = tz_finder.timezone_at(
            lng=float(GEO_COORDS[iso]["lon"]), lat=float(GEO_COORDS[iso]["lat"])
        )
        iso_data.index = pd.DatetimeIndex(iso_data["timestamp"])
        iso_data.index = iso_data.index.tz_convert(tz_name)
        historical_load[iso] = iso_data["load_MW"]
        predictions[iso] = iso_data["predicted_load"]
        temperature[iso] = iso_data["temperature"]
    return predictions, historical_load, temperature
