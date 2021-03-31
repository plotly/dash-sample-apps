import calendar
import datetime as dt
from datetime import timedelta
import holidays
import math
import os

from dateutil.relativedelta import relativedelta
import json
import numpy as np
import pandas as pd
import pickle
from pyiso import client_factory
from pyiso.eia_esod import EIAClient
import requests
from sklearn import preprocessing
from timezonefinderL import TimezoneFinder


GEO_COORDS = {
    "NYISO": {"lat": "40.7128", "lon": "-73.935242"},
    "ISONE": {"lat": "42.3601", "lon": "-71.0589"},
    "CAISO": {"lat": "34.0522", "lon": "-118.2437"},
    "PJM": {"lat": "39.9526", "lon": "-75.1652"},
    "MISO": {"lat": "44.9778", "lon": "-93.2650",},
}

MONTH_TO_SEASON = {
    1: "Winter",
    2: "Winter",
    3: "Spring",
    4: "Spring",
    5: "Spring",
    6: "Summer",
    7: "Summer",
    8: "Summer",
    9: "Fall",
    10: "Fall",
    11: "Fall",
    12: "Winter",
}


BASE_URL = "https://api.darksky.net/forecast"
EXCLUDE = "flags, minutely, daily, alerts"

LOAD_COLS = ["load_MW", "timestamp"]
EASTERN_TZ = "US/Eastern"

US_HOLIDAYS = holidays.UnitedStates()
CATEGORICAL_FEATURES = ["weekday", "hour_of_day", "holiday"]
NUMERICAL_FEATURES = ["temperature", "load (t-24)"]


class LoadCollector:
    def __init__(self, iso: str, start_date: str, end_date: str):
        self.start = start_date
        self.end = end_date
        self.iso_name = iso
        self.lat = GEO_COORDS[iso]["lat"]
        self.lon = GEO_COORDS[iso]["lon"]
        self.iso = self._set_iso(iso)
        self.holidays = holidays.UnitedStates()
        self.load = self.get_historical_load()
        self.model_input = None

    def get_historical_load(self) -> pd.DataFrame:
        if self.iso_name == "CAISO":
            load = self.get_caiso_load()
        elif (
            self.iso_name == "MISO"
            or self.iso_name == "PJM"
            or self.iso_name == "ERCOT"
        ):
            load = self.get_eia_load()
        else:
            load = pd.DataFrame(
                self.iso.get_load(
                    latest=False, yesterday=False, start_at=self.start, end_at=self.end
                )
            )[LOAD_COLS].set_index("timestamp")
        tz_finder = TimezoneFinder()
        tz_name = tz_finder.timezone_at(lng=float(self.lon), lat=float(self.lat))
        load.index = load.index.tz_convert(tz_name)
        return load.resample("H").mean()

    def get_historical_peak_load(self) -> pd.DataFrame:
        daily_peak = self.load.resample("D").max()
        holiday_bool = dict()
        for date, _ in daily_peak.iterrows():
            holiday_bool[date] = self._check_for_holiday(date)
        daily_peak["month"] = daily_peak.index.month_name()
        daily_peak["season"] = daily_peak.index.month.map(MONTH_TO_SEASON)
        daily_peak["weekday"] = daily_peak.index.day_name()
        daily_peak["holiday"] = daily_peak.index.map(holiday_bool)
        return daily_peak

    def get_eia_load(self):
        load = pd.DataFrame(
            self.iso.get_load(
                latest=True, yesterday=False, start_at=self.start, end_at=self.end
            )
        )
        load = load.iloc[::-1]
        return load[LOAD_COLS].set_index("timestamp")

    def get_caiso_load(self):
        if pd.Timestamp(self.start).month == pd.Timestamp(self.end).month:
            months = [pd.Timestamp(self.start)]
        else:
            months = pd.date_range(self.start, self.end, freq="MS").tolist()
        monthly_load = []
        for month in months:
            start, end = self.get_month_day_range(month)
            start = start.strftime("%Y-%m-%d")
            end = end.strftime("%Y-%m-%d")
            load = pd.DataFrame(
                self.iso.get_load(
                    latest=False, yesterday=False, start_at=start, end_at=end
                )
            )[LOAD_COLS].set_index("timestamp")
            monthly_load.append(load)
        return pd.concat(monthly_load)

    @staticmethod
    def get_month_day_range(date):
        """
        Returns the start and end date for the month of 'date'.
        """
        last_day = date + relativedelta(day=1, months=+1, days=-1)
        first_day = date + relativedelta(day=1)
        return first_day, last_day

    def engineer_features(self):
        temperatures = dict()
        holiday_bool = dict()
        for date, _ in self.load.iterrows():
            temperatures[date] = self._get_temperature(date)
            holiday_bool[date] = self._check_for_holiday(date)
        self.load["weekday"] = self.load.index.dayofweek
        self.load["hour_of_day"] = self.load.index.hour
        self.load["temperature"] = self.load.index.map(temperatures)
        self.load["holiday"] = self.load.index.map(holiday_bool)
        self.load["load (t-24)"] = self.load.load_MW.shift(24)

    def engineer_features_lite(self, weather_dict: dict):
        holiday_bool = dict()
        for date, _ in self.load.iterrows():
            holiday_bool[date] = self._check_for_holiday(date)
        self.load["weekday"] = self.load.index.dayofweek
        self.load["hour_of_day"] = self.load.index.hour
        self.load["temperature"] = self.load.index.map(pd.Series(weather_dict))
        self.load["holiday"] = self.load.index.map(holiday_bool)
        self.load["load (t-24)"] = self.load.load_MW.shift(24)

    def build_model_input(self):
        featurized_df = self.dummify_categorical_features(self.load.copy())
        self.model_input = featurized_df[featurized_df.notna()]

    @staticmethod
    def _set_iso(iso_name: str):
        if iso_name == "NYISO":
            iso_engine = client_factory("NYISO")
        elif iso_name == "ISONE":
            iso_engine = client_factory("ISONE", timeout_seconds=60)
        elif iso_name == "CAISO":
            iso_engine = client_factory("CAISO")
        elif iso_name == "ERCOT":
            iso_engine = EIAClient(timeout_seconds=60)
            iso_engine.BA = "ERCOT"
        elif iso_name == "PJM":
            iso_engine = EIAClient(timeout_seconds=60)
            iso_engine.BA = "PJM"
        elif iso_name == "MISO":
            iso_engine = EIAClient(timeout_seconds=60)
            iso_engine.BA = "MISO"
        else:
            print(f"Peaky Finders does not support {iso_name} yet!")
        return iso_engine

    def _get_temperature(self, date):
        date_input = date.strftime("%s")
        full_url = f"{self.weather_url}{date_input}?exclude={EXCLUDE}"
        response = requests.get(full_url)
        if response.status_code == 200:
            print(response.status_code)
        else:
            raise ValueError(
                f"Error getting data from DarkSky API: "
                f"Response Code {response.status_code}"
            )
        info = response.json()
        current_info = info["currently"]
        try:
            temp = current_info["temperature"]
        except KeyError:
            temp = None
        return temp

    @staticmethod
    def _check_for_holiday(day):
        if day in US_HOLIDAYS:
            return True
        else:
            return False

    @staticmethod
    def dummify_categorical_features(load_df: pd.DataFrame):
        for feature in CATEGORICAL_FEATURES:
            dummies = pd.get_dummies(load_df[feature], prefix=feature, drop_first=True)
            load_df = load_df.drop(feature, axis=1)
            load_df = pd.concat([load_df, dummies], axis=1)
            load_df = load_df.dropna(axis=0, how="any")
        return load_df
