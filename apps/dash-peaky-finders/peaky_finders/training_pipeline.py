import argparse
import os
import pickle

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error

from peaky_finders.data_acquisition.train_model import LoadCollector

"""
CLI demo command:
python peaky_finders/training_pipeline.py --iso NYISO --model xgboost --start_date 01-01-2019 --end_date 07-28-2020 --save_model_input True --save_model_output True
"""


MODEL_INPUT_DIR = os.path.join(os.path.dirname(__file__), "training_data")
MODEL_OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "models")
"""Directory to hold all pipeline output csv files."""


class Pipeline:
    def __init__(
        self,
        iso: str,
        model: str,
        start_date: str,
        end_date: str,
        save_model_input: bool,
        save_model_output: bool,
    ) -> None:
        """
        Args:
            iso: the iso to forecast ('nyiso', 'ercot', etc.)
            model: logistic regression 'log' or xgboost
        """
        self.iso_name = iso
        self.model = model
        self.start = start_date
        self.end = end_date
        self.iso = LoadCollector(iso, start_date, end_date)
        self.save_model_input = save_model_input
        self.save_model_output = save_model_output

    def phase_one(self):
        """Feature engineering + model preparation """
        self.iso.engineer_features()
        self.iso.build_model_input()
        if self.save_model_input:
            self.iso.model_input.to_csv(self.model_input_filepath)

    def phase_two(self):
        """Model training & serialization"""
        y = self.iso.model_input["load_MW"]
        X = self.iso.model_input.drop("load_MW", axis=1)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
        reg = xgb.XGBRegressor()
        reg.fit(X_train, y_train)
        training_preds = reg.predict(X_train)
        val_preds = reg.predict(X_test)
        print("Mean Absolute Error:", mean_absolute_error(y_test, val_preds))
        print(
            "Root Mean Squared Error:", np.sqrt(mean_squared_error(y_test, val_preds))
        )
        if self.save_model_output:
            pickle.dump(reg, open(self.model_output_filepath, "wb"))

    def execute(self):
        self.phase_one()
        self.phase_two()

    @property
    def model_input_filepath(self):
        """
        Creates a dew filename depending on data version and feeder selected.
        """
        return os.path.join(
            MODEL_INPUT_DIR, (f"{self.iso_name}_{self.start}_{self.end}.csv")
        )

    @property
    def model_output_filepath(self):
        """
        Creates a dew filename depending on data version and feeder selected.
        """
        return os.path.join(
            MODEL_OUTPUT_DIR, (f"xg_boost_{self.iso_name}_load_model.pkl")
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--iso",
        default="NYISO",
        type=str,
        help="Select ISO for model prediction.",
    )
    parser.add_argument(
        "-m",
        "--model",
        type=str,
        help="Model to train (xgboost for forecast, log regression for peak.",
    )
    parser.add_argument(
        "-s", "--start_date", type=str, help="Start date range for model training."
    )
    parser.add_argument(
        "-e", "--end_date", type=str, help="End date range for model training."
    )
    parser.add_argument(
        "-mi",
        "--save_model_input",
        type=bool,
        help="Save model input (w/ features and scaled.",
    )
    parser.add_argument(
        "-mo",
        "--save_model_output",
        type=bool,
        help="Save trained model as pickle file.",
    )

    args = parser.parse_args()
    pipeline = Pipeline(
        iso=args.iso,
        model=args.model,
        start_date=args.start_date,
        end_date=args.end_date,
        save_model_input=args.save_model_input,
        save_model_output=args.save_model_output,
    )
    pipeline.execute()
