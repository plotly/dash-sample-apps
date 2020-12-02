import os
import requests


def download_dataset():
    # Download and cache dataset
    data_path = "./data/nyc_taxi_small.parq"
    if not os.path.exists(data_path):
        print("Downloading Dataset")
        os.makedirs("./data", exist_ok=True)
        response = requests.get(
            "https://github.com/plotly/dash-holoviews-taxi/releases/download/v0.0.1a1/nyc_taxi_small.parq"
        )
        with open(data_path, "wb") as f:
            f.write(response.content)
    else:
        print("Found cached dataset")

    return data_path


if __name__ == "__main__":
    download_dataset()
