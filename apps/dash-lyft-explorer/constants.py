from lyft_dataset_sdk.lyftdataset import LyftDataset
import colorlover as cl
from utils.helper_functions import get_token_list

## Define constants used across the app
CAMERAS = [
    "CAM_FRONT",
    "CAM_BACK",
    "CAM_FRONT_ZOOMED",
    "CAM_FRONT_LEFT",
    "CAM_FRONT_RIGHT",
    "CAM_BACK_RIGHT",
    "CAM_BACK_LEFT",
]
LIDARS = ["LIDAR_TOP", "LIDAR_FRONT_RIGHT", "LIDAR_FRONT_LEFT"]

NAME2COLOR = dict(
    zip(
        ["bus", "car", "other_vehicle", "pedestrian", "truck"],
        cl.to_numeric(cl.scales["5"]["div"]["Spectral"]),
    )
)

## Create Lyft object
lv5 = LyftDataset(data_path="./data", json_path="./data/train_data", verbose=True)

# Load a single scene
scene = lv5.scene[0]
token_list = get_token_list(scene, lv5)
INITIAL_TOKEN = scene["first_sample_token"]

