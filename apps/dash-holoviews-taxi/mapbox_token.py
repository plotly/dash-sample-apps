import os


def get_mapbox_token():
    if os.environ.get("MAPBOX_TOKEN", None):
        return os.environ["MAPBOX_TOKEN"]
    elif os.path.exists(".mapbox_token"):
        with open(".mapbox_token", 'rt') as f:
            return f.read().strip()
    else:
        raise ValueError(
            "Could not load mapbox token\n"
            "Either set the MAPBOX_TOKEN environment variable to your token string \n"
            "or create a file named .mapbox_token containing your token in the\n"
            "top-level project directory"
        )
