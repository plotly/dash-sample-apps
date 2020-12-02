import os


def get_mapbox_token():
    if os.environ.get("MAPBOX_TOKEN", None):
        return os.environ["MAPBOX_TOKEN"]
    else:
        return "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNrOWJqb2F4djBnMjEzbG50amg0dnJieG4ifQ.Zme1-Uzoi75IaFbieBDl3A"
