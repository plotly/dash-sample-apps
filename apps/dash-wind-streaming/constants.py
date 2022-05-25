import os

app_color = {"graph_bg": "#082255", "graph_line": "#007ACE"}

DB_FILE = "data/wind-data.db"

GRAPH_INTERVAL = int(os.environ.get("GRAPH_INTERVAL", 5000))