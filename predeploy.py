import os
import shutil
import subprocess

from apps_directory_mapping import APPNAME_TO_DIRECTORY


pyfiles = ["requirements.txt", "Procfile"]
rfiles = ["Aptfile", "init.R", "Procfile"]

app_name = os.environ["DASH_APP_NAME"]

app_path = os.path.join("apps", APPNAME_TO_DIRECTORY.get(app_name, app_name))

if "DOKKU_SCALE" in os.listdir(app_path):
    shutil.copyfile(os.path.join(app_path, "DOKKU_SCALE"), "DOKKU_SCALE")

try:
    with open(
        os.path.join(os.path.join(app_path, "assets"), "plotly_ga.js"), "w+"
    ) as f:
        f.write(os.environ["PLOTLY_GA_CODE"])
except FileNotFoundError:
    print("No assets/ folder found.")
    raise SystemExit

file_list = rfiles if "-" in app_name and app_name.split("-")[0] == "dashr" else pyfiles

for f in file_list:
    shutil.copyfile(os.path.join(app_path, file_list), file_list)

subprocess.run("python -m pip install -r requirements.txt".split(" "))
