import os
import shutil
import subprocess

from apps_directory_mapping import APPNAME_TO_DIRECTORY


files = ["requirements.txt", "Procfile"]
rfiles = ["Aptfile", "init.R"]

app_name = os.environ["DASH_APP_NAME"]

app_path = os.path.join("apps", APPNAME_TO_DIRECTORY.get(app_name, app_name))

for f in files:
    shutil.copyfile(os.path.join(app_path, f), f)

if "DOKKU_SCALE" in os.listdir(app_path):
    shutil.copyfile(os.path.join(app_path, "DOKKU_SCALE"), "DOKKU_SCALE")

if "-" in app_name and app_name.split("-")[0] == "dashr":
    for rfile in rfiles:
        shutil.copyfile(os.path.join(app_path, rfile), rfile)

subprocess.run("python -m pip install -r requirements.txt".split(" "))
