import os
import shutil
import subprocess

from apps_directory_mapping import APPNAME_TO_DIRECTORY


files = ["requirements.txt", "Procfile"]

app_path = os.path.join(
    "apps",
    APPNAME_TO_DIRECTORY.get(os.environ["DASH_APP_NAME"], os.environ["DASH_APP_NAME"]),
)

for f in files:
    shutil.copyfile(os.path.join(app_path, f), f)

if "DOKKU_SCALE" in os.listdir(app_path):
    shutil.copyfile(os.path.join(app_path, "DOKKU_SCALE"), "DOKKU_SCALE")

subprocess.run("python -m pip install -r requirements.txt".split(" "))
