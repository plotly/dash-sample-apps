import os
import shutil
import subprocess

from apps_directory_mapping import APPNAME_TO_DIRECTORY


pyfiles = ["requirements.txt", "Procfile"]

app_name = os.environ["DASH_APP_NAME"]

app_path = os.path.join("apps", APPNAME_TO_DIRECTORY.get(app_name, app_name))

if "DOKKU_SCALE" in os.listdir(app_path):
    shutil.copyfile(os.path.join(app_path, "DOKKU_SCALE"), "DOKKU_SCALE")

for f in pyfiles:
    shutil.copyfile(os.path.join(app_path, f), f)

subprocess.run("python -m pip install -r requirements.txt".split(" "))
