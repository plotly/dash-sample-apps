from os import listdir, environ
from os.path import isdir, join
import shutil
import subprocess

app_path = "apps"

# Read apps directory
app_dir_names = [f for f in listdir(app_path) if isdir(join(app_path, f))]

# DDS provides DASH_APP_NAME variable
app_name = environ["DASH_APP_NAME"]

# It should match one of the directory names
# We assume directory name == App Name
if app_name not in app_dir_names:
    raise Exception(
        "App name {} not found in {} directory - the names must match".format(
            app_name, app_path
        )
    )

files = ["requirements.txt", "Procfile", "DOKKU_SCALE"]
for f in files:
    shutil.copyfile(join("apps", app_name, f), f)

subprocess.run("python -m pip install -r requirements.txt".split(" "))
