import os
import shutil
import subprocess

files = ["requirements.txt", "Procfile", "DOKKU_SCALE"]
for f in files:
    shutil.copyfile(os.path.join("apps", os.environ["DASH_APP_NAME"], f), f)

subprocess.run("python -m pip install -r requirements.txt".split(" "))
