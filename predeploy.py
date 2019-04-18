import os
import shutil
import subprocess

from apps_directory_mapping import APPNAME_TO_DIRECTORY


files = ["requirements.txt", "Procfile", "DOKKU_SCALE"]
for f in files:
    shutil.copyfile(
        os.path.join(
            "apps",
            APPNAME_TO_DIRECTORY.get(
                os.environ["DASH_APP_NAME"], os.environ["DASH_APP_NAME"]
            ),
            f,
        ),
        f,
    )

subprocess.run("python -m pip install -r requirements.txt".split(" "))
