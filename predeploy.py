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

if "-" in app_name and app_name.split("-")[0] == "dashr":
    for rfile in rfiles:
        shutil.copyfile(os.path.join(app_path, rfile), rfile)
else:
    for f in pyfiles:
        shutil.copyfile(os.path.join(app_path, f), f)
    lines = []
    ga_line = 0
    with open(os.path.join(app_path, "app.py"), "r") as f:
        lines = f.readlines()

    try:
        ga_line = [i for i in range(len(lines)) if "app.layout" in lines[i]]
        lines.insert(
            ga_line[0],
            "app.index_string = '''\n{}\n'''\n\n".format(os.environ["DASH_GA_CODE"]),
        )
    except IndexError as e:
        pass

    with open(os.path.join(app_path, "app.py~"), "w+") as f:
        f.write("".join(lines))

    shutil.copyfile(os.path.join(app_path, "app.py~"), os.path.join(app_path, "app.py"))

subprocess.run("python -m pip install -r requirements.txt".split(" "))
