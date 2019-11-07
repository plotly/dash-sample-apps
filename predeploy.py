import os
import shutil
import subprocess

from apps_directory_mapping import APPNAME_TO_DIRECTORY


pyfiles = ["requirements.txt", "Procfile"]

app_name = os.environ["DASH_APP_NAME"]

app_path = os.path.join("apps", APPNAME_TO_DIRECTORY.get(app_name, app_name))

app_file_name = ''

with open(os.path.join(app_path, 'Procfile'), 'r') as f:
    contents = f.read().split(' ')
    for item in contents:
        if 'server' in item:
            app_file_name = item.split(':')[0]

full_app_path = os.path.join(app_path, app_file_name + '.py')

lines = []

with open(full_app_path, 'r') as f:
    lines = f.readlines()
    name_main_index = 0
    for line in lines:
        if '__name__ == "__main__"' in line:
            name_main_index = lines.index(line)
            break
    lines.insert(name_main_index, '''app.index_string = \'\'\'
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        <div>My Custom header</div>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
        <div>My Custom footer</div>
    </body>
</html>
\'\'\'\n\n''')

if len(lines) > 0:
    with open(full_app_path, 'w') as f:
        f.writelines(lines)


if "DOKKU_SCALE" in os.listdir(app_path):
    shutil.copyfile(os.path.join(app_path, "DOKKU_SCALE"), "DOKKU_SCALE")

for f in pyfiles:
    shutil.copyfile(os.path.join(app_path, f), f)

subprocess.run("python -m pip install -r requirements.txt".split(" "))
