import os
import shutil
import subprocess

APPNAME_TO_DIRECTORY = {
    'lastodash': 'lastodash',
}

files = ['requirements.txt', 'Procfile', 'DOKKU_SCALE']
for f in files:
    shutil.copyfile(
        os.path.join(
            'apps',
            APPNAME_TO_DIRECTORY[os.environ['DASH_APP_NAME']],
            f
        ),
        f
    )

subprocess.run('python -m pip install -r requirements.txt'.split(' '))
