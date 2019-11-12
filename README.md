# Dash Sample Apps

This is a monorepo designed to host all of the apps that have been
created for the Python Dash Gallery.

## Running an example app

You will need to run applications, and specify filenames, from the
root directory of the repository. e.g., if the name of the app you
want to run is `my_dash_app` and the app filename is `app.py`, you
would need to run `python apps/my_dash_app/app.py` from the root
of the repository.

Each app has a requirements.txt, install the dependecies in a virtual
environment.

## Contributing to the sample apps repo

_"DDS app" below refers to the deployed application. For example, if
the deployment will eventually be hosted at
https://dash-gallery.plotly.host/my-dash-app, "DDS app name" is
`my-dash-app`._

### Branches

Each app has its own branch off of `master` that has the _exact same_
name as the DDS app. This is an effective `master` branch _for that
app only_. This is because we sync the apps in this repository with
our staging deployment server, and the automatic deploys sync the app
name with the github branch name. So, for automatic deploys to work,
any changes for a particular app should be done on a branch that has
the same name as the app.

### Adding a new app

Create an app on Dash Playground. This will be the location of the
auto-deployment. To do this, log into the app manager on
[dash-playground.plotly.host](https://dash-playground.plotly.host)
and click "initialize app".

Create a branch from `master` that has the _exact same_ name as the
Dash app name. Switch to this branch, then navigate to the `apps/`
directory and add a directory for your app.

There are two options when you are naming the folder:

1. Make the folder have the _exact same_ name as the Dash app name.

2. (Python apps only) Select any other name, but _update the file
   [`apps_mapping.py`](apps_directory_mapping.py)_ with the Dash app
   name and the folder name you have selected.

Navigate to the directory you just created, and write a small README
that only contains the name of the app. Stage the README and commit it
to your app branch.

See [project boilerplate!](https://github.com/plotly/dash-sample-apps#project-boilerplate)

### Notes on adding a new Dash for R app

Contributing an app written with Dash for R is very similar to the steps outlined above. 

1. Make the folder have the _exact same_ name as the Dash app name.

2. Ensure that the file containing your app code is named `app.R`.

3. The `Procfile` should contain 

```
web: R -f /app/apps/"$DASH_APP_NAME"/app.R
```

4. Routing and request pathname prefixes should be set. One approach might be to include

```
appName <- Sys.getenv("DASH_APP_NAME")
pathPrefix <- sprintf("/%s/", appName)

Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
           DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
```

at the head of your `app.R` file.

5. `run_server()` should be provided the host and port information explicitly, e.g.

``
app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
``

6. For convenience, it is probably easiest to set the working directory in `app.R` as well:

``
setwd(sprintf("/app/apps/%s", appName))
``

### Making changes to an existing app

Switch to the branch that has the same name as the DDS app (the "app
branch"). Then, navigate to the directory that has the same name as
the DDS app.

When you are finished, make a pull request from the app branch to the master
branch. Once you have passed your code review, you can merge your PR.

## Dash app project structure

#### Data
- All data (csv, json, txt, etc) should be in a data folder
- `/apps/{DASH_APP_NAME}/data/`

#### Assets
- All stylesheets and javascript should be in an assets folder
- `/apps/{DASH_APP_NAME}/assets/`

####  These files will still need to be present within the app folder.

- **`Procfile`** gets run at root level for deployment
    - Make sure python working directory is at the app level
    - Ex. `web: gunicorn --pythonpath apps/{DASH_APP_NAME} app:server`
- **`requirements.txt`**
    - Install project dependecies in a virtual environment
- **`runtime.txt`**
    - App python version

#### Project boilerplate

    apps
    ├── ...
    ├── {DASH_APP_NAME}         # app project level
    │   ├── assets/             # all stylesheets and javascript files
    │   ├── data/               # all data (csv, json, txt, etc)
    │   ├── app.py              # dash application entry point
    │   ├── Procfile            # used for heroku deployment (how to run app)
    │   ├── requirements.txt    # project dependecies
    │   ├── runtime.txt         # used for heroku deployment (python version)
    │   └── ...                 
    └── ...

#### Handle relative path

Since deployment happens at the root level `/` and not at the app level (`/apps/{DASH_APP_NAME}`), we need to make sure our application is able to run at both levels for flexibility.

Reading from assets and data folder
```Python
Img(src="./assets/logo.png") will fail at root level
```

Tips
 
-  Use [get_asset_url()](https://dash.plot.ly/dash-deployment-server/static-assets)
-  Use [Pathlib](https://docs.python.org/3/library/pathlib.html) for more flexibility

```Python
import pathlib
import pandas as pd

# get relative assets
html.Img(src=app.get_asset_url('logo.png'))       # /assets/logo.png

# get relative data

DATA_PATH = pathlib.Path(__file__).parent.joinpath("data")  # /data
df = pd.read_csv(DATA_PATH.joinpath("sample-data.csv"))    # /data/sample-data.csv

with open(DATA_PATH.joinpath("sample-data.csv")) as f:  # /data/sample-data.csv
    some_string = f.read()
```

## Developer Guide

#### Creating a new project

```
# branch off master
git checkout -b "{DASH_APP_NAME}"

# create a new folder in apps/
mkdir /apps/{DASH_APP_NAME}

# push new app branch
git push -u origin {DASH_APP_NAME}
```

#### Before committing

```
# make sure your code is linted (we use black)
black . --exclude=venv/ --check

# if black is not installed
pip install black
```


#### App is ready to go!
```
# once your app branch is ready, make a PR into master!

PR has two checkers.
1. make sure your code passed the black linter
2. make sure your project is deployed on dns playground
```

