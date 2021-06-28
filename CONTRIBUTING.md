# CONTRIBUTING

## Running an example app after cloning the repo

You will need to run applications, and specify filenames, from the
root directory of the repository. e.g., if the name of the app you
want to run is `my_dash_app` and the app filename is `app.py`, you
would need to run `python apps/my_dash_app/app.py` from the root
of the repository.

Each app has a requirements.txt, install the dependencies in a virtual
environment.


## Contributing to the sample apps repo

_"DDS app" below refers to the deployed application. For example, if
the deployment will eventually be hosted at
https://dash-gallery.plotly.host/my-dash-app, "DDS app name" is
`my-dash-app`._

### Adding a new app

Create an app on Dash Playground. This will be the location of the
auto-deployment. To do this, log into the app manager on
[dash-playground.plotly.host](https://dash-playground.plotly.host)
and click "initialize app".

Create a branch from `master` to work on your app, the name is not required
to be anything specific. Switch to this branch, then navigate to the `apps/`
directory and add a directory for your app.

There are two options when you are naming the folder:

1. Make the folder have the _exact same_ name as the Dash app name.

2. (Python apps only) Select any other name, but _update the file
   [`apps_mapping.py`](apps_directory_mapping.py)_ with the Dash app
   name and the folder name you have selected.

Navigate to the directory you just created, and write a small README
that only contains the name of the app. Stage the README and commit it
to your branch.

See [project boilerplate!](https://github.com/plotly/dash-sample-apps#project-boilerplate)

### Notes on adding a new Dash for R app

Contributing an app written with Dash for R is very similar to the steps outlined above. 

1. Make the folder have the _exact same_ name as the Dash app name.

2. Ensure that the file containing your app code is named `app.R`.

3. The `Procfile` should contain 

```
web: R -f /app/app.R
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

### Making changes to an existing app

Create a new branch - of any name - for your code changes.
Then, navigate to the directory that has the same name as
the DDS app.

When you are finished, make a pull request from your branch to the master
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
    - Ex. `web: gunicorn app:server`
- **`requirements.txt`**
    - Install project dependecies in a virtual environment

#### Project boilerplate

    apps
    ├── ...
    ├── {DASH_APP_NAME}         # app project level
    │   ├── assets/             # all stylesheets and javascript files
    │   ├── data/               # all data (csv, json, txt, etc)
    │   ├── app.py              # dash application entry point
    │   ├── Procfile            # used for heroku deployment (how to run app)
    │   ├── requirements.txt    # project dependecies
    │   └── ...                 
    └── ...

#### Handle relative path

Assets should never use a relative path, as this will fail when deployed to Dash Enterprise due to use of subdirectories for serving apps.

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
git checkout -b "{YOUR_CUSTOM_BRANCH}"

# create a new folder in apps/
mkdir /apps/{DASH_APP_NAME}

# push new branch
git push -u origin {YOUR_CUSTOM_BRANCH}
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
# once your branch is ready, make a PR into master!

PR has two checkers.
1. make sure your code passed the black linter
2. make sure your project is deployed on dns playground
```

