# Dash Sample Apps

This is a monorepo designed to host all of the apps that have been
created for the Python Dash Gallery.

## Running an example app

You will need to run applications, and specify filenames, from the
root directory of the repository. e.g., if the name of the app you
want to run is `my_dash_app` and the app filename is `app.py`, you
would need to run `python my_dash_app/app.py` from the root of the
repository.

## Contributing to the sample apps repo

_"DDS app" below refers to the deployed application. For example, if
the deployment will eventually be hosted at
https://dash-gallery.plotly.host/my-dash-app, "DDS app name" is
`my-dash-app`._

### Branches

Each app has its own branch off of `master` that has the _exact same_
name as the DDS app. This is an effective `master` branch _for that
app only_.

### Adding a new app

Create a branch from `master` that has the _exact same_ name as the
Dash app name. Switch to this branch, then navigate to the `apps/`
directory and add a directory for your app.

There are two options when you are naming the folder:

1. Make the folder have the _exact same_ name as the Dash app name.

2. Select any other name, but _update the file [`apps_mapping.py`](#)_
   with the Dash app name and the folder name you have selected.

Navigate to the directory you just created, and write a small README
that only contains the name of the app. Stage the README and commit it
to your app branch.

### Making changes to an existing app

Switch to the branch that has the same name as the DDS app (the "app
branch"). Then, navigate to the directory that has the same name as
the DDS app.

Create a branch off of the app branch, and do all of your development
in this branch (the "feature branch"). When you are finished, make a
pull request from the feature branch to the app branch. Once you have
passed your code review, you can merge your PR.

*Merges from the app branch to the `master` branch of the repository
 are handled automatically, and will only happen if the code in the
 app branch passes the tests run by `black`.*


#### Reading files from within an app

File paths must be specified relative to the root directory. So, for
example, if you have a file 'data.csv' in your app folder
`my_dash_app/` you will have to refer to this data file as
`apps/my_dash_app/data.csv` within your `app.py` file.

#### `assets/`, `Procfile`, `requirements.txt`, `runtime.txt`

These files will still need to be present within the app folder.
