# Dash cuML Regression Demo

This app shows how to build a simple Dash app that lets you train super fast ML models using Rapids.ai's cuML library.

![demo](demo.gif)

## Instructions

To get started, first clone this repo:
```
git clone https://github.com/plotly/dash-sample-apps.git
cd dash-sample-apps/apps/dash-cuml-umap
```

Create a conda env and install the requirements:
```
conda create --name dash-cuml-umap --file requirements.txt
```

You can now run the app:
```
python app.py
```

and visit http://127.0.0.1:8050/.


## Note about exporting

Since RAPIDS.ai only support a conda environment, `pip freeze` the requirements would not work. Instead, do:
```
conda list -e > requirements.txt
```