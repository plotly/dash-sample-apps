# Dash cuML UMAP

A demo of RAPIDS.ai's cuML UMAP functionality.

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