wget https://github.com/plotly/datasets/releases/download/v0.0.1-data/dash-tsne.zip
unzip dash-tsne.zip
gunicorn app:server --workers 4