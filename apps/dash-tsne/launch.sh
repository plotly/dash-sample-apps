wget https://github.com/plotly/datasets/releases/download/dash-sample-apps/dash-tsne.zip
unzip dash-tsne.zip
gunicorn app:server --workers 4