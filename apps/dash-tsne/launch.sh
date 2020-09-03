wget -N https://github.com/plotly/datasets/raw/master/dash-sample-apps/dash-tsne/data.zip
wget -N https://github.com/plotly/datasets/raw/master/dash-sample-apps/dash-tsne/demo_embeddings.zip
unzip -nq ./data.zip
unzip -nq ./demo_embeddings.zip
gunicorn app:server --workers 4