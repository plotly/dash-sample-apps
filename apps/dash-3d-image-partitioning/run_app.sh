gunicorn \
--pythonpath plotly-common \
-e LOAD_SUPERPIXEL=assets/BraTS19_2013_10_1_flair_superpixels.npz.gz \
app:server
