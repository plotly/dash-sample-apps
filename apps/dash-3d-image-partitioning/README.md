# Partitioning 3D images with `Dash` and `scikit-image`

## Setup

run `./setup` to create the virtualenv and download the necessary libraries.
This will also set up the git submodules.

To run programs in this directory, do `source venv/bin/activate` after running
`./setup`.

## Running

A recommended command line is

```bash
LOAD_SUPERPIXEL=assets/BraTS19_2013_10_1_flair_superpixels.npz.gz \
PYTHONPATH=plotly-common \
python app.py
```

Then you can navigate to the displayed link in your browser. You can also run
without specifying `LOAD_SUPERPIXEL`, in which case the segmentation will happen
when the app loads.

To generate a new superpixel file, you can specify a path to the environment
variable `SAVE_SUPERPIXEL` in which case the app will run, compute the
superpixels, save them, and exit.

## Resources

To learn more about [Dash](https://plot.ly/dash).
To learn more about [scikit-image](https://scikit-image.org/).
