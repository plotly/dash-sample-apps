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

The data used in this demo app have been taken from the [MICCAI Brain Tumor Segmentation Challenge](http://braintumorsegmentation.org/) [1-3].

[1] B. H. Menze, A. Jakab, S. Bauer, J. Kalpathy-Cramer, K. Farahani, J. Kirby, et al. "The Multimodal Brain Tumor Image Segmentation Benchmark (BRATS)", IEEE Transactions on Medical Imaging 34(10), 1993-2024 (2015) DOI: 10.1109/TMI.2014.2377694

[2] S. Bakas, H. Akbari, A. Sotiras, M. Bilello, M. Rozycki, J.S. Kirby, et al., "Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features", Nature Scientific Data, 4:170117 (2017) DOI: 10.1038/sdata.2017.117

[3] S. Bakas, M. Reyes, A. Jakab, S. Bauer, M. Rempfler, A. Crimi, et al., "Identifying the Best Machine Learning Algorithms for Brain Tumor Segmentation, Progression Assessment, and Overall Survival Prediction in the BRATS Challenge", arXiv preprint arXiv:1811.02629 (2018)
