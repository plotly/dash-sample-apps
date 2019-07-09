# dashr-moledule3d
# Moledule3d

This is a dash for R version of the moledule3d app written in python [moledule3d](https://github.com/plotly/dash-bio/blob/master/tests/dashbio_demos/app_molecule3d.py)

## Screenshots
![assets/molecule3d.gif](assets/molecule3d.gif)

## About this app:

Molecule3D is a visualizer that allows you to view biomolecules in multiple representations: sticks, spheres, and cartoons. "pdb" files 
are required for the database.

### Using the demo

#### Running the app locally

Clone the git repo and change to the root directory 

```
git clone https://github.com/plotly/dash-sample-apps
cd dash-sample-apps/apps/dashr-molecule3d
```
Install the requirements. From the terminal, run the following to install the required packages in the default location:

```
remotes::install_github("plotly/dashR", ref="0.1.0-cran")
remotes::install_github("plotly/dash-html-components", ref="1.0.0-cran")
remotes::install_github("plotly/dash-core-components", ref="1.0.0-cran")
remotes::install_github("plotly/dash-daq", ref="update-for-dash-0.1.0")
remotes::install_github('plotly/dash-bio', ref="update-for-dash-0.1.0")

install.packages("bio3d")
install.packages("jsonlite")
```

Run the app. From the terminal, run:

```
Rscript app.R
```

Open a browser at http://127.0.0.1:8050
