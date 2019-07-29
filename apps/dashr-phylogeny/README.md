# Phylogeny

## About this app

This demo lets you interactively explore the phylogeny and spread of selected diseases. 

## Use this app

Select a disease from the "Dataset" dropdown, and the date range you want from the "Data Range" slider, and it will display:

- A line plot of the number of cases per country, per year
- The phylogenetic tree of the disease
- A map showing the number of cases per country
- A histogram of the distribution of cases by country

![demo/demo.gif](demo/demo.gif)

## Run this app locally

Clone the repository:

```
$ git clone https://github.com/plotly/dash-sample-apps.git
```

Redirect to the respective app directory:

```
$ cd dash-sample-apps/apps/dashr-phylogeny
```

Install the requirements:

```
$ Rscript init.R
```

Run the app:

```
$ Rscript app.R
```

Open a browser at http://127.0.0.1:8050.

## Resources/References

* [Dash R](https://dashr.plot.ly/) - Main server and interactive components 
* [Plotly R](https://plot.ly/r/) - Used to create the interactive plots
