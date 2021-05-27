# LAStoDash

## About this app

LAStoDash is a sample Dash project that takes a [Log ASCII Standard (LAS) file](http://www.cwls.org/las/) and builds a web app to view its content and print in PDF format.

As indicated in the [LAS 2.0 Specifications](http://www.cwls.org/wp-content/uploads/2017/02/Las2_Update_Feb2017.pdf), LAS files contain sections that are marked by a ~. The [LAS file](data/alcor2.las) used for this app contains four of these sections, including version and wrap mode information (~V), well identification (~W), curve information (~C), and ASCII log data (~A), all of which are displayed in some format (e.g. graph, table) in the demo app and the [printable report in PDF format](demo/alcor2.pdf).

## How to run this app locally

Clone the repository:

```
$ git clone https://github.com/plotly/dash-sample-apps.git
```

Redirect to the respective app directory:

```
$ cd dash-sample-apps/apps/dashr-lastodash
```

Install the requirements:

```
$ Rscript init.R
```

Run the app:

```
$ Rscript app.R
```

View in your browser at http://127.0.0.1:8050.

## Screenshots

![demo/demo.gif](demo/demo.gif)

## Resources

* [Dash documentation for R](https://dashr-docs.herokuapp.com/)
* [Log ASCII Standard (LAS) file](http://www.cwls.org/las/)
