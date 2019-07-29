# DashR - Circos

## About this app:

Circos is a circular visualization of data, and can be used to highlight relationships between objects in a dataset (e.g., genes that are located on different chromosomes in the genome of an organism).

A Dash Circos graph consists of two main parts: the layout and the tracks. The layout sets the basic parameters of the graph, such as radius, ticks, labels, etc; the tracks are graph layouts that take in a series of data points to display.

The visualizations supported by Dash Circos are: heatmaps, chords, highlights, histograms, line, scatter, stack, and text graphs.

### Using the demo

In the "Data" tab, you can opt to use preloaded datasets; additionally, you can download sample data that you would use with a Dash Circos component, upload that sample data, and render it with the "Render" button.

In the "Graph" tab, you can choose the type of Circos graph to display, control the size of the graph, and access data that are generated upon hovering over parts of the graph.

In the "Table" tab, you can view the datasets that define the parameters of the graph, such as the layout, the highlights, and the chords. You can interact with Circos through this table by selecting the "Chords" graph in the "Graph" tab, then viewing the "Chords" dataset in the "Table" tab.
Reference: [Seminal Paper](http://www.doi.org/10.1101/gr.092759.109)
For a look into Circos and the Circos API, please visit the original repository [here](https://github.com/nicgirault/circosJS).

![screenshot](assets/dashr-circos-screenshot.png)

### Running the app locally

Clone the git repo.

Run `app.R` within the repo working directory.

Open a browser at http://127.0.0.1:8050
