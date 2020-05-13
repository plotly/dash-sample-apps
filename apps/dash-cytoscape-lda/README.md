# Dash-cytoscape NLP demo for Plotly

This demo app is built using plotly-dash and dash-cytoscape, as an 
example of a network diagram that may be displayed.

The dataset used it the 'Commercial use subset' of the CORD-19 dataset from 
[Semantic Scholar](https://pages.semanticscholar.org/coronavirus-research), 
downloaded on 2/Apr/2020. 

Simply run `python app.py` to see the network visualisation as below.

![app screenshot](screenshot1.png "Screenshot of the app")

The scripts used to produce the data are labelled by number:
1_lda_analysis.py
2_add_biblio_data.py
3_reduce_num_nodes.py

Simply downloading the dataset (including metadata.csv) 
and running these scripts should reproduce the same outputs.

* Created by JP Hwang for Plotly 
