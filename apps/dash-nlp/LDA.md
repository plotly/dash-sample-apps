## LDA precomputed data for the whole dataset

All data needed for displaying LDA plot and table is located in [this json file](https://github.com/vildly/dash-sample-apps/blob/master/apps/dash-nlp/data/precomputed.json) in "data" folder of this app. 

For each bank it contains the following data in json format:

* Dataframe  "tsne_df" 
Main dataframe containing results of LDA analysis. LDA dimensionality is reduced using t-SNE and clusters are visiualised only using two first principal components. Hense the saved data: 

Columns: 
"tsne_x" - PC0

"tsne_y" - PC1

"topic_num" - topic id

"doc_num" -  df_dominant_topic["Document_No"], document id
                

* Dataframe "df_dominant_topic"

Additional data of LDA analysis. Columns are:
  "Document_No",
  
  "Dominant_Topic",
  
  "Topic_Perc_Contrib",
  
  "Keywords",
  
  "Text",
  
  "Date"

This df is mainly used for displaying data in the table and contains original text of the complaint ("Text") as well as mappings to document id (  "Document_No" ) and topic id ("Dominant_Topic")

* Dataframe "df_top3words" 

3 most popular words per topic found in this bank complaints. Columns: 
columns are "topic_id" and "words". Helper dataframe for diplaying some information about the topic
without displaying all the corresponding Keywords. 


## How LDA was precomputed? 
See [this python script](https://github.com/plotly/dash-sample-apps/blob/master/apps/dash-nlp/precomputing.py) to learn how precomputed data was generated 

## How LDA is visualised 
LDA plot is a scatter with a separate trace for each cluster. For each trace:

x=tsne_df_f["tsne_x"],

y=tsne_df_f["tsne_y"],

mode="markers",

hovertext=tsne_df_f["doc_num"],

marker=dict(
    size=6,
    color=mycolors[tsne_df_f["topic_num"]],  # set color equal to a variable
    colorscale="Viridis",
    showscale=False,
),
)

Color of the cluster is based on topic id. 

LDA table is populated when bank is selected, but all records are hidden until user clicks on selected document in the scatter plot. See filter_table_on_scatter_click method for the click callback. Table is implemented using dash table component.

Filter query used:

filter_query = (
    "({Document_No} eq "
    + str(selected_complaint)
    + ") || ("
    + current_filter
    + ")"
)

More on filtering syntax in R: https://dashr.plot.ly/datatable/filtering

       
