library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(jsonlite)
library(dashBio)
library(GEOquery)


appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf("/app/apps/%s", appName))
}

app <- Dash$new()


source("assets/CircosFunctions.R")

json_file <- "assets/circos_graph_data.json"


circos_graph_data <- read_json(path = json_file)


applayout <- layout()



app$layout(htmlDiv(list(
  applayout
)))





app$callback(output = list( id= 'circos-hold',property = 'children'),
             params = list(input(id = 'circos-graph-type', property = 'value'),
                           input(id = 'circos-size', property = 'value')
             ),
             
             update_graph <- function(key, size) {
               return(get_circos_graph(key, size = size))
             }
)



app$callback(
  output(id = "chords-text", property = "children"),
  params = list(
    input(id = 'circos-graph-type', property = 'value')
  ),
  update_desc <- function(option) {
    if (option == 'select-dataset-chords') {
      return(htmlP("Highlight chords on the graph by selecting the 
                   'Table' tab and selecting rows in the 'Chords' dataset."))
    }
    
    else {
      return("    ")
    }
  }
)



if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server(debug = TRUE)
}
