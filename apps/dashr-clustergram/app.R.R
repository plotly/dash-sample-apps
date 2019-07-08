library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashBio)
library(heatmaply)
library(data.table)
library(dashDaq)



appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf("/app/apps/%s", appName))
}


source("utils/dash_bio_utils.R")



app <- Dash$new()

# Load the datasets.

iris <- read.table(file = "assets/sample_data/clustergram_iris.tsv", sep = 
                     "\t", skip = 4, row.names = 1,  header = TRUE)

mtcars<- read.table(file = "assets/sample_data/clustergram_mtcars.tsv", sep = 
                      "\t", skip = 4, row.names = 1, header = TRUE)

datasets = list("iris" = iris, "mtcars" = mtcars)


options_tabs <- htmlDiv(id = 'clustergram-body', className = 'app-body', children = list(
  dccLoading(className = 'dashbio-loading', children = htmlDiv(
    id = 'clustergram-wrapper',
    children = dccGraph(id = 'clustergram', style = list('display' = 'none'))
  )),
  
  htmlDiv(id = 'clustergram-control-tabs', className = 'control-tabs', children = list(
    dccTabs(id = 'clustergram-tabs', value = 'what-is', children = list(
      dccTab(
        label = 'About',
        value = 'what-is',
        children = htmlDiv(className = 'control-tab', children = list(
          htmlH4(className = 'what-is', children = "What is Clustergram?"),
          htmlP('Clusterman is a combination of a heatmap and dendrograms that
                allows you to display hierarchical clustering data. Clusters on
                the dendrogram are highlighted in one color if they comprise
                 data points that share some minimal level of correlation.'),
          htmlP('In the "Data" tab, you can select a preloaded dataset or, 
                alternatively, upload one of your own. A sample dataset is also
                available for download in the tab.'),
          htmlP('In the "Graph" tab, you can choose the dimension(s) along which
                clustering will be performed (row or column). You can also change
                the threshold that determines the point at which clusters are 
                highlighted for the tow and column dendrograms, and choose which rows
                and columns are used to compute the clustering.'),
          htmlP('In addition, you can highlight specific clusters by adding annotations
                to the clustergram, and choose whether to show or hide the labels for
                the rows and/or columns.')
        ))
      ),
      dccTab(
        label = 'Data',
        value = 'datasets',
        children = htmlDiv(className = 'control-tab', children = list(
          
          htmlDiv(
            'Preloaded dataset',
            title = 'Choose from some pre-loaded datasets to
            view them on the heatmap.',
            className = 'fullwidth-app-controls-name'
          ),
          
          dccDropdown(
            id = 'clustergram-datasets',
            options = list(
              list('label' = 'Andersons Iris Data', 'value' = 'iris'),
              list('label' = 'mtcars', 'value' = 'mtcars'),
              list('label' = 'Prostate Cancer', 'value' = 'prostatecancer'),
              list('label' = 'Lung Cancer Subtypes', 'value' = 'lungcancer')
            ),
            value = 'prostatecancer'
          ),
          
          htmlBr(),
          
          htmlDiv(
            'Upload dataset',
            title = 'Upload your own dataset below',
            className = 'app-controls-name'
          ),
          
          htmlDiv(
            id = 'file-upload-name'
          ),
          
          htmlDiv(
            id = 'clustergram-file-upload-container',
            title = 'Upload your own dataset here.',
            children = list(
              dccUpload(
                id = 'file-upload',
                className = 'control-upload',
                children = htmlDiv(list(
                  "Drag and drop .soft files, or click to select files."
                )),
                accept = '.soft'
              )
            )
          ),
          
          htmlDiv(list(
            htmlA(
              htmlButton(
                'Download sample .soft data',
                id = 'clustergram-download-sample-data',
                className = 'control-download'
              ),
              href = ("assets/sample_data/clustergram_GDS5826.soft"),
              download = 'clustergram_GDS5826.soft'
            )
          )),
          htmlHr(),
          htmlDiv(
            id = 'clustergram-info'
          )
        ))
      ),
      
      dccTab(
        label = 'Graph',
        value = 'graph',
        children = list(htmlDiv(className = 'control-tab', children = list(
          htmlDiv(className = 'app-controls-block', children= list(
            htmlDiv(
              'Cluster by:',
              title = 'Calculate dendrogram for row data, column data, or both.',
              className = 'app-controls-name'
            ),
            dccDropdown(
              id = 'cluster-checklist',
              options = list(
                list('label' = 'Row', 'value' = 'row'),
                list('label' = 'Column', 'value' = 'col')
              ),
              value = list('row', 'col'),
              multi = TRUE
            )
          )),
          
          htmlDiv(className = 'app-controls-block', children = list(
            htmlDiv(
              'Hide labels',
              title = 'Hide labels for the row and/or column dendrograms.',
              className = 'app-controls-name'
            ),
            dccDropdown(
              id = 'hide-labels',
              options = list(
                list('label' = 'Row', 'value' = 'row'),
                list('label' = 'Column', 'value' = 'col')
              ),
              value = list('row'),
              multi = TRUE
            )
          )),
          
          htmlHr(),
          
          htmlDiv(className = 'app-controls-block', children = list(
            htmlDiv(
              'Color threshold',
              className = 'app-controls-name'
            ),
            
            htmlDiv(
              className = 'app-controls-desc',
              children = 'Change the threshold level that is used to determine
              seperate clusters.'
            )
          )),
          
          htmlBr(),
          
          htmlDiv(
            id = 'threshold-wrapper',
            children = list(
              'Column: ',
              dccSlider(
                id = 'column-threshold',
                className = 'control-slider',
                min = 0,
                max = 20,
                step = 0.5,
                value = 10
              ),
              htmlBr(),
              
              'Row: ',
              dccSlider(
                id = 'row-threshold',
                className = 'control-slider',
                min = 0,
                max = 20,
                step = 0.5,
                value = 10
              )
            )
          ),
          htmlBr(),
          
          htmlHr(),
          
          htmlDiv(
            id = 'add-group-markers',
            children = list(
              htmlDiv(className = 'app-controls-block', children = list(
                htmlDiv(
                  className = 'app-controls-name',
                  children = 'Annotations:'
                ),
                
                htmlButton(
                  id = 'remove-all-group-markers',
                  children = 'Remove all',
                  n_clicks = 0,
                  n_clicks_timestamp = 0
                ),
                
                htmlDiv(className = 'app-controls-desc', children = list(
                  'Annotate your heatmap by labeling clusers; below you can
                  choose a color for the annotation as well as text for the 
                  annotation. Then, click on the row cluster or column cluster
                  that you wish to annotate.'
                ))
              )),
              
              daqColorPicker(
                id = 'clustergram-annot-color',
                value = list('hex' = 'rgb(128, 0, 96)'),
                size = 315
              ),
              
              dccInput(
                id = 'annotation',
                placeholder = 'annotation text',
                type = 'text',
                value = ''
              )
            )
          ),
          
          htmlBr(),
          
          htmlHr(),
          
          htmlDiv(className = 'app-controls-block', children = list(
            htmlDiv(
              'Rows to display:',
              title = 'Select a subset of rows from the uploaded or preloaded dataset
              to compute clustering on.',
              className = 'fullwidth-app-controls-name'
            ),
            
            dccDropdown(
              id = 'selected-rows',
              multi = TRUE,
              value = list()
            )
          )),
          
          htmlDiv(className = 'app-controls-block', children = list(
            htmlDiv(
              'Columns to display:',
              title = 'Select a subset of columns from the uploaded or preloaded dataset
              to compute clustering on.',
              className = 'fullwidth-app-controls-name'
            ),
            
            dccDropdown(
              id = 'selected-columns',
              multi = TRUE,
              value = list()
            )
          ))
          
        )))
      )
    )),
    
    dccStore(
      id = 'data-meta-storage'
    ),
    
    dccStore(
      id = 'fig-options-storage'
    ),
    
    dccStore(
      id = 'computed-traces',
    ),
    
    dccStore(
      id = 'curves-dict'
    ),
    
    dccStore(
      id = 'group-markers'
    )
  ))
))



app$layout(htmlDiv(list(
  options_tabs
)))


# Store file and dataset callback.

# app$callback(
#   output(id = 'data-meta-storage', property = 'data'),
#   params = list(
#     input(id = 'file-upload', property = 'contents'),
#     input(id = 'file-upload', property = 'filename'),
#     input(id = 'clustergram-datasets', property = 'value')
#   ),
#   
#   store_file_meta_data <- function(contents, filename, dataset_name) {
#     desc = ''
#     subsets = list()
#     row_options = list()
#     col_options = list()
#     
#     if (is.null(dataset_name) == FALSE) {
#       dataset = datasets[[dataset_name]]
#     }
#   }
# )

# 
# app$callback(
#   output(id = "row-threshold", property = "value"),
#   params = list(
#     input(id = "clustergram-datasets", property = "value"),
#     input(id = "file-upload", property = "contents")
#   ),
#   
#   update_row_threshold_value <- function(dataset_name, contents) {
#     
#   }
# )



if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server(debug = TRUE, showcase = TRUE)
}









