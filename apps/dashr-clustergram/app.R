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
}

source("utils/dash_bio_utils.R")

app <- Dash$new()

# Load the datasets.

gg_back_box <- theme(
  panel.background = element_rect(fill = "#191A1C"),
  plot.background = element_rect(fill = "#191A1C"),
  legend.background = element_rect(fill = "#191A1C")
)

# Code to generate header.

header <- htmlDiv(
  id = "app-page-header",
  style = list(
    width = "100%",
    background = "#232323" ,
    color = "#FFFFFF"
  ),
  children = list(
    htmlA(
      id = "dashbio-logo",
      children = list(
        htmlImg(src='assets/plotly-dash-bio-logo.png', height = '36', width = '180',
                style = list('top' = '10', 'margin-left' = '10px'))
      ),
      href = "/Portal"
    ),
    htmlH2("Clustergram"),
    htmlA(
      id = "gh-link",
      children = list("View on GitHub"),
      href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-clustergram",
      style = list(color = "white", border = "solid 1px white")
    ),
    htmlImg(
      src = "assets/GitHub-Mark-Light-64px.png"
    )
  )
)

# Code to generate the options container and adjustable components within.

options_tabs <- htmlDiv(id = 'clustergram-body', className = 'app-body', children = list(
  header,
  htmlBr(),
  htmlBr(),
  htmlBr(),
  dccLoading(className = 'dashbio-loading', type = "circle",  children = htmlDiv(
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
          htmlP('Clustergram is a combination of a heatmap and dendrograms that
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
            value = 'mtcars'
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
                id='upload-data',
                children=htmlDiv(list(
                  'Drag and Drop .soft files or ',
                  htmlA('Select Files')
                )),
                style=list(
                  'width'= '100%',
                  'height'= '60px',
                  'lineHeight'= '60px',
                  'borderWidth'= '1px',
                  'borderStyle'= 'dashed',
                  'borderRadius'= '5px',
                  'textAlign'= 'center'
                ),
                # Allow multiple files to be uploaded
                multiple=FALSE,
                accept = '.SOFT'
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
              href = ("assets/sample_data/clustergram_GDS3627.soft"),
              download = 'clustergram_GDS3627.soft'
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
                list('label' = 'Column', 'value' = 'column')
              ),
              value = 'row',
              multi = FALSE
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
                list('label' = 'Rows', 'value' = "row"),
                list('label' = 'Columns', 'value' = "col"),
                list('label' = 'Both', 'value' = "both"),
                list('label' = 'None', 'value' = "none")
              ),
              value = "row",
              multi = FALSE
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
                daqToggleSwitch(
                  id = 'add-all-group-markers',
                  color = '#191A1C',
                  label = list('Enable', 'Disable'),
                  size = 35,
                  labelPosition = 'bottom',
                  value = FALSE
                ),
                htmlBr(),
                htmlBr(),
                
                htmlDiv(className = 'app-controls-desc', children = list(
                  'Annotate your heatmap by labeling clusers.'
                ))
              ))
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
      id = 'dataset-storage'
    ),
    
    dccStore(
      id = 'data-meta-storage'
    )
  ))
))



app$layout(htmlDiv(list(
  options_tabs,
  htmlDiv(list(
    "                "
  ), className = "blankspace")
)))


# Callback to load data storage based on selected dataset.

app$callback(
  output(id = 'dataset-storage', property = 'data'),
  params = list(
    input(id = 'clustergram-datasets', property = 'value'),
    input(id = "upload-data", property = "contents"),
    input(id = 'upload-data', property = 'filename')
  ),
  
  update_data <- function(selection, uploaded_file, filename) {
    if (selection == "iris") {
      return(
        list("Name" = "Iris Data Set",
             "Description" =  "This is perhaps the best known database to be found
             in the pattern recognition literature. Fisher's paper is a classic in
             the field and is referenced frequently to this day. (See Duda & Hart, 
             for example.) The data set contains 3 classes of 50 instances each, where 
             each class refers to a type of iris plant. One class is linearly 
             separable from the other 2; the latter are NOT linearly separable 
             from each other. [from https://archive.ics.uci.edu/ml/datasets/iris]",
             
             "Source" = "https://archive.ics.uci.edu/ml/machine-learning-databases/iris/",
             
             "filepath" = "assets/sample_data/clustergram_iris.tsv"
        )
      )
    }
    else if (selection == "mtcars") {
      return(
        list("Name" = "Motor Trend Car Road Tests",
             "Description" =  "The data was extracted from the 1974 Motor Trend US
             magazine, and comprises fuel consumption and 10 aspects of automobile design
             and performance for 32 automobiles (1973-74 models).
             [from https://www.rdocumentation.org/packages/datasets/versions/3.5.1/topics/mtcars]",
             
             "Source" = "https://gist.github.com/seankross/a412dfbd88b3db70b74b",
             
             "filepath" = "assets/sample_data/clustergram_mtcars.tsv"
        )
      )
    }
    
    else if (selection == "prostatecancer") {
      return(
        list("Name" = "miR-221 expression effect on prostate cancer cell line",
             "Description" =  "Analysis of PC-3 prostate cancer cells expressing pre-miR-221. 
             miR-221 is frequently downregulated in primary prostate cancer. Results provide 
             insight into the role of miR-221 in the pathogenesis of prostate cancer.",
             
             "Meta" = Meta(getGEO(filename = "assets/sample_data/clustergram_GDS5373.SOFT")),
             
             "filepath" = "assets/sample_data/clustergram_GDS5373.SOFT"
        )
      )
    }
    
    else if (selection == "lungcancer") {
      return(
        list("Name" = "Non-small lung cancer subtypes: adenocarcinoma and squamous cell carcinoma",
             "Description" =  "Comparison of two non-small cell lung cancer histological subtypes: 
             adenocarcinomas (AC) and squamous cell carcinomas (SCC). Results provide insight 
             into the molecular differences between AC and SCC.",
             
             "Meta" = Meta(getGEO(filename = "assets/sample_data/clustergram_GDS3627.SOFT")),
             
             "filepath" = "assets/sample_data/clustergram_GDS3627.SOFT"
        )
      )
    }
  }
  
)



# Callback to change the depicted clustergram.


app$callback(
  output(id = "clustergram-wrapper", property = "children"),
  params = list(
    input(id = "dataset-storage", property = "data"),
    input(id = "hide-labels", property = "value"),
    input(id = "cluster-checklist", property = "value"),
    input(id = "column-threshold", property = "value"),
    input(id = "row-threshold", property = "value"),
    input(id = "selected-rows", property = "value"),
    input(id = "selected-columns", property = "value"),
    input(id = "add-all-group-markers", property = "value")
  ),
  
  update_clustergram <- function(data, labels, cluster, columnslider, rowslider, selected_rows,
                                 selected_columns, annotations) {
    
    if (endsWith(data$filepath, ".tsv")) {
      df <- read.table(file = data$filepath, sep = 
                         "\t", skip = 4, row.names = 1,  header = TRUE)
      
    }
    
    else if (endsWith(data$filepath, ".SOFT")) {
      df <- importSOFT(data$filepath)
      
      df <- df[1:10, 1:5]
      
    }
    
    if (labels == "row") {
      showlabels = c(T,F)
    } 
    
    else if (labels == "col") {
      showlabels = c(F,T)
    }   
    
    else if (labels == "both") {
      showlabels = c(F,F)
    }
    
    else if (labels == "none") {
      showlabels = c(T,T)
    }
    
    if (length(selected_rows) < 2) {
      return("Please select at least two rows to display.")
    }
    else if (length(selected_rows) > 2){
      
      df <- subset(df, rownames(df) %in% selected_rows)
    }
    
    
    if (length(selected_rows) < 2) {
      return("Please select at least two rows to display.")
    }
    else if (length(selected_rows) > 2){
      
      df <- subset(df, colnames(df) %in% selected_columns)
    }
    
    if (annotations == TRUE) {
      heatmap <- dccGraph(figure = heatmaply(df,
                                             heatmap_layers = gg_back_box,
                                             side_color_layers = gg_back_box,
                                             colors = plasma,
                                             showticklabels = showlabels,
                                             scale = cluster,
                                             k_col = 4,
                                             k_row = 4,
                                             column_labels = as.list(colnames(df)),
                                             color_threshold = list(
                                               "row" = rowslider,
                                               "col" = columnslider
                                             ),
                                             seriate = "mean",
                                             row_side_colors = df[, 2:3],
                                             col_side_colors = c(rep(0, length(colnames(df))-2), rep(1,2))
      ))  
    }
    
    else if (annotations == FALSE) {
      heatmap <- dccGraph(figure = heatmaply(df,
                                             heatmap_layers = gg_back_box,
                                             side_color_layers = gg_back_box,
                                             colors = plasma,
                                             showticklabels = showlabels,
                                             scale = cluster,
                                             k_col = 4,
                                             k_row = 4,
                                             column_labels = as.list(colnames(df)),
                                             color_threshold = list(
                                               "row" = rowslider,
                                               "col" = columnslider
                                             )
      ))  
    }
    
    
    return(heatmap)
  }
)


# Callback to change the dataset information and details.

app$callback(
  output(id = "clustergram-info", property = "children"),
  params = list(
    input(id = "dataset-storage", property = "data")
  ),
  
  update_info <- function(value) {
    if (endsWith(value$filepath, ".tsv")) {
      return(htmlDiv(list(
        htmlH2(children = "Dataset Information"),
        htmlBr(),
        sprintf("Name: %s", value$Name),
        htmlBr(),
        htmlBr(),
        sprintf("Description: %s", value$Description),
        htmlBr(),
        htmlBr(),
        sprintf("Source: %s", value$Source)
      )))
    }
    
    else if (endsWith(value$filepath, ".SOFT")) {
      return(htmlDiv(list(
        htmlH2(children = "Dataset Information"),
        htmlBr(),
        sprintf("Name: %s", value$Name),
        htmlBr(),
        htmlBr(),
        sprintf("Description: %s", value$Description),
        paste(names(value$Meta), value$Meta, sep = ": ", collapse = "\n")
      )))
    }
  }
)

# Callback to update options for selected rows.


app$callback(
  output(id = "selected-rows", property = "options"),
  params = list(
    input(id = "dataset-storage", property = "data")
  ),
  
  update_rows_options <- function(data) {
    if (endsWith(data$filepath, ".tsv")) {
      df <- read.table(file = data$filepath, sep = 
                         "\t", skip = 4, row.names = 1,  header = TRUE)
      
    }
    
    else if (endsWith(data$filepath, ".SOFT")) {
      df <- importSOFT(data$filepath)
      
      df <- df[1:10,]
      
    }
    all_options <- rownames(df)
    
    options_list <- list()
    
    for (x in 1:length(all_options)){
      option = list("label" = all_options[x], "value" = all_options[x])
      options_list[[x]] = option
    }
    
    return(options_list)
  }
)


#Callback to update values for selected rows.

app$callback(
  output(id = "selected-rows", property = "value"),
  params = list(
    input(id = "dataset-storage", property = "data")
  ),
  
  update_rows_options <- function(data) {
    if (endsWith(data$filepath, ".tsv")) {
      df <- read.table(file = data$filepath, sep = 
                         "\t", skip = 4, row.names = 1,  header = TRUE)
      
    }
    
    else if (endsWith(data$filepath, ".SOFT")) {
      df <- importSOFT(data$filepath)
      
      df <- df[1:10,]
      
    }
    all_options <- rownames(df)
    
    
    
    return(all_options)
  }
)

# Callback, same as above, but for updating options for selected columns. 

app$callback(
  output(id = "selected-columns", property = "options"),
  params = list(
    input(id = "dataset-storage", property = "data")
  ),
  
  update_rows_options <- function(data) {
    if (endsWith(data$filepath, ".tsv")) {
      df <- read.table(file = data$filepath, sep = 
                         "\t", skip = 4, row.names = 1,  header = TRUE)
      
    }
    
    else if (endsWith(data$filepath, ".SOFT")) {
      df <- importSOFT(data$filepath)
      
      df <- df[1:10,]
      
    }
    all_options <- colnames(df)
    
    options_list <- list()
    
    for (x in 1:length(all_options)){
      option = list("label" = all_options[x], "value" = all_options[x])
      options_list[[x]] = option
    }
    
    return(options_list)
  }
)



#Callback to update values for selected columns.

app$callback(
  output(id = "selected-columns", property = "value"),
  params = list(
    input(id = "dataset-storage", property = "data")
  ),
  
  update_rows_options <- function(data) {
    if (endsWith(data$filepath, ".tsv")) {
      df <- read.table(file = data$filepath, sep = 
                         "\t", skip = 4, row.names = 1,  header = TRUE)
      
    }
    
    else if (endsWith(data$filepath, ".SOFT")) {
      df <- importSOFT(data$filepath)
      
      df <- df[1:10,]
      
    }
    all_options <- colnames(df)
    
    
    
    return(all_options)
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server(debug = TRUE, showcase = TRUE)
}
