appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

library(plotly)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashDaq)
library(redux)
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(rrq)

redux::hiredis()$FLUSHALL()

# Load and subset sample data into usable format. Subsetting necessary or the server may crash.

# sample_data_frame <- data.frame(fread(file = "data/sample_RNA.csv",
#                                       header = TRUE), row.names = 1)
# 
# 
# sample_data <- as.sparse(sample_data_frame)
# sample_data <- CollapseSpeciesExpressionMatrix(sample_data)
# 
# cbmc.adt <- as.sparse(read.csv(file = "data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",",
#                                header = TRUE, row.names = 1)[, 1:1000])
# 
# cbmc.adt <- cbmc.adt[setdiff(rownames(cbmc.adt), c("CCR5", "CCR7", "CD10")), ]

# Code to initialize and generate queue.

my_queue <- rrq::rrq_controller("jobs")

# Layout and Design

header <- htmlDiv(
  id = "app-page-header",
  style = list(
    width = "100%",
    background = "#003c54" ,
    color = "white"
  ),
  children = list(
    htmlA(
      id = "dash-logo",
      children = list(
        htmlImg(src='assets/plotly-dash-logo.png', height = '36', width = '180',
                style = list('top' = '10', 'margin' = '10px'))
      ),
      href = "/Portal"
    ),
    htmlH2("Asynchronous Machine Learning"),
    htmlA(
      id = "gh-link",
      children = list("View on GitHub"),
      href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-uber-rasterizer",
      style = list(color = "white", border = "solid 1px white")
    ),
    htmlImg(
      src = "assets/GitHub-Mark-Light-64px.png"
    )
  )
)


app <- Dash$new()

app$layout(htmlDiv(list(
  header,
  htmlDiv(list(
    htmlDiv(list(
      dccTabs(id = "tabs", value = "tab-1", children = list(
        dccTab(label = "Description", value = "tab-1", children=list(
          dccMarkdown("
                  The Asynchronous Machine Learning Dash app uses a simple K-Nearest Neighbors machine
                  learning algorithm for cell-clustering and visualization of high-throughput single-cell R
                  NA-seq (scRNA-seq) data. Exploring a dataset of 8,617 cord blood mononuclear cells 
                  sequenced with [CITE-seq](https://cite-seq.com/) we will jointly analyze single-cell gene and 
                  protein expression using the [Seurat](https://satijalab.org/seurat/v3.0/multimodal_vignette.html)
                  package adapted for multi-modal data.
                  
                  
                  This app uses an asynchronous task queue with a Redis database backend to deliver concurrent 
                  execution and prioritization of 'tasks'. Machine learning and other computationally intensive
                  processes can be streamlined and optimized for efficiency with multiple simultaneous workers.
                  
                  
                  To get started, queue up one or multiple tasks, then spawn workers to process them.
                  ", style = list("margin" = "10px"))
        ), selected_style = list("borderTop" = "3px solid #0077A7", "borderBottom" = "3px solid #0077A7")),
        dccTab(label = "Controls", value = "tab-2", children=list(
          htmlDiv(list(
            dccMarkdown("
                  Queue up the tasks below by selecting the corresponding buttons. The tasks in this application are KNN clustering 
                  using transcriptome or proteome data indexed with CITE-seq. You must have the Redis server configured with basic settings
                  and running before proceeding.
                  
                  
                  You can test that your Redis db is responding by typing `redis-cli PING` in your command line. 
                  ", style = list("margin" = "10px")),
            htmlBr(),
            htmlButton(
              id = 'first-button',
              n_clicks = 0,
              children = 'RNA  Clustering'
            ),
            htmlButton(
              id = 'second-button',
              n_clicks = 0,
              children = 'Protein Clustering'
            ),
            htmlButton(
              id = 'worker-spawn',
              n_clicks = 0,
              children = 'Spawn a Worker'),
            htmlButton(
              id = 'generate-results',
              n_clicks = 0,
              children = "Generate Results"
            ),
            htmlButton(
              id = 'clear-queue',
              n_clicks = 0,
              children = "Clear Redis Keys",
              style = list('color' = '#8b0000')
            )
          ), className = 'buttons-container')
        ), selected_style = list("borderTop" = "3px solid #0077A7", "borderBottom" = "3px solid #0077A7"))
      ))
    ), className = 'item-a'),
    htmlDiv(list(
      htmlH2("Results", style = list("font-weight" = "200",
                                     "letter-spacing" = "1px",
                                     "text-decoration" = "underline",
                                     "text-decoration-color" = "#5096b3",
                                     "margin-top" = "0px",
                                     "margin-top" = "5px"
      )),
      htmlDiv(id = 'task-1-status'),
      daqIndicator(
        id="task-1-indicator",
        color="#B53737"
      ),
      htmlDiv(id = 'task-1-results',
              style = list(
                "height" = "40%",
                "width" = "70%",
                "align-self" = "center"
              )),
      
      htmlDiv(id = 'task-2-status'),
      daqIndicator(
        id="task-2-indicator",
        color="#B53737"
      ),
      htmlDiv(id = 'task-2-results',
              style = list(
                "height" = "40%",
                "width" = "70%",
                "align-self" = "center"
              ))
    ), className = 'item-b'),
    htmlDiv(list(
      htmlH2("Worker Logs", style = list("font-weight" = "200",
                                         "letter-spacing" = "1px",
                                         "text-decoration" = "underline",
                                         "text-decoration-color" = "#5096b3",
                                         "padding-left" = "2px",
                                         "margin-top" = "0px",
                                         "margin-top" = "5px"
      )),
      dccRadioItems(
        id = 'workers-list',
        options = list(
          list(label = "No Workers", value = "No Workers")
        ),
        value = "No Workers"
      ),
      htmlDiv(id = 'worker-output')
    ), className = 'item-c'),
    dccInterval(
      id = 'update-interval',
      interval = 3*1000,
      n_intervals = 0
    ),
    dccStore(id = 'task-1-store'),
    dccStore(id = 'task-2-store'),
    dccStore(id = 'collection'),
    dccStore(id = 'collection-2')
  ), className = 'container')
)))


# Callback to initialize Tasks A and B

app$callback(output = list(
  output(id = 'task-1-status', property = 'children'),
  output(id = 'task-1-store', property = 'data')
),
params = list(
  input(id = 'first-button', property = 'n_clicks')
),
task_1 <- function(n_clicks) {
  if (n_clicks > 0) {
    task <- my_queue$enqueue({
      library(plotly)
      library(dash)
      library(dashCoreComponents)
      library(dashHtmlComponents)
      library(caret)
      library(redux)
      library(Seurat)
      library(ggplot2)
      library(dplyr)
      library(data.table)
      library(readr)
      library(rrq)
      
      setwd("/home/hammadtheone/Documents/redis_test")
      
      sample_data_frame <- data.frame(fread(file = "data/sample_RNA.csv",
                                            header = TRUE), row.names = 1)
      
      
      sample_data <- as.sparse(sample_data_frame)
      sample_data <- CollapseSpeciesExpressionMatrix(sample_data)
      
      cbmc <- CreateSeuratObject(counts = sample_data)
      
      # normalize data for each cell by total expression
      cbmc <- NormalizeData(cbmc)
      
      # identify highly variable features
      cbmc <- FindVariableFeatures(cbmc)
      
      # scale data so that the mean is 0 and variance is 1 for each gene
      cbmc <- ScaleData(cbmc)
      
      # PCA 
      cbmc <- RunPCA(cbmc, verbose = FALSE)
      
      # Construct a KNN graph
      cbmc <- FindNeighbors(cbmc, dims = 1:25)
      
      # group cells together
      cbmc <- FindClusters(cbmc, resolution = 0.8)
      
      cbmc <- RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")
      
      # identify markers for each cluster
      cbmc.rna.markers <- suppressMessages(
        FindAllMarkers(cbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE))
      
      markers_table <- cbmc.rna.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
      
      # rename clusters with cell type
      new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono", "Mouse", "B", 
                           "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", 
                           "Mouse", "DC", "pDCs")
      
      names(new.cluster.ids) <- levels(cbmc)
      cbmc <- suppressMessages(RenameIdents(cbmc, new.cluster.ids))
      
      final_plot <- ggplotly(DimPlot(cbmc))
    })
    return(list(my_queue$task_status(task),
                names(my_queue$task_status(task))
    )
    )
  }
  return(list("Start a task to view the result.",
              NULL))
}
)

app$callback(output = list(
  output(id = 'task-2-status', property = 'children'),
  output(id = 'task-2-store', property = 'data')
),
params = list(
  input(id = 'second-button', property = 'n_clicks')
),
task_2 <- function(n_clicks) {
  if (n_clicks > 0) {
    taskB <- my_queue$enqueue({
      library(plotly)
      library(dash)
      library(dashCoreComponents)
      library(dashHtmlComponents)
      library(caret)
      library(redux)
      library(Seurat)
      library(ggplot2)
      library(dplyr)
      library(data.table)
      library(readr)
      library(rrq)
      
      setwd("/home/hammadtheone/Documents/redis_test")
      
      sample_data_frame <- data.frame(fread(file = "data/sample_RNA.csv",
                                            header = TRUE), row.names = 1)
      
      sample_data <- as.sparse(sample_data_frame)
      sample_data <- CollapseSpeciesExpressionMatrix(sample_data)
      
      cbmc <- CreateSeuratObject(counts = sample_data)
      
      
      cbmc.adt <- as.sparse(read.csv(file = "data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",",
                                     header = TRUE, row.names = 1)[, 1:1000])
      
      cbmc.adt <- cbmc.adt[setdiff(rownames(cbmc.adt), c("CCR5", "CCR7", "CD10")), ]
      
      # normalize data for each cell by total expression
      cbmc <- NormalizeData(cbmc)
      
      # identify highly variable features
      cbmc <- FindVariableFeatures(cbmc)
      
      # scale data so that the mean is 0 and variance is 1 for each gene
      cbmc <- ScaleData(cbmc)
      
      # PCA 
      cbmc <- RunPCA(cbmc, verbose = FALSE)
      
      # Construct a KNN graph
      cbmc <- FindNeighbors(cbmc, dims = 1:25)
      
      # group cells together
      cbmc <- FindClusters(cbmc, resolution = 0.8)
      
      cbmc <- RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")
      
      # identify markers for each cluster
      cbmc.rna.markers <- suppressMessages(
        FindAllMarkers(cbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE))
      
      cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)
      
      # normalization using the CLR method
      cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
      cbmc <- ScaleData(cbmc, assay = "ADT")
      
      DefaultAssay(cbmc) <- "ADT"
      
      # calculate euclidean distance matrix
      adt.data <- GetAssayData(cbmc, slot = "data")
      adt.dist <- dist(t(adt.data))
      
      # RNA clustering ids
      cbmc[["rnaClusterID"]] <- Idents(cbmc)
      
      # tSNE
      cbmc[["tsne_adt"]] <- RunTSNE(adt.dist, assay="ADT", reduction.key="adtTSNE_")
      
      # clustering
      cbmc[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
      cbmc <- FindClusters(cbmc, resolution=0.2, graph.name="adt_snn")
      
      
      # label clusters
      new.cluster.ids <- c("CD4 T", "CD14+ Mono", "NK", "B", "CD8 T", "NK", "CD34+", "T/Mono doublets", 
                           "CD16+ Mono", "pDCs", "B")
      
      names(new.cluster.ids) <- levels(cbmc)
      cbmc <- RenameIdents(cbmc, new.cluster.ids)
      
      tsne_adtClusters <- ggplotly(DimPlot(cbmc, reduction="tsne_adt"))
    })
    return(list(my_queue$task_status(taskB),
                names(my_queue$task_status(taskB))
    )
    )
  }
  return(list("Start Task B to view the result.",
              NULL))
}
)

# Callbacks to update statuses of Tasks A and B

app$callback(
  output(id = 'task-1-status', property = 'children'),
  params = list(
    input(id = 'first-button', property = 'n_clicks'),
    input(id = 'task-1-store', property = 'data'),
    input(id = 'update-interval', property = 'n_intervals')
  ),
  update_status <- function(n_clicks, data, n_intervals) {
    if (n_clicks > 0) {
      return (my_queue$task_status(data[[1]]))
    }
    else {
      return("Start Task A to view the result.")
    }
  }
)

app$callback(
  output(id = 'task-1-indicator', property = 'color'),
  params = list(
    input(id = 'task-1-store', property = 'data'),
    input(id = 'update-interval', property = 'n_intervals')
  ),
  update_status <- function(data, n_intervals) {
    if (length(my_queue$task_status(data[[1]])) > 0) {
      if (my_queue$task_status(data[[1]])[[1]] == "PENDING") {
        return("#1DA1F2")
      }
      else if (my_queue$task_status(data[[1]])[[1]] == "RUNNING") {
        return("#FFCB05")
      }
      else if (my_queue$task_status(data[[1]])[[1]] == "COMPLETE") {
        return("#00cc96")
      }
      else {
        return("#B53737")
      }
    }
    else {
      return("#777777")
    }
  }
)

app$callback(
  output(id = 'task-2-status', property = 'children'),
  params = list(
    input(id = 'second-button', property = 'n_clicks'),
    input(id = 'task-2-store', property = 'data'),
    input(id = 'update-interval', property = 'n_intervals')
  ),
  update_status <- function(n_clicks, data, n_intervals) {
    if (n_clicks > 0) {
      return (my_queue$task_status(data[[1]]))
    }
    else {
      return("Start Task B to view the result.")
    }
  }
)

app$callback(
  output(id = 'task-2-indicator', property = 'color'),
  params = list(
    input(id = 'task-2-store', property = 'data'),
    input(id = 'update-interval', property = 'n_intervals')
  ),
  update_status <- function(data, n_intervals) {
    if (length(my_queue$task_status(data[[1]])) > 0) {
      if (my_queue$task_status(data[[1]])[[1]] == "PENDING") {
        return("#1DA1F2")
      }
      else if (my_queue$task_status(data[[1]])[[1]] == "RUNNING") {
        return("#FFCB05")
      }
      else if (my_queue$task_status(data[[1]])[[1]] == "COMPLETE") {
        return("#00cc96")
      }
      else {
        return("#B53737")
      }
    }
    else {
      return("#777777")
    }
  }
)

# Callback to spawn workers to complete tasks.

app$callback(
  output(id = 'collection', property = "data"),
  params = list(
    input(id = 'worker-spawn', property = "n_clicks")
  ),
  spawn_worker <- function(n_clicks) {
    if (n_clicks > 0) {
      rrq::worker_spawn(my_queue)
    }
    return("Spawned!")
  }
)

# Callback to display results of Tasks

app$callback(
  output(id = 'task-1-results', property = 'children'),
  params = list(
    input(id = 'generate-results', property = 'n_clicks'),
    input(id = 'task-1-store', property = 'data')
  ),
  output_results <- function(n_clicks, data) {
    if (n_clicks > 0) {
      test <- try(my_queue$task_result(data[[1]]))
      
      if ("try-error" %in% class(test)) {
        return("Task is pending or hasn't been started yet.")
      }
      else { 
        my_figure = my_queue$task_result(data[[1]])
        return(dccGraph(
          id = 'rna_graph',
          figure = my_figure,
          style = list("height" = "100%", "width" = "100%")
        ))
      }
    }
  }
)

app$callback(
  output(id = 'task-2-results', property = 'children'),
  params = list(
    input(id = 'generate-results', property = 'n_clicks'),
    input(id = 'task-2-store', property = 'data')
  ),
  output_results <- function(n_clicks, data) {
    if (n_clicks > 0) {
      test <- try(my_queue$task_result(data[[1]]))
      
      if ("try-error" %in% class(test)) {
        return("Task is pending or hasn't been started yet.")
      }
      else { 
        my_figure = my_queue$task_result(data[[1]])
        return(dccGraph(
          id = 'protein_graph',
          figure = my_figure[[1]],
          style = list("height" = "100%", "width" = "100%")
        ))
      }
    }
  }
)


# Test Callbacks

app$callback(
  output(id = 'worker-output', property = 'children'),
  params = list(
    input(id = 'workers-list', property = 'value'),
    input(id = 'update-interval', property = 'n_intervals')
  ),
  output_worker <- function(worker_name, n_intervals) {
    if (length(my_queue$worker_list()) > 0 && !is.null(worker_name)) {
      return(dccMarkdown(sprintf(
        "
        ```
        %s
        ```
        ", as.character(my_queue$worker_process_log(worker_name))
      )))
    }
  }
)

app$callback(
  output(id = 'workers-list', property = 'options'),
  params = list(
    input(id = 'update-interval', property = 'n_intervals')
  ),
  update_workers_list <- function(n_intervals) {
    if (length(my_queue$worker_list()) > 0) {
      list_options <- lapply(my_queue$worker_list(), function(x) {
        list(label = x, value = x)
      })
    }
    else {
      list_options <- list(
        list(label = "No Workers", value = "No Workers")
      )
    }
    return(list_options)
  }
)


# Callback to wipe the Redis keys, erasing the queue, all workers, and tasks.

app$callback(
  output(id = 'collection-2', property = 'data'),
  params = list(
    input(id = 'clear-queue', property = 'n_clicks')
  ),
  clear_queue <- function(n_clicks) {
    if (n_clicks > 0) {
      redux::hiredis()$FLUSHALL()
    }
  }
)



if(appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(showcase = TRUE)
}
