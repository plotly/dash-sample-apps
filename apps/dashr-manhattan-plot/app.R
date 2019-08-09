library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashDaq)
library(dashBio)

dataset <- read.csv("sample_data/manhattan_data.csv",
                    check.names = FALSE,
                    stringsAsFactors = FALSE)

description <- function(){
  paste(
    "Create an interactive Manhattan plot with multiple",
    "annotation options. This plot can be used to display",
    "the results of genomic studies sorted out by chromosome."
  )
}

header_colors <- function(){
  list(
    bg_color = "#0D76BF",
    font_color = "#fff",
    "light_logo" = TRUE
  )
}

app <- Dash$new()

app$layout(
  htmlDiv(
    children = list(
      htmlDiv(
        id = "app-page-header",
        style = list(
          width = "100%",
          background = header_colors()[["bg_color"]],
          color = header_colors()[["font_color"]]
        ),
        children = list(
          htmlA(
            id = "dashbio-logo",
            children = list(
              htmlImg(src='assets/plotly-dash-bio-logo.png', height = '36', width = '190',
                      style = list('top' = '10', 'margin-left' = '10px'))
            ),
            href = "/Portal"
          ),
          htmlH2("Manhattan Plot"),
          htmlA(
            id = "gh-link",
            children = list("View on GitHub"),
            href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-manhattan-plot",
            style = list(color = "white", border = "solid 1px white")
          ),
          htmlImg(
            src = "assets/GitHub-Mark-Light-64px.png"
          )
        )
      ),
      htmlBr(),
      htmlDiv(
        id = "mhp-page-content",
        style = list(paddingTop = "50px", minHeight = "calc(100vh - 70px)"),
        className = "app-body",
        children = list(
          dccLoading(
            className = "dashbio-loading",
            children = list(
              htmlDiv(
                id = "mhp-graph-div",
                children = dccGraph(
                  id = "mhp-graph"
                )
              )
            )
          ),
          htmlDiv(
            id = "mhp-control-tabs",
            className = "control-tabs",
            children = list(
              dccTabs(
                id = "mhp-tabs",
                value = "what-is",
                children = list(
                  dccTab(
                    label = "About",
                    value = "what-is",
                    children = htmlDiv(
                      className = "control-tab",
                      children = list(
                        htmlH4(
                          className = "what-is", 
                          children = "What is Manhattan Plot?"
                        ),
                        htmlP(
                          paste0(
                            "Manhattan Plot allows you to visualize ",
                            "genome-wide association studies (GWAS) ",
                            "efficiently. Using WebGL under the hood, you ",
                            "can interactively explore overviews of massive ",
                            "datasets comprising hundreds of thousands of ",
                            "points at once, or take a closer look at a small ",
                            "subset of your data."
                          )
                        ),
                        htmlP(
                          paste0(
                            "You can adjust the threshold level and the ",
                            "suggestive line in the \"Graph\" tab."                            
                          ) 
                        )
                      )
                    )
                  ),
                  dccTab(
                    label = "Graph",
                    value = "graph",
                    children = htmlDiv(
                      className = "control-tab",
                      children = list(
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Threshold value (red)"
                            ),
                            dccSlider(
                              id = "mhp-threshold-val",
                              min = 1,
                              max = 4,
                              value = 2,
                              step = 0.01,
                              marks = as.list(
                                setNames(
                                  as.character(1:4), 
                                  as.character(1:4)
                                  )
                                )
                            )
                          )
                        ),
                        htmlBr(),
                        htmlBr(),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Suggestive line (purple)"
                            ),
                            dccSlider(
                              id = "mhp-suggestive-line-val",
                              value = 3,
                              min = 1,
                              max = 4,
                              step = 0.1,
                              marks = as.list(
                                setNames(
                                  as.character(1:4), 
                                  as.character(1:4)
                                  )
                                )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

app$callback(
  output("mhp-graph", "figure"),
  list(
    input("mhp-threshold-val", "value"),
    input("mhp-suggestive-line-val", "value")
  ),
  function(threshold_vals, suggestiveline_value){
    arguments <- list( 
        dataframe = dataset,
        suggestiveline_value = as.numeric(suggestiveline_value),
        genomewideline_value = as.numeric(threshold_vals),
        showlegend = TRUE
      )
    do.call(
      dashbioManhattan, 
      args = arguments
    )
  }
)

app$run_server()
