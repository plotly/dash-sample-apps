library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashDaq)
library(dashBio)

DATAPATH <- "sample_data/volcano_" 

DATASETS <- list(
  SET1 = list(
    label = "Set1",
    dataframe = NULL,
    datafile = sprintf("%sdata1.csv", DATAPATH),
    datasource = "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/",
    dataprops = list()
  ),
  SET2 = list(
    label = "Set2",
    dataframe = NULL,
    datafile = sprintf("%sdata2.csv", DATAPATH),
    datasource = "https://doi.org/10.1371/journal.pntd.0001039.s001",
    dataprops = list(
      effect_size = "log2_.L3i.L1._signal_ratio",
      p = "p.value",
      #SNP column is required, although the volcanor documentation is unclear 
      snp = "X",
      gene = "PFAM_database_id",
      annotation1 = "annotation"
    )
  )
)

for (dataset in names(DATASETS)){
  DATASETS[[dataset]][["dataframe"]] <- read.csv(
    DATASETS[[dataset]][["datafile"]], comment.char = "#"
  )
}

# Right now the only way is to rename the column to `EFFECTSIZE`
colnames(DATASETS[["SET2"]][["dataframe"]])[[3]] <- "EFFECTSIZE"

description <- function(){
  paste(
    "Interactively identify clinically meaningful markers",
    "in genomic experiments with this volcano plot"
  )
}

header_colors <- function(){
  list(
    bg_color = "#19d3f3",
    font_color = "white",
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
          htmlH2("Volcano Plot"),
          htmlA(
            id = "gh-link",
            children = list("View on GitHub"),
            href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-volcanoplot",
            style = list(color = "white", border = "solid 1px white")
          ),
          htmlImg(
            src = "assets/GitHub-Mark-Light-64px.png"
          )
        )
      ),
      htmlBr(),
      htmlDiv(
        id = "vp-page-content",
        style = list(paddingTop = "50px", minHeight = "calc(100vh - 70px)"),
        className = "app-body",
        children = list(
          dccLoading(
            className = "dashbio-loading",
            children = list(
              htmlDiv(
                id = "vp-graph-div",
                children = dccGraph(
                  id = "vp-graph"
                )
              )
            )
          ),
          htmlDiv(
            id = "vp-control-tabs",
            className = "control-tabs",
            children = list(
              dccTabs(
                id = "vp-tabs",
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
                          children = "What is Volcano Plot?"
                        ),
                        htmlP(
                          paste0(
                            "You can use Volcano Plot to interactively ",
                            "identify clinically meaningful markers in ",
                            "genomic experiments, i.e., markers that are ",
                            "statistically significant and have an effect ",
                            "size greater than some threshold. ",
                            "Specifically, volcano plots depict the negative ",
                            "log-base-10 p-values plotted against their ",
                            "effect size."
                          )
                        ),
                        htmlP(
                          paste0(
                            "In the \"Data\" tab, you can select a dataset ",
                            "to view on the plot. In the \"View\" tab, you ",
                            "can control the color of the highlighted ",
                            "points, as well as the threshold lines that ",
                            "define which values are significant. You can ",
                            "also access metadata from hovering and ",
                            "clicking on the graph."
                          ) 
                        )
                      )
                    )
                  ),
                  dccTab(
                    label = "Data",
                    value = "data",
                    children = htmlDiv(
                      className = "control-tab",
                      children = list(
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Dataset: "
                            ),
                            dccDropdown(
                              id = "vp-dataset-dropdown",
                              options = lapply(
                                names(DATASETS),
                                function(dset){
                                  list(
                                    label = DATASETS[[dset]][["label"]],
                                    value = dset
                                  )
                                }
                              ),
                              value = "SET2"
                            )
                          )
                        )
                      )
                    )
                  ),
                  dccTab(
                    label = "View",
                    value = "view",
                    children = htmlDiv(
                      className = "control-tab",
                      children = list(
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Effect size bounds"
                            ),
                            dccRangeSlider(
                              id = "vp-bound-val",
                              min = -4,
                              max = 4,
                              value = list(-1, 1),
                              step = 0.01,
                              marks = as.list(setNames( -4:4, -4:4))
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Threshold"
                            ),
                            dccSlider(
                              id = "vp-genomic-line-val",
                              value = 4,
                              max = 10,
                              min = 0,
                              step = 0.1,
                              marks = as.list(
                                setNames(
                                  seq(0, 10, by = 2), 
                                  seq(0, 10, by = 2)
                                )
                              )
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            daqColorPicker(
                              id = "vp-color-picker",
                              value = list(hex = "#0000FF"),
                              size = 150
                            ),
                            htmlDiv(
                              id = "vp-num-points-display",
                              children = list(
                                htmlDiv(
                                  title = "Number of points in the upper left",
                                  children = list(
                                    daqLEDDisplay(
                                      className = "vp-input-like",
                                      label = "Upper left points",
                                      id = "vp-upper-left",
                                      size = 25,
                                      color = "#19D3F3"
                                    ),
                                    htmlDiv(
                                      className = "vp-test-util-div",
                                      id = "vp-upper-left-val"
                                    )
                                  )
                                ),
                                htmlBr(),
                                htmlDiv(
                                  className = "vp-vertical-style",
                                  title = "Number of points in the upper right",
                                  children = list(
                                    daqLEDDisplay(
                                      className = "vp-input-like",
                                      label = "Upper right points",
                                      id = "vp-upper-right",
                                      size = 25,
                                      color = "#19D3F3"
                                    ),
                                    htmlDiv(
                                      className = "vp-test-util-div",
                                      id = "upper-right-val"
                                    )
                                  )
                                )
                              )
                            )
                          )
                        ),
                        htmlHr(),
                        htmlDiv(id = "vp-event-data")
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
  output("vp-graph", "figure"),
  list(
    input("vp-bound-val", "value"),
    input("vp-genomic-line-val", "value"),
    input("vp-dataset-dropdown", "value"),
    input("vp-color-picker", "value")
  ),
  function(effect_lims, genomic_line, dataset_id, color){
    l_lim = effect_lims[[1]]
    u_lim = effect_lims[[2]]
    color <- ifelse("hex" %in% names(color), color[["hex"]], "red")
    arguments <- c(
      list(
        dataframe = DATASETS[[dataset_id]][["dataframe"]],
        genomewideline_value = as.numeric(genomic_line),
        effect_size_line = c(l_lim, u_lim),
        highlight_color = color
      ),
      DATASETS[[dataset_id]][["dataprops"]]
    )
    do.call(
      dashbioVolcano, 
      args = arguments
    )
  }
)

app$callback(
  output("vp-event-data", "children"),
  list(
    input("vp-graph", "hoverData"),
    input("vp-graph", "clickData")
  ),
  function(hover, click){
    hover_data_div <- list(
      htmlDiv(className = "app-controls-name", children = "Hover data")
    )
    hover_data <- "Hover over a data point to see it here."
    if (!is.null(unlist(hover))){
      hovered_point <- hover[["points"]][[1]]
      hovered_text <- unlist(strsplit(
        gsub("<br />", "<br>", 
        gsub("<br>+$", "", hovered_point["text"])), 
        "<br>"
      ))
      idx <- length(hovered_text)
      hover_data <- list(
        sprintf("x: %s", hovered_point[["x"]]),
        htmlBr(),
        sprintf("y: %s", hovered_point[["y"]]),
        htmlBr(),
        sprintf(
          "%s, (%s)", 
          hovered_text[[idx - 1]], 
          hovered_text[[idx]]
        )
      )
    }
    hover_data_div <- c(
      hover_data_div,
      list(htmlDiv(className = "vp-event-data-display", children = hover_data))
    )

    click_data_div <- list(
      htmlDiv(className = "app-controls-name", children = "Click data")
    )
    click_data <- "Click on a data point to see it here."
    if (!is.null(unlist(click))){
      clicked_point <- click[["points"]][[1]]
      clicked_text <- unlist(strsplit(
        gsub("<br />", "<br>", 
        gsub("<br>+$", "", clicked_point["text"])), 
        "<br>"
      ))
      idx <- length(clicked_text)
      click_data <- list(
        sprintf("x: %s", clicked_point[["x"]]),
        htmlBr(),
        sprintf("y: %s", clicked_point[["y"]]),
        htmlBr(),
        sprintf(
          "%s, (%s)", 
          clicked_text[[idx - 1]], 
          clicked_text[[idx]]
        )
      )
    }
    click_data_div <- c(
      click_data_div,
      list(htmlDiv(className = "vp-event-data-display", children = click_data))
    )

    htmlDiv(
      list(
        htmlDiv(hover_data_div),
        htmlDiv(click_data_div)
      )
    )
  }
)

app$callback(
  output("vp-upper-right", "value"),
  list(
    input("vp-graph", "figure"),
    input("vp-genomic-line-val", "value"),
    state("vp-bound-val", "value")
  ),
  function(fig, thresh, bounds){
    u_lim <- bounds[[2]]
    number = 0
    if (length(fig[["data"]]) > 1){
      x <- unlist(fig[["data"]][[1]]["x"])
      y <- unlist(fig[["data"]][[1]]["y"])
      number <- sum((x > u_lim) & (y > thresh))
    }
    number
  }
)

app$callback(
  output("vp-upper-left", "value"),
  list(
    input("vp-graph", "figure"),
    input("vp-genomic-line-val", "value"),
    state("vp-bound-val", "value")
  ),
  function(fig, thresh, bounds){
    u_lim <- bounds[[1]]
    number = 0
    if (length(fig[["data"]]) > 1){
      x <- unlist(fig[["data"]][[1]]["x"])
      y <- unlist(fig[["data"]][[1]]["y"])
      number <- sum((x < u_lim) & (y > thresh))
    }
    number
  }
)

app$run_server()

