library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashBio)
library(dashDaq)
library(jsonlite)


appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

DATAPATH <- "data/oncoprint_"

dataset1 <- fromJSON(readLines(sprintf("%sdataset1.json", DATAPATH)))
dataset2 <- fromJSON(readLines(sprintf("%sdataset2.json", DATAPATH)))
dataset3 <- fromJSON(readLines(sprintf("%sdataset3.json", DATAPATH)))
cBioPortalData <- fromJSON(readLines(sprintf("%scBioPortalData.json", DATAPATH)))

DATASETS <- list(
  dataset1 = dataset1,
  dataset2 = dataset2,
  dataset3 = dataset3,
  cBioPortalData = cBioPortalData
)

TRACKS_COLORS_OPT <- list(
  "#aaaaaa",
  "#440154",
  "#472d7b",
  "#3b528b",
  "#2c728e",
  "#21918c",
  "#28ae80",
  "#5ec962",
  "#addc30",
  "#fde725"
)

COLORSCALE_MUTATIONS_OPT <- list(
  "",
  "MISSENSE",
  "INFRAME",
  "FUSION",
  "AMP",
  "GAIN",
  "HETLOSS",
  "HMODEL",
  "UP",
  "DOWN"
)

COLORSCALE_COLORS_OPT <- list(
  "",
  "#440154",
  "#472d7b",
  "#3b528b",
  "#2c728e",
  "#21918c",
  "#28ae80",
  "#5ec962",
  "#addc30",
  "#fde725"
)

TRIGGER_KEY = "trigger"
PADDING_KEY = "padding"
COLORSCALE_KEY = "colorscale"
COLORSCALE_MUT_KEY = "colorscale-mut"
COLORSCALE_COL_KEY = "colorscale-col"

text_style <- list(
  color = "#506784",
  fontFamily = "Open Sans"
)

header_colors <- list(
  bg_color = "#0F5BA7",
  font_color = "white"
)

app <- Dash$new()

app$layout(
  htmlDiv(
    children = list(
      htmlDiv(
        id = "app-page-header",
        style = list(
          width = "100%",
          background = header_colors[["bg_color"]],
          color = header_colors[["font_color"]]
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
          htmlH2("Onco Print"),
          htmlA(
            id = "gh-link",
            children = list("View on GitHub"),
            href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-oncoprint",
            style = list(color = "white", border = "solid 1px white")
          ),
          htmlImg(
            src = "assets/GitHub-Mark-Light-64px.png"
          )
        )
      ),
      htmlBr(),
      htmlDiv(
        id = "oncoprint-body",
        className = "app-body",
        children = list(
          dccLoading(
            className = "dashbio-loading",
            children = dashbioOncoPrint(
              id = "oncoprint-chart",
              height = 550,
              data = list()
            )
          ),
          htmlDiv(
            id = "oncoprint-control-tabs",
            className = "control-tabs",
            children = list(
              dccTabs(
                id = "oncoprint-tabs",
                value = "what-is",
                children = list(
                  dccTab(
                    label = "About",
                    value = "what-is",
                    children = list(
                      htmlDiv(
                        className = "control-tab",
                        children = list(
                          htmlH4(className = "what-is", children = "What is OncoPrint?"),
                          htmlP(
                            paste(
                              "The OncoPrint component is used to view multiple genomic",
                              "alteration events through an interactive and zoomable",
															"heatmap. It is a React/Dash port of the popular",
															"oncoPrint() function from the Bioconductor R",
															"package. Under the hood, the rendering is done with",
															"D3 via Plotly.js. Plotly's interactivity allows",
															"you to bind clicks and hovers to genetic events,",
															"letting you create complex bioinformatics apps",
															"or workflows that leverage crossfiltering."
                            )
                          ),
													htmlP(
														paste(
															"Read more about the component here:",
															"https://github.com/plotly/react-oncoprint"
														)
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
                              children = "Select dataset"
                            ),
                            dccDropdown(
                              id = "oncoprint-dropdown",
                              className = "oncoprint-select",
                              options = lapply(
                                names(DATASETS),
                                function(ds){
                                  list(
                                    label = sprintf("%s.json", ds),
                                    value = ds
                                  )
                                }
                              ),
                              value = "cBioPortalData"
                            )
                          )
                        ),
                        htmlHr(className = "oncoprint-separator"),
                        htmlDiv(
                          list(
                            htmlH4("Event metadata"),
                            htmlDiv(
                              id = "oncoprint-events"
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
                        htmlH4("Layout"),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Overview"
                            ),
                            daqToggleSwitch(
                              id = "oncoprint-show-overview",
                              label = list("hide", "show"),
                              color = "#009DFF",
                              size = 35,
                              value = TRUE
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Legend"
                            ),
                            daqToggleSwitch(
                              id = "oncoprint-show-legend",
                              label = list("hide", "show"),
                              color = "#009DFF",
                              size = 35,
                              value = TRUE
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Padding",
                              style = list(width = "35%")
                            ),
                            dccSlider(
                              className = "control-slider",
                              id = "oncoprint-padding-input",
                              value = 0.05,
                              min = 0,
                              max = 0.1,
                              step = 0.01,
                              marks = as.list(
                                setNames(
                                  as.character(seq(0, 0.1, by = 0.02)),
                                  as.character(seq(0, 0.1, by = 0.02))
                                )
                              )
                            )
                          ) 
                        ),
                        htmlDiv(
                          className = "app-controls-desc",
                          children = "Adjust padding (as percentage) between two tracks."
                        ),
                        htmlHr(className = "oncoprint-separator"),
                        htmlDiv(
                          list(
                            htmlH4("Colors"),
                            htmlP(
                              "Change default background color for the tracks."
                            ),
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "fullwidth-app-controls-name",
                                  children = "Track color"
                                ),
                                daqColorPicker(
                                  id = "oncoprint-tracks-color",
                                  value = list(hex = "#AAAAAA")
                                )
                              )
                            ),
                            htmlHr(className = "oncoprint-separator"),
                            htmlH6("Mutation colors"),
                            htmlP(
                              "Select a mutation type and a color to customize its look."
                            ),
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Mutation type"
                                ),
                                dccDropdown(
                                  id = "oncoprint-colorscale-mutation-dropdown",
                                  options = lapply(
                                    COLORSCALE_MUTATIONS_OPT,
                                    function(mut_type){
                                      list(
                                        label = mut_type,
                                        value = mut_type
                                      )
                                    }
                                  ),
                                  value = COLORSCALE_MUTATIONS_OPT[[1]]
                                )
                              )
                            ),
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Mutation Color"
                                ),
                                daqColorPicker(
                                  id = "oncoprint-mutation-color",
                                  value = list(hex = COLORSCALE_COLORS_OPT[[1]])
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
      ),
      dccStore(id = "oncoprint-store")
    )
  )
)


app$callback(
  output("oncoprint-chart", "data"),
  list(input("oncoprint-dropdown", "value")),
  function(dropdown){
    DATASETS[[dropdown]]
  }
)

app$callback(
  output("oncoprint-store", "data"),
  list(
    input("oncoprint-padding-input", "value"),
    input("oncoprint-colorscale-mutation-dropdown", "value"),
    input("oncoprint-mutation-color", "value"),
    state("oncoprint-store", "data")
  ),
  function(padding_val, mut_type, mut_col, stored_data){
    if (is.null(unlist(stored_data))){
      stored_data <- list(
        padding = "",
        colorscale = list(),
        trigger = ""
      )
    }
    if ((is.null(unlist(mut_col))) | (!"hex" %in% names(mut_col))){
      mut_col <- list(hex = stored_data[[COLORSCALE_KEY]][[mut_type]])
    }
    mut_col = mut_col[["hex"]]
    if (padding_val != stored_data[[PADDING_KEY]]){
      stored_data[[PADDING_KEY]] <- padding_val
      stored_data[[TRIGGER_KEY]] <- PADDING_KEY
    }
    if (!mut_type %in% names(stored_data[[COLORSCALE_KEY]])){
      stored_data[[COLORSCALE_KEY]][[mut_type]] <- mut_col
    } else {
      if (mut_col != stored_data[[COLORSCALE_KEY]][[mut_type]]){
        stored_data[[COLORSCALE_KEY]][[mut_type]] <- mut_col
        stored_data[[TRIGGER_KEY]] <- COLORSCALE_COL_KEY
      } else {
        stored_data[[TRIGGER_KEY]] <- COLORSCALE_MUT_KEY
      }
    }
    stored_data
  }
)

app$callback(
  output("oncoprint-events", "children"),
  list(input("oncoprint-chart", "eventDatum")),
  function(data){
    if ((!is.null(data)) & (nchar(data) > 0)){
      data <- tryCatch(
        {
          fromJSON(data)
        },
        error = function(cond){
          cond
          return("")
        }
      )
      return(
        lapply(
          names(data),
          function(key){
            relabelledKey <- gsub("eventType", "event Type", key, fixed = TRUE)
            htmlDiv(
              sprintf("%s: %s",
                relabelledKey, gsub("<br>", "-", data[[key]], fixed = TRUE)
              )
            )
          }
        )
      )
    }
    return("Hover over or click on a data point on the graph to see it here")
  }
)

app$callback(
  output("oncoprint-chart", "showlegend"),
  list(input("oncoprint-show-legend", "value")),
  function(val){
    val
  }
)

app$callback(
  output("oncoprint-chart", "showoverview"),
  list(input("oncoprint-show-overview", "value")),
  function(val){
    val
  }
)

app$callback(
  output("oncoprint-chart", "backgroundcolor"),
  list(input("oncoprint-tracks-color", "value")),
  function(val){
    if ((!is.null(unlist(val))) & ("hex" %in% names(val))){
      return(val[["hex"]])
    }
    return("#AAAAAA")
  }
)

app$callback(
  output("oncoprint-chart", "padding"),
  list(input("oncoprint-store", "data")),
  function(data){
    data[[PADDING_KEY]]
  }
)

app$callback(
  output("oncoprint-chart", "colorscale"),
  list(input("oncoprint-store", "data")),
  function(data){
    data[[COLORSCALE_KEY]]
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
