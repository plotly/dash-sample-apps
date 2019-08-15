library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(data.table)

appName <- Sys.getenv("DASH_APP_NAME")

if (!appName == "") {
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("/app/apps/%s", appName))
}

source("utils/helperFunctions.R")

# Load pre-generated data
portfolioData1 <- readRDS("data/portfolioData1.rds")
portfolioData2 <- readRDS("data/portfolioData2.rds")
portfolioData3 <- readRDS("data/portfolioData3.rds")
# load list of all available symbols
allSymbols <- readRDS("data/allSymbols.rds")

app <- Dash$new()

app$layout(
  htmlDiv(
    className = "container scalable",
    children = list(
      htmlDiv(
        id = "banner-div",
        className = "banner",
        children = list(
          htmlH1(id = "appTitle", children = "Portfolio Optimization"),
          htmlImg(src = "assets/plotly_logo_white.png")
        )
      ),
      htmlBr(),
      htmlDiv(
        id = "app-body",
        #className = "",
        children = list(
          htmlDiv(
            id = "left-div",
            children = list(
              htmlDiv(
                id = "controls-div",
                children = list(
                  dccTabs(
                    children = list(
                      dccTab(
                        label = "Simulate Portfolio",
                        children = list(
                          htmlDiv(
                            className = "control",
                            style = list(
                              paddingBottom = "0px"
                            ),
                            children = list(
                              htmlH5("Generate a portfolio sample:"),
                              htmlP(
                                "Note that generating a portfolio sample may take a long time",
                                style = list(color = "gray")
                              )
                            )
                          ),
                          htmlDiv(
                            className = "control",
                            children = list(
                              htmlH6("Number of Permutations"),
                              dccSlider(
                                id = "nPermutationsSlider",
                                value = 5000,
                                min = 100,
                                max = 50000,
                                step = 1000,
                                marks = list(
                                  "10" = 100,
                                  "10000" = 10000,
                                  "20000" = 20000,
                                  "30000" = 30000,
                                  "40000" = 40000,
                                  "50000" = 50000
                                )
                              )
                            )
                          ),
                          htmlDiv(
                            className = "control",
                            children = list(
                              htmlH6("Symbols"),
                              dccDropdown(
                                id = "symbolsDropdown",
                                options = lapply(
                                  allSymbols,
                                  function(symbol){
                                    list(label = symbol, value = symbol)
                                  }
                                ),
                                multi=TRUE,
                                value = colnames(portfolioData1$rportfolios)
                              )
                            )
                          ),
                          htmlDiv(
                            className = "control",
                            children = list(
                              htmlH6("Random Permutation Method"),
                              dccDropdown(
                                id = "pMethodDropdown",
                                options = list(
                                  list(label = "Sample", value = "sample"),
                                  list(label = "Simplex", value = "simplex"),
                                  list(label = "Grid", value = "grid")
                                ),
                                value = "sample"
                              )
                            )
                          ),
                          htmlDiv(
                            className = "control-submit-store",
                            children = list(
                              htmlButton(
                                id = "resample-button", 
                                children = "resample", 
                                n_clicks = 0
                              ),
                              dccLoading(
                                dccStore(
                                  id = "data-store", data = portfolioData1
                                ),
                              )
                            )
                          )
                        )
                      ),
                      dccTab(
                        label = "Load Portfolio",
                        children = list(
                          htmlDiv(
                            className = "control",
                            children = list(
                              htmlH5("Load a pre-generated portfolio sample:") 
                            ) 
                          ),
                          htmlDiv(
                            className = "control",
                            children = list(
                              dccDropdown(
                                id = "loadDataDropdown",
                                options = list(
                                  list(
                                    label = "MSFT-SBUX-IBM-AAPL-GSPC-AMZN",
                                    value = "portfolioData1"
                                  ),
                                  list(
                                    label = "FB-WORK-AAPL-MSFT-TSLA",
                                    value = "portfolioData2"
                                  ),
                                  list(
                                    label = "SBUX-NXPI-FB-SFIX-JNJ-CNC",
                                    value = "portfolioData3"
                                  )
                                ),
                                value = "portfolioData1"
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              ),
              htmlDiv(
                className = "plot-outer",
                id = "diversification-plot-outer",
                dccLoading(
                  id = "diversification-loading",
                  dccGraph(
                    id = "diversification-plot",
                    figure = generateEmpty()
                  )
                )
              )
            )
          ),
          htmlDiv(
            id = "plots-div",
            children = list(
              htmlDiv(
                className = "plot-outer",
                id = "frontier-plot-outer",
                dccLoading(
                  id = "frontier-loading",
                  dccGraph(
                    id = "frontier-plot",
                    figure = generateEmpty()
                  )
                )
              ),
              htmlDiv(
                className = "plot-outer",
                id = "history-plot-outer",
                dccLoading(
                  id = "history-loading",
                  children = dccGraph(
                    id = "history-plot",
                    figure = generateEmpty()
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
  output("data-store", "data"),
  list(
    input("resample-button", "n_clicks"),
    state("nPermutationsSlider", "value"),
    state("symbolsDropdown", "value"),
    state("pMethodDropdown", "value")
  ),
  function(n_clicks, n_permutations, symbolList, rp_method){
    if (unlist(n_clicks) > 0){
      d <- getSymbolData(unlist(symbolList))
      return(generatePortfolios(d, n_permutations, rp_method))
    }
    return(portfolioData1)
    #return(dashNoUpdate())
  }
)

app$callback(
  output("frontier-plot", "figure"),
  list(
    input("data-store", "data"),
    input("loadDataDropdown", "value")
  ),
  function(d, dropdown){
    ctx <- app$callback_context()
    prop_id <- strsplit(
      as.character(ctx[["triggered"]]["prop_id"]), "\\."
      )[[1]][1]
    if (prop_id == "data-store"){
      d <- list(
        feasible.sd = as.numeric(d$feasible.sd),
        feasible.means = as.numeric(d$feasible.means),
        feasible.sr = as.numeric(d$feasible.sr),
        eff.frontier = as.data.frame(rbindlist(d$eff.frontier, fill = TRUE)),
        eff.frontier.wc = as.data.frame(
          rbindlist(d$eff.frontier.wc, fill = TRUE)
        )
      )
      return(generateFrontierPlot(d))
    } else if (prop_id == "loadDataDropdown"){
      d <- get(dropdown)
      return(generateFrontierPlot(d))
    }
  }
)

app$callback(
  output("diversification-plot", "figure"),
  list(
    input("data-store", "data"),
    input("loadDataDropdown", "value"),
    state("symbolsDropdown", "value")
  ),
  function(d, dropdown, symbolnames){
    ctx <- app$callback_context()
    prop_id <- strsplit(
      as.character(ctx[["triggered"]]["prop_id"]), "\\."
      )[[1]][1]
    if (prop_id == "data-store"){
      fw <- as.data.frame(rbindlist(d$frontier.weights))
      colnames(fw) <- unlist(symbolnames)
      d <- list(
        frontier.weights = fw 
      )
      return(generateDiversificationPlot(d))
    } else if (prop_id == "loadDataDropdown"){
      d <- get(dropdown)
      return(generateDiversificationPlot(d))
    }
  }
)

app$callback(
  output("history-plot", "figure"),
  list(
    input("data-store", "data"),
    input("loadDataDropdown", "value"),
    state("symbolsDropdown", "value")
  ),
  function(d, dropdown, symbolnames){
    ctx <- app$callback_context()
    prop_id <- strsplit(
      as.character(ctx[["triggered"]]["prop_id"]), "\\."
      )[[1]][1]
    if (prop_id == "data-store"){
      pdata <- as.data.frame(rbindlist(d$price.data$price.data))
      pdata <- pdata[2:nrow(pdata),]
      colnames(pdata) <- unlist(symbolnames)
      rownames(pdata) <- unlist(d$dates)
      d <- list(
        price.data = list(
          price.data = pdata 
        )
      )
      generateHistoryPlot(d)
    } else if (prop_id == "loadDataDropdown"){
      d <- get(dropdown)
      return(generateHistoryPlot(d))
    }
  }
)

if (!appName == ""){
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}

