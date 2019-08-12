library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(data.table)

source("helperFunctions.R")

# Load default data
load("data/MSFT-SBUX-IBM-AAPL-GSPC-AMZN_500000.RData")
# load list of all symbols
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
                                    value = "MSFT-SBUX-IBM-AAPL-GSPC-AMZN"
                                  ),
                                  list(
                                    label = "FB-WORK-AAPL-MSFT-TSLA",
                                    value = "FB-WORK-AAPL-MSFT-TSLA"
                                  )
                                ),
                                value = "MSFT-SBUX-IBM-AAPL-GSPC-AMZN"
                              )
                            )
                          )
                        )
                      ),
                      dccTab(
                        label = "Simulate Portfolio",
                        children = list(
                          htmlDiv(
                            className = "control",
                            children = htmlH5("Generate a portfolio sample:")
                          ),
                          htmlDiv(
                            className = "control",
                            children = list(
                              htmlH6("Number of Permutations"),
                              dccSlider(
                                id = "nPermutationsSlider",
                                value = 5000,
                                min = 10,
                                max = 50000,
                                step = 1000,
                                marks = list(
                                  "10" = 10,
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
                                value = colnames(portfolioData$rportfolios)
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
                              htmlButton(id = "resample-button", children = "resample"),
                              dccLoading(
                                dccStore(id = "data-store", data = portfolioData),
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
                    figure = generateDiversificationPlot(portfolioData)
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
                    figure = generateFrontierPlot(portfolioData)
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
                    figure = generateHistoryPlot(portfolioData)
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
    #if (!is.null(unlist(n_clicks))){
      d <- getSymbolData(unlist(symbolList))
      generatePortfolios(d, n_permutations, rp_method)
      }
    }
)

#app$callback(
  #output("frontier-plot", "figure"),
  #list(
    #input("loadDataDropdown", "value")
  #),
  #function(dataset){
    #print(dataset)
    ##fname <- sprintf("data/%s_500000.RData", unlist(dataset))
    ##print(fname)
    ##load(fname)
    ##generateFrontierPlot(portfolioData)     
  #}
#)

app$callback(
  output("frontier-plot", "figure"),
  list(
    input("data-store", "data")
  ),
  function(d){
    d <- list(
      feasible.sd = as.numeric(d$feasible.sd),
      feasible.means = as.numeric(d$feasible.means),
      feasible.sr = as.numeric(d$feasible.sr),
      eff.frontier = as.data.frame(rbindlist(d$eff.frontier)),
      eff.frontier.wc = as.data.frame(rbindlist(d$eff.frontier.wc))
    )
    generateFrontierPlot(d)
  }
)

app$callback(
  output("diversification-plot", "figure"),
  list(
    input("data-store", "data"),
    state("symbolsDropdown", "value")
  ),
  function(d, symbolnames){
    fw <- as.data.frame(rbindlist(d$frontier.weights))
    colnames(fw) <- unlist(symbolnames)
    d <- list(
      frontier.weights = fw 
    )
    generateDiversificationPlot(d)
  }
)

app$callback(
  output("history-plot", "figure"),
  list(
    input("data-store", "data"),
    state("symbolsDropdown", "value")
  ),
  function(d, symbolnames){
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
  }
)

app$run_server()

