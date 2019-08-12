library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(data.table)

source("helperFunctions.R")

# Load default data
load("data/portfolioData.RData")
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
          htmlH1(id = "appTitle", "Portfolio Optimization"),
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
                list(
                  htmlDiv(
                    className = "control",
                    children = list(
                      htmlH6("Number of Permutations"),
                      dccSlider(
                        id = "nPermutationsSlider",
                        value = 500000,
                        min = 0,
                        max = 500000,
                        step = 1000,
                        marks = list(
                          "1" = 1,
                          "100000" = 100000,
                          "200000" = 200000,
                          "300000" = 300000,
                          "400000" = 400000,
                          "500000" = 500000
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
                  htmlButton(id = "resample-button", children = "resample"),
                  dccLoading(
                    dccStore(id = "data-store", data = portfolioData),
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

#app$callback(
  #output("frontier-plot", "figure"),
  #list(
    #input("resample-button", "n_clicks"),
    #state("nPermutationsSlider", "value"),
    #state("symbolsDropdown", "value"),
    #state("pMethodDropdown", "value")
  #),
  #function(n_clicks, n_permutations, symbolList, rp_method){
    #if (unlist(n_clicks) > 0){
      #returns.data <- getSymbolData(unlist(symbolList))
      #portfolioData <- generatePortfolios(returns.data, n_permutations, rp_method)
      #generateFrontierPlot(portfolioData)
    #}
  #}
#)

app$callback(
  output("data-store", "data"),
  list(
    input("resample-button", "n_clicks"),
    state("nPermutationsSlider", "value"),
    state("symbolsDropdown", "value"),
    state("pMethodDropdown", "value")
    #input("nPermutationsSlider", "value"),
    #input("symbolsDropdown", "value"),
    #input("pMethodDropdown", "value")
  ),
  function(n_clicks, n_permutations, symbolList, rp_method){
    if (unlist(n_clicks) > 0){
      d <- getSymbolData(unlist(symbolList))
      generatePortfolios(d, n_permutations, rp_method)      
    }
  }
)

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
    colnames(fw) <- symbolnames
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
    colnames(pdata) <- symbolnames
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

