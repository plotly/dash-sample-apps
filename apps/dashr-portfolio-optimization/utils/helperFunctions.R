library(PortfolioAnalytics)
library(quantmod)
library(PerformanceAnalytics)
library(zoo)
library(plotly)
library(parallel)

# Get data
getSymbolData <- function(symbolList){
  getSymbols(symbolList, auto.assign = TRUE)
  # Assign to dataframe
  # Get adjusted prices
  p.data <- zoo(
    get(gsub("^", "", symbolList[1], fixed = TRUE))[, 6]
  )
  for (i in 2:length(symbolList)){
    p.data <- merge.zoo(
      p.data, 
      get(gsub("^", "", symbolList[i], fixed = TRUE))[, 6]
    )
  }
  colnames(p.data) <- symbolList
  # calculate returns
  returns.data <- na.omit(CalculateReturns(p.data))
  colnames(returns.data) <- symbolList
  list(returns.data = returns.data, price.data = p.data)
}

# Generate random portfolios and efficient frontiers
generatePortfolios <- function(
    d, 
    n_permutations, 
    rp_method = c("sample", "simplex", "grid")
  ){
  if (class(d) == "list"){
    returns.data <- d$returns.data
  } else {
    returns.data <- d
  }
  # Save mean return vector and sample covariance matrix
  meanReturns <- colMeans(returns.data)
  covMat <- cov(returns.data)
  # Create portfolio specification
  # Start with the names of the assets
  port <- portfolio.spec(assets = colnames(returns.data))
  # Add constraints to portfolio
  # Box
  port <- add.constraint(port, type = "box", min = 0.05, max = 0.8)
  # Leverage
  port <- add.constraint(portfolio = port, type = "full_investment")
  # Generate random portfolios
  rportfolios <- random_portfolios(
    port, permutations = n_permutations, rp_method = rp_method
  )

  # Get minimum variance portfolio
  minvar.port <- add.objective(port, type = "risk", name = "var")
  # Optimize
  minvar.opt <- optimize.portfolio(
    returns.data, minvar.port, optimize_method = "random", rp = rportfolios
  )
  # Generate maximum return portfolio
  maxret.port <- add.objective(port, type = "return", name = "mean")
  # Optimize
  maxret.opt <- optimize.portfolio(
    returns.data, maxret.port, optimize_method = "random", rp = rportfolios
  )
  # Generate vector of returns
  minret <- 0.06/100
  maxret <- maxret.opt$weights %*% meanReturns
  vec <- seq(minret, as.numeric(maxret), length.out = 100)

  # Generate efficient frontiers
  eff.frontier <- data.frame(Risk = rep(NA, length(vec)),
                             Return = rep(NA, length(vec)), 
                             SharpeRatio = rep(NA, length(vec)))

  frontier.weights <- mat.or.vec(nr = length(vec), nc = ncol(returns.data))
  colnames(frontier.weights) <- colnames(returns.data)

  for(i in 1:length(vec)){
    eff.port <- add.constraint(
      port, type = "return", name = "mean", return_target = vec[i]
    )
    eff.port <- add.objective(eff.port, type = "risk", name = "var")
    eff.port <- optimize.portfolio(
      returns.data, eff.port, optimize_method = "ROI"
    )
    eff.frontier$Risk[i] <- sqrt(
      t(eff.port$weights) %*% covMat %*% eff.port$weights
    )
    eff.frontier$Return[i] <- eff.port$weights %*% meanReturns
    eff.frontier$Sharperatio[i] <- eff.port$Return[i] / eff.port$Risk[i]
    frontier.weights[i,] = eff.port$weights
  }

  # repeat w/ weight concentration
  eff.frontier.wc <- data.frame(Risk = rep(NA, length(vec)),
                             Return = rep(NA, length(vec)), 
                             SharpeRatio = rep(NA, length(vec)))

  frontier.weights.wc <- mat.or.vec(nr = length(vec), nc = ncol(returns.data))
  colnames(frontier.weights.wc) <- colnames(returns.data)

  for(i in 1:length(vec)){
    eff.port <- add.constraint(
      port, type = "return", name = "mean", return_target = vec[i]
    )
    eff.port <- add.objective(eff.port, type = "risk", name = "var")
    eff.port <- add.objective(
      eff.port, type = "weight_concentration", 
      name = "HHI", conc_aversion = 0.001
    )
    eff.port <- optimize.portfolio(
      returns.data, eff.port, optimize_method = "ROI"
    )
    eff.frontier.wc$Risk[i] <- sqrt(
      t(eff.port$weights) %*% covMat %*% eff.port$weights
    )
    eff.frontier.wc$Return[i] <- eff.port$weights %*% meanReturns
    eff.frontier.wc$Sharperatio[i] <- eff.port$Return[i] / eff.port$Risk[i]
    frontier.weights.wc[i,] = eff.port$weights
  }

  feasible.sd <- apply(rportfolios, 1, function(x){
    return(sqrt(matrix(x, nrow = 1) %*% covMat %*% matrix(x, ncol = 1)))
  })

  feasible.means <- apply(rportfolios, 1, function(x){
    return(x %*% meanReturns)
  })

  feasible.sr <- feasible.means / feasible.sd
  dates <- rownames(as.data.frame(returns.data))

  list(
    feasible.means = feasible.means,
    feasible.sd = feasible.sd,
    feasible.sr = feasible.sr,
    eff.frontier = eff.frontier,
    eff.frontier.wc = eff.frontier.wc,
    rportfolios = rportfolios,
    frontier.weights = frontier.weights,
    frontier.weights.wc = frontier.weights.wc,
    price.data = d,
    dates = dates 
  )
}


# Generate plots:
generateFrontierPlot <- function(portfolioData){
  feasible.sd <- as.numeric(portfolioData$feasible.sd)
  feasible.means <- as.numeric(portfolioData$feasible.means)
  feasible.sr <- as.numeric(portfolioData$feasible.sr)
  eff.frontier <- as.data.frame(portfolioData$eff.frontier)
  eff.frontier.wc <- as.data.frame(portfolioData$eff.frontier.wc)
  n_simulations <- length(portfolioData$feasible.means)
  p <- plot_ly(
    x = feasible.sd, y = feasible.means, color = feasible.sr, 
    mode = "markers", type = "scattergl", showlegend = FALSE,
    hoverinfo = "none",
    marker = list(size = 3, opacity = 0.5)
  ) %>% 
    colorbar(title = "Sharpe Ratio", len = 1) %>%
    add_trace(inherit = FALSE, data = eff.frontier, 
              x = eff.frontier$Risk, y = eff.frontier$Return, 
              mode = "markers", 
              type = "scattergl", 
              name = "Efficient Frontier (no weight concentration)",
              marker = list(color = "#F7C873", size = 5)) %>% 
    add_trace(inherit = FALSE, data = eff.frontier.wc,
              x = eff.frontier.wc$Risk, y = eff.frontier.wc$Return, 
              mode = "markers", 
              type = "scattergl",
              name = "Efficient Frontier (with weight concentration)",
              marker = list(color = "#ff471a", size = 5)) %>% 
    layout(
      title = sprintf(
        "Efficient Frontier based on <br> %s Randomly Generated Portfolios",
        n_simulations
      ),
      yaxis = list(
        title = "Mean Returns", tickformat = ".2%", showgrid = FALSE
      ),
      xaxis = list(
        title = "Standard Deviation (Risk)", 
        tickformat = ".2%", showgrid = FALSE
      ),
      plot_bgcolor = "#F8F8F8",
      paper_bgcolor = "#F8F8F8",
      hovermode = "colsest",
      legend = list(
        x = 0,
        y = 1
      )
    )
  p
}

generateDiversificationPlot <- function(portfolioData){
  frontier.weights <- portfolioData$frontier.weights
  colnamefactor <- as.factor(gsub("\\..*$", "", colnames(frontier.weights)))
  q <- plot_ly()
  for (column in dim(frontier.weights)[2]:1){
    q <- q %>%
      add_trace(
        y = frontier.weights[, column],
        type = "bar",
        width = 1,
        color = colnamefactor[column],
        colors = "Paired",
        name = gsub("\\..*$", "", colnames(frontier.weights)[column])
      )
  }
  q <- q %>%
    layout(title = "Portfolio weights across frontier", barmode = "stack",
           plot_bgcolor = "#F8F8F8",
           paper_bgcolor = "#F8F8F8",
           hovermode = "closest",
           xaxis = list(title = "Index"),
           yaxis = list(title = "Weights(%)", tickformat = ".0%"))
  q
}

generateHistoryPlot <- function(portfolioData){
  d <- as.data.frame(portfolioData$price.data$price.data)
  colnamefactor <- as.factor(gsub("\\..*$", "", colnames(d)))
  p <- plot_ly()
  for (column in 1:dim(d)[2]){
    p <- p %>%
      add_trace(
        x = rownames(d),
        #x = as.Date(portfolioData$dates),
        y = d[,column],
        type = "scatter",
        mode = "lines",
        color = colnamefactor[column],
        colors = "Paired",
        name = gsub("\\..*$", "", colnames(d)[column])
      )
  }
  p <- p %>%
    layout(
      title = "Historical Performance",
      xaxis = list(
        showgrid = FALSE,
        type = "date",
        format = "%Y"
        ),
      yaxis = list(
        showgrid = FALSE, 
        title = "Adjusted Price", 
        tickformat = "$"
      ),
      plot_bgcolor = "#F8F8F8",
      paper_bgcolor = "#F8F8F8"
    )
}

getAllSymbols <- function(){
  c(stockSymbols()[,1], list("^GSPC"))
}

# Generate Default data
#allSymbols <- getAllSymbols()
#saveRDS(allSymbols, file = "allSymbols.rds")

# 500000 permutations
#symbolList <- c("MSFT", "SBUX", "IBM", "AAPL", "^GSPC", "AMZN")
#x <- getSymbolData(symbolList)
#portfolioData1 <- generatePortfolios(
  #x, n_permutations = 500000, rp_method = "sample"
#)

#symbolList <- c("FB", "WORK", "AAPL", "MSFT", "TSLA")
#x <- getSymbolData(symbolList)
#portfolioData2 <- generatePortfolios(
  #x, n_permutations = 500000, rp_method = "sample"
#)

#symbolList <- c("SBUX", "NXPI", "FB", "SFIX", "JNJ", "CNC")
#x <- getSymbolData(symbolList)
#portfolioData3 <- generatePortfolios(
  #x, n_permutations = 500000, rp_method = "sample"
#)

#saveRDS(portfolioData1, "data/portfolioData1.rds")
#saveRDS(portfolioData2, "data/portfolioData2.rds")
#saveRDS(portfolioData3, "data/portfolioData3.rds")
