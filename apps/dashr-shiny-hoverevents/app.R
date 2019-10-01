appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("/app/apps/%s", appName))
}

library(plotly)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(reshape2)
library(dplyr)


# Load & Prep Data ------------------------------

# Read data
stockdata <- read.csv("assets/stockdata.csv")

# Create dates
stockdata$Date <- as.Date(stockdata$Date)

# Reshape
ds <- reshape2::melt(stockdata, id = "Date")
ds <- filter(ds, variable != "GSPC")

# Set some colors
plotcolor <- "#F5F1DA"
papercolor <- "#E3DFC8"

# Add SP500
timeseriesFigure <- plot_ly(source = "source")
timeseriesFigure <- timeseriesFigure %>%
  add_lines(
    data = ds,
    x = ~Date,
    y = ~value,
    color = ~variable,
    mode = "lines",
    line = list(width = 3))

timeseriesFigure <- timeseriesFigure %>%
  add_lines(data = stockdata,
            x = ~Date,
            y = ~GSPC,
            mode = "lines",
            yaxis = "y2",
            name = "SP500",
            opacity = 0.3,
            line = list(width = 5)) %>%
  layout(title = "Stock prices for different stocks overlaid with SP500",
         xaxis = list(title = "Dates", gridcolor = "#bfbfbf", domain = c(0, 0.98)),
         yaxis = list(title = "Stock Price", gridcolor = "#bfbfbf"),
         plot_bgcolor = plotcolor,
         paper_bgcolor = papercolor,
         yaxis2 = list(title = "SP500", side = "right", overlaying = "y"))


# Create Layout Variables ------------------------------

plotlyLogo <-
  htmlA(list(htmlImg(id = "banner-image", src = "assets/image.png")),
        className = "logo",
        href = "https://dashr.plot.ly")

windowOptions <- lapply(
  c(10, 20, 30, 60, 90),
  function(x) {
    list(label = x, value = x)
  })

windowLength <- dccDropdown(
  id = "window-length",
  options = windowOptions,
  value = 10,
  style = list("width" = "300px")
)


# App Start ------------------------------

# Initiate application
app <- Dash$new(name = "dashr-shiny-hoverevents")


# Create Layout ------------------------------

app$layout(
  htmlDiv(list(
    plotlyLogo,
    htmlH2("Coupled hover-events in plotly charts using Dash"),
    htmlDiv(list(
      htmlH4("This Shiny app showcases coupled hover-events using Dash's ",
             style = list("display" = "inline")),
      htmlCode(list("hoverData"),
               style = list("display" = "inline")),
      htmlH4(" property.", style = list("display" = "inline"))
    )),
    htmlHr(),
    htmlDiv(list(
      htmlStrong("Select Window Length"),
      windowLength
    )),
    htmlBr(),
    htmlDiv(
      className = "row",
      list(
        htmlDiv(list(
          dccGraph(id = "timeseries", figure = timeseriesFigure)
        ), className = "six columns"),
        htmlDiv(list(
          dccGraph(id = "correlation", style = list("display" = "none"))
        ), className = "six columns")
      )),
    htmlHr(),
    htmlBlockquote("Hover over time series chart to fix a specific date.
      Correlation chart will update with historical
      correlations (time span will be hover date +/- selected window length)")
  ))
)


# Callbacks Start ------------------------------

app$callback(output = list(
  output(id = "correlation", property = "figure"),
  output(id = "correlation", property = "style")
),
params = list(
  input(id = "timeseries", property = "hoverData"),
  input(id = "window-length", property = "value")

),

  function(hover, window) {

    # Ensure hover is not empty
    if (!is.null(hover[[1]])) {

      datapoint <- as.numeric(hover$points[[1]]$pointNumber)
      rng <- (datapoint - window):(datapoint + window)
      cormat <- round(cor(stockdata[rng, 1:5]), 2)

      return(list(
        plot_ly(x = rownames(cormat),
                y = colnames(cormat),
                z = cormat,
                type = "heatmap",
                colors = colorRamp(c('#e3dfc8', '#808c6c')))%>%
          layout(title = "Correlation heatmap",
                 xaxis = list(title = ""),
                 yaxis = list(title = "")),
        list("display" = "inline")
        ))
    }
})


if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}
