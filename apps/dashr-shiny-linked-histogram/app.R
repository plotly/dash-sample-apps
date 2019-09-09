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


# Load data ---------------------------

data(cars)


# Helper function(s) ---------------------------

# Generates min, max, size for xbins/ybins attributes of plot_ly histogram
compute_bins <- function(x, n) {
  list(
    start = min(x),
    end = max(x),
    size = (max(x) - min(x)) / n
  )
}


# Helper variables ---------------------------

# Marker objects
m <- list(color = toRGB("black"))
m2 <- list(color = toRGB("black", 0.2))


# App Start ---------------------------

# Initiate application
app <- Dash$new(name = "dashr-shiny-linked-histogram")


# Create Layout Variables ---------------------------

plotlyLogo <- htmlA(
  list(htmlImg(
  id = 'banner-image', src = 'assets/image.png')),
  className = 'logo',
  href = 'https://dashr.plot.ly')

pageTitle <- htmlH2("Linked highlighting with plotly and Dash for R")

firstP <- htmlB("Number of x bins")

secondP <- htmlB("Number of y bins")

xSlider <- dccSlider(
  id = "x-slider",
  min = 1,
  max = 50,
  marks = setNames(c(as.list(seq(1, 50, 5)), as.list(50)),
                  c(seq(1, 50, 5), 50)),
  value = 20
)

ySlider <- dccSlider(
  id = "y-slider",
  min = 1,
  max = 50,
  marks = setNames(c(as.list(seq(1, 50, 5)), as.list(50)),
                   c(seq(1, 50, 5), 50)),
  value = 20
)

xHistogram <- dccGraph(
  id = 'x-histogram',
  p <- plot_ly(x = cars$speed, type = "histogram", autobinx = FALSE,
               xbins = compute_bins(cars$speed, length(unique(cars$speed))), marker = m2) %>%
         layout("yaxis" = list("title" = "count")),
  figure = p
)

yHistogram <- dccGraph(
  id = 'y-histogram',
  p <- plot_ly(y = cars$dist, type = "histogram", autobiny = F,
               ybins = compute_bins(cars$dist, length(unique(cars$dist))), marker = m2) %>%
         layout("xaxis" = list("title" = "count")),
  figure = p
  )

scatterGraph <- dccGraph(
  id = 'scatter-graph',
  figure = plot_ly(cars, x = ~speed, y = ~dist,
                   mode = "markers", marker = m) %>%
    layout(dragmode = "select")
)


# Create Layout ---------------------------

app$layout(htmlDiv(
  list(
    plotlyLogo,
    pageTitle,
    htmlBr(),
    htmlDiv(
      className = "row",
      style = list("width" = "650px"),
      list(
        htmlDiv(xHistogram,
                className = "eight columns" ),
        htmlDiv(list(firstP, xSlider, htmlBr(), secondP, ySlider),
                className = 'four columns')
      )
    ),
    htmlDiv(
      className = "row",
      style = list("width" = "650px"),
      list(
        htmlDiv(scatterGraph,
                className = "seven columns"),
        htmlDiv(
          yHistogram,
          #style = list('marginBottom'= 50),
          className = 'five columns')
      )
    )
  )
))


# Callbacks Start ---------------------------

# the 'x' histogram
app$callback(
  output = list(id = 'x-histogram', property = 'figure'),
  params = list(input(id = 'x-slider', property = 'value'),
                input(id = 'scatter-graph', property = 'selectedData')),

  function(xbins, selected) {
    p <- plot_ly(x = cars$speed, type = "histogram", autobiny = F,
                 xbins = compute_bins(cars$speed, xbins), marker = m2)

    # if points are selected, subset the data, and highlight
    if (!is.null(selected[[1]])) {
      p <- add_trace(p, x = fromJSON(toJSON(selected$points))$x, type = "histogram", autobiny = F,
                     xbins = compute_bins(cars$speed, xbins), marker = m)
      p <- config(p, displayModeBar = F, showLink = F) %>%
           layout(showlegend = F, barmode = "overlay",
             yaxis = list(title = "count"),
             xaxis = list(title = "", showticklabels = F))
    }
    return(p)
})

# the 'y' histogram
app$callback(
  output = list(id = 'y-histogram', property = 'figure'),
  params = list(input(id = 'y-slider', property = 'value'),
                input(id = 'scatter-graph', property = 'selectedData')),

  function(ybins, selected) {
    p <- plot_ly(y = cars$dist, type = "histogram", autobiny = F,
                 ybins = compute_bins(cars$dist, ybins), marker = m2)

    # if points are selected, subset the data, and highlight
    if (!is.null(selected[[1]])) {
      p <- add_trace(p, y = fromJSON(toJSON(selected$points))$y, type = "histogram", autobiny = F,
                     ybins = compute_bins(cars$dist, ybins), marker = m)
      p <- config(p, displayModeBar = F, showLink = F) %>%
        layout(showlegend = F, barmode = "overlay",
               yaxis = list(title = "count"),
               xaxis = list(title = "", showticklabels = F))
    }
    return(p)
  })


if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}
