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


#################################### LOAD DATA & CREATE GLOBAL OBJECTS #############################

#cars {dataset - speed and stopping distances of cars}
data(cars)
head(cars)
summary(cars)

#convenience function for computing xbin/ybin object given a number of bins

compute_bins <- function(x, n) {
  list(
    start = min(x),
    end = min(x),
    size = (max(x) - min(x)) / n)
}


x <- cars$speed

y <- cars$dist

#marker objects 
m <- list(color = toRGB("black"))
m2 <- list(color = toRGN("black", 0.2))

####################################################################################################



#################################### APP START #####################################################

app <- Dash$new(name = "dashr-shiny-linked-histogram")
# Initiate application

####################################################################################################

#################################### CREATE LAYOUT VARIABLES #######################################

plotlyLogo <-
  htmlA(list(htmlImg(id = 'banner-image', src = 'assets/image.png')), className = 'logo',
        href = 'https://dashr.plot.ly')

pageTitle <- htmlH2("Linked highlighting with plotly and DashR")

firstP <- htmlB("Number of x bins")

secondP <- htmlB("Number of y bins")

xSlider <- dccSlider(id = "x-slider",
                     min = 1,
                     max = 50,
                     marks = list(
                       '1' = '1',
                       '6' = '6',
                       '11' = '11',
                       '16' = '16',
                       '21' = '21',
                       '26' = '26',
                       '31' = '31',
                       '36' = '36',
                       '41' = '41',
                       '46' = '46',
                       '50' = '50'
                       ),
                     value = 20)

ySlider <- dccSlider(
  id = "y-slider",
  min = 1,
  max = 50,
  marks = list(
    '1' = '1',
    '6' = '6',
    '11' = '11',
    '16' = '16',
    '21' = '21',
    '26' = '26',
    '31' = '31',
    '36' = '36',
    '41' = '41',
    '46' = '46',
    '50' = '50'
  ),
  value = 20
)

xHistogram -> dccGraph(
  id = 'x-histogram',
  figure = list(
    data = cars,
    type = 'bar',
  )
)

scatterGraph -> dccGraph(
  id = 'scatter-graph',
  figure = list(
    data = cars,
    layout = list(
      xaxis = list ('title' = 'speed'),
      yaxis = list('title' = 'dist'),
    )
  )
  
)

yHistogram -> dccGraph(
  id = 'y-histogram',
  figure = list(
    data = cars,
    type = 'bar',
  )
  
)
#################################### CREATE LAYOUT###################################################

app$layout(htmlDiv(
  list(
    plotlyLogo,
    pageTitle,
    xHistogram,
    firstP,
    htmlBr(),
    htmlDiv(list(xSlider), style=list('marginBottom'= 50, className = 'five columns')),
    secondP,
    htmlBr(),
    htmlDiv(list(ySlider), style=list('marginBottom'= 50, className = 'seven columns')),
    scatterGraph,
    yHistogram,
  )
))

#################################### CALLBACKS START ###############################################

app$callback(
  output = list(id = 'x-histogram', property = 'figure'),
  params = list(
        input(id = 'x-slider', property = 'value')),

function(xbins) {
xbins <- compute_bins(x)
p <- plot_ly(x = x, type = "histogram", autobinx = F,
             xbins = xbins, marker = m2)

# obtain plotlyjs selection

s <- event_data("plotly_selected")

# if points are selected, subset the data, and highlight
if (length(s) > 0) {
  p <- add_trace(p, x = s$x, type = "histogram", autobinx = F,
                 xbins = xbins, marker = m)
}
p %>%
  config(displayModeBar = F, showLink = F) %>%
  layout(showlegend = F, barmode = "overlay", yaxis = list(title = "count"),
         xaxis = list(title = "", showticklabels = F))
})

app$callback(
    output = list(id = 'y-histogram', property = 'figure'),
    input = list(id = 'y-slider', property = 'value'),

function(ybins)  {

ybins <- compute_bins(y, input$ybins)
p <- plot_ly(y = y, type = "histogram", autobiny = F,
             ybins = ybins, marker = m2)
s <- event_data("plotly_selected")
if (length(s$y) > 0) {
  p <- add_trace(p, y = y, type = "histogram", autobiny = F,
                 ybins = ybins, marker = m)
}
p %>%
  config(displayModeBar = F, showLink = F) %>%
  layout(showlegend = F, barmode = "overlay", xaxis = list(title = "count"),
         yaxis = list(title = "", showticklabels = F))
})

 

####################################################################################################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
