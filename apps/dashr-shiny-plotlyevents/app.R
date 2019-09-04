appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("/app/apps/%s", appName))
}

#setwd("/Users/Caner/Desktop/plotly/dashr-shiny-examples/dashr-shiny-plotlyevents")

library(plotly)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(jsonlite)

#################################### LOAD DATA & CREATE GLOBAL OBJECTS #############################

#Load dataset - mtcars
data("mtcars")

####################################################################################################



#################################### APP START #####################################################

app <- Dash$new(name = "dashr-shiny-plotlyevents")
# Initiate application

####################################################################################################




#################################### CREATE LAYOUT VARIABLES #######################################

plotlyLogo <- htmlA(
  list(htmlImg(
    id = 'banner-image', src = 'assets/image.png')),
  className = 'logo',
  href = 'https://dashr.plot.ly')

pageTitle <- htmlB("Plot Type")

#radiobuttons
radioButton <- dccRadioItems(
  id = "radiobutton-selector",
  options = list(
  list("label" = "ggplotly", "value" = "ggplotly"),
  list("label" = "plotly", "value" = "plotly")
), value = "ggplotly")

#################################### CREATE LAYOUT###################################################

app$layout(
  htmlDiv(list(
    plotlyLogo,
    pageTitle,
    radioButton,
    dccGraph(id = "scatter-plot", style = list("width" = "calc(100vw - 50px)")),
    htmlP(id = "hover-data", className = "shiny-container"),
    htmlP(id = "click-data", className = "shiny-container"),
    htmlP(id = "selected-data", className = "shiny-container"),
    htmlP(id = "relayout-data", className = "shiny-container")
    )),
    htmlDiv(
      dccMarkdown(children = paste("```r", "\n",
                                   readChar("app.R", file.info("app.R")$size), "\n",
                                   "```", sep = "")
      ),
      className = "container",
      style = list("maxWidth" = "650px", "borderLeft" = "thin solid lightgrey")
    )
)


#################################### CALLBACKS START ###############################################

app$callback(
  output = list(id = "scatter-plot", property = "figure"),
  params = list(input(id = "radiobutton-selector", property = "value")),

    function(selection) {
      key <- row.names(mtcars)
      if (identical(selection, "ggplotly")) {
        p <- ggplot(mtcars,
                    aes(x = mpg,
                        y = wt,
                        colour = factor(vs),
                        key = key)) + geom_point()

        ggplotly(p) %>% layout(dragmode = "select")
      } else {
        plot_ly(mtcars, x = ~mpg, y = ~wt, key = ~key) %>%
          layout(dragmode = "select")
      }
    }
)


app$callback(output = list(id = 'hover-data', property = 'children'),
             params = list(input(id = 'scatter-plot', property = 'hoverData')),
             function(hoverData) {
               d <- hoverData[[1]]
               if (is.null(d)) {
                 "Hover events appear here (unhover to clear)"
               } else {
               return(prettify(toJSON(hoverData), indent = 2))
               }
             })


app$callback(output = list(id = 'click-data', property = 'children'),
             params = list(input(id = 'scatter-plot', property = 'clickData')),

             function(clickData) {
               d <- clickData[[1]]
               if (is.null(d)) {
                 "Click events appear here (double-click to clear)"
               } else {
                 return(prettify(toJSON(clickData), indent = 2))
               }
             })


app$callback(output = list(id = 'selected-data', property = 'children'),
             params = list(input(id = 'scatter-plot', property = 'selectedData')),
             function(selectedData) {
               d <- selectedData[[1]]
               if (is.null(d)) {
                 "Click and drag events (i.e., select/lasso) appear here (double-click to clear)"
               } else {
                 return(prettify(toJSON(selectedData), indent = 2))
               }
             })


app$callback(output = list(id = 'relayout-data', property = 'children'),
             params = list(input(id = 'scatter-plot', property = 'relayoutData')),
             function(relayoutData) {
               d <- relayoutData[[1]]
               if (is.null(d) || (length(relayoutData) == 1 && isTRUE(relayoutData$autosize))) {
               "Relayout (i.e., zoom) events appear here"
               } else {
                 return(prettify(toJSON(relayoutData), indent = 2))
               }
             })


####################################################################################################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}
