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

#Load dataset - mtcars
data("mtcars")

####################################################################################################



#################################### APP START #####################################################

app <- Dash$new(name = "dashr-shiny-plotlyevents")
# Initiate application

####################################################################################################




#################################### CREATE LAYOUT VARIABLES #######################################

plotlyLogo <-
  htmlA(list(htmlImg(id = 'banner-image', src = 'assets/image.png')), className = 'logo',
        href = 'https://dashr.plot.ly')

pageTitle <- htmlB("Plotly Type")

#radiobuttons
radioButton <- dccRadioItems(
  id = "radiobutton-selector",
  options = list(
  list("label" = "ggplotly", "value" = "ggp"),
  list("label" = "plotly", "value" = "plty")
), value = "ggp")

#################################### CREATE LAYOUT###################################################

app$layout(
  htmlDiv(list(
    plotlyLogo,
    pageTitle,
    radioButton
  ), className = "twelve columns"),
  
  htmlDiv(list(
    dccGraph( id = "scatter-plot",)
  )
  )
  
)



#################################### CALLBACKS START ###############################################

app$callback(
  output = list(id = "scatterplot", property = "figure"),
  params = list(input(id = "radiobutton-selector", property = "value"))
)

####################################################################################################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
