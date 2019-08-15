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


####################################################################################################



#################################### APP START #####################################################

app <- Dash$new(name = "dashr-shiny-coupledclickevents")
# Initiate application

####################################################################################################




#################################### CREATE LAYOUT VARIABLES #######################################

plotlyLogo <-
  htmlA(list(htmlImg(id = 'banner-image', src = 'assets/image.png')), className = 'logo',
        href = 'https://dashr.plot.ly')

pageTitle <- g

#################################### CREATE LAYOUT###################################################

app$layout(
  htmlDiv(list(
    plotlyLogo,
    pageTitle,
  ), className = "twelve columns"),
  
)



#################################### CALLBACKS START ###############################################

app$callback(
  output = list(id = , property =),
  params = list(input(id = , property = ))
)

####################################################################################################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
