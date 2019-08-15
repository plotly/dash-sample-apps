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

firstP <- htmlDiv(htmlLabel("Number of x bins"), htmlBr())

secondP <- htmlDiv(htmlLabel("Number of y bins"), htmlBr())

xSlider <- dccSlider(id = "x-slider",
                     min = 1,
                     max = 50,
                     value = 20,
                     width = 250)

ySlider <- dccSlider(id = "y-slider",
                     min = 1,
                     max = 50,
                     value = 20,
                     width = 250)


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
