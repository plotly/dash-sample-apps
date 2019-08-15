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

firstP <- htmlDiv(htmlLabel("Number of x bins"), htmlBr())

secondP <- htmlDiv(htmlLabel("Number of y bins"), htmlBr())

xSlider <- dccSlider(id = "x-slider",
                     min = 1,
                     max = 50,
                     value = 20)

ySlider <- dccSlider(id = "y-slider",
                     min = 1,
                     max = 50,
                     value = 20)

#################################### CREATE LAYOUT###################################################

app$layout(htmlDiv(
  list(
    plotlyLogo,
    pageTitle,
    firstP,
    htmlDiv(list(xSlider), style = ('marginBottom' = 50), className = 'five columns'),
    secondP,
    htmlDiv(list(ySlider), style = ('marginBottom' = 50), className = 'seven columns')
    )
  )
)

#################################### CALLBACKS START ###############################################

app$callback(
  output = list(id = "histogram", property = "figure"),
  params = list(input(id = , property = ))
)

####################################################################################################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
