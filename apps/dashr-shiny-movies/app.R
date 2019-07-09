appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf("/app/apps/%s", appName))
}

library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)

######### LOAD DATA, CREATE FUNCTIONS & GLOBAL OBJECTS ##############

install.packages("ggplot2movies")  #moved its own package to reduce the download size of ggplot2.
library(ggplot2movies)

minx <- min(movies$rating)
maxx <- max(movies$rating)
# size of the bins depend on the input 'bins'
size <- (maxx - minx) / input$bins

############################################# APP START ##########################

app <- Dash$new()

################################### LAYOUT VARIABLES ##########################

pageTitle <- htmlH1("Movie Ratings!")

plotlyLogo <- htmlA(
  list(
    htmlImg(
      src = "assets/dash-logo-new.png"
),
  href = "https://dashr-docs.herokuapp.com/"))

firstP <- htmlDiv(htmlLabel("Number of bins:"), htmlBr())

slider <- htmlDiv(
  dccSlider(
  id = "movies-slider",
  min = 1,
  max = 50,
  marks = c("1", "6", "11", "16", "21", "26", "31", "36", "41", "46", "50"),
  value = 10
)
  
)))

##################################################################################################
app$layout(

htmlDiv(
  list(
  plotlyLogo,
  firstP,
  dccSlider(id = "movies-slider")
  ),
  htmlDiv(
    (list(
    dccGraph(id = "histogram")
  )))))

################## CALLBACKS ##################

app$callbacks(

p <- plot_ly(movies, x = rating, autobinx = F, type = "histogram",
             xbins = list(start = minx, end = maxx, size = size))
)

################## CONDITIONAL STATEMENT FOR APP RUNNING ON CLOUD SERVER & LOCAL #########################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}



