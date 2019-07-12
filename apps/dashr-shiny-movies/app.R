appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf("/app/apps/%s", appName))
}

setwd("/Users/milaroisin/Documents/GitHub/plot.ly/dash-sample-apps/apps/dashr-shiny-movies")

library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)

######### LOAD DATA, CREATE FUNCTIONS & GLOBAL OBJECTS ##############

#install.packages("ggplot2movies") 
#include in its own init.R file
#moved its own package to reduce the download size of ggplot2.

library(ggplot2movies)
minx <- min(movies$rating)
maxx <- max(movies$rating)

############################################# APP START ##########################

app <- Dash$new()

################################### LAYOUT VARIABLES ##########################

pageTitle <- htmlH1("Movie Ratings!")

plotlyLogo <- htmlA(
  list(
    htmlImg(
      src = "assets/image.png")), className = "logo", 
    href = "https://dashr-docs.herokuapp.com/")

firstP <- htmlDiv(htmlLabel("Number of bins:"), htmlBr())

slider <- dccSlider(
  id = "movies-slider",
  min = 1,
  max = 50,
  marks = list('1' = '1','6'='6', '11' = '11', '16' = '16', '21' = '21', '26' = '26', '31' = '31', '36' = '36', '41' = '41', '46' = '46', '50' = '50'),
  value = 10
)

##################################################################################################
app$layout(htmlDiv(list(
  plotlyLogo,
  pageTitle,
  firstP,
  slider,
  # dccSlider(id = "movies-slider"),
  htmlBr(),
  
  htmlDiv(list(dccGraph(id = "histogram")))
  
)))

################## CALLBACKS ##################

# app$callbacks(
# 
#   
# p <- plot_ly(movies, x = rating, autobinx = F, type = "histogram",
#              xbins = list(start = minx, end = maxx, size = size))
# )

app$callback(output = list(id = "histogram", property = "figure"),
             params = list(
               input(id = "movies-slider", property = "value")
              ),
             
             function(sliderval){
               size <- (maxx - minx) / input$bins
               
               # a simple histogram of movie ratings
               p <- plot_ly(movies, x = rating, autobinx = F, type = "histogram",
                            xbins = list(start = minx, end = maxx, size = size)),
               list(
                 layout = list(
                   xaxis = list(title = "Ratings",
                                range = c(minx, maxx), 
                                autotick = F,
                                tick0 = minx, 
                                dtick = size)
                 )
               )
             }
)

################## CONDITIONAL STATEMENT FOR APP RUNNING ON CLOUD SERVER & LOCAL #########################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}



