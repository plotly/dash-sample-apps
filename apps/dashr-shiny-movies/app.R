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
      src = "assets/dashR-logo-stripe.png",
      className = "logo")),
  href = "https://dashr-docs.herokuapp.com/")

firstP <- htmlP("Number of bins:")

slider <- dccSlider(
  
)

##################################################################################################
app$layout(
  
  
)


################## CONDITIONAL STATEMENT FOR APP RUNNING ON CLOUD SERVER & LOCAL #########################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}



