appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(
    DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
    DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix
  )
  
  setwd(sprintf("/app/apps/%s", appName))
}

library(shiny)
library(plotly)
library(shinythemes)
library(dplyr)

######### LOAD DATA, CREATE FUNCTIONS & GLOBAL OBJECTS ##############

#Load Data

stockdata <- read.csv("https://cdn.rawgit.com/plotly/datasets/master/stockdata.csv")


################# LAYOUT VARIABLES ############

pageTitle <- htmlH2("Coupled hover-events in plotly charts using Shiny")
subTitle <- htmlH4("This Shiny app showcases coupled hover-events using Plotly's event_data() function.")


plotlyLogo <-
  htmlA(list(htmlImg(id = "banner-image", src = "assets/image.png")), 
        className = "logo", href = "https://dashr-docs.herokuapp.com/")

windowDropdown <- dccDropdown(
  id = "windowlength"
)

############## APP START ##################

app <- Dash$new()

################## CALLBACKS ##################




######### CONDITIONAL STATEMENT FOR APP RUNNING ON CLOUD SERVER & LOCAL ########

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}
