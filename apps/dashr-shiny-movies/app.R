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

#################################### LOAD DATA & CREATE GLOBAL OBJECTS #############################
install.packages("ggplot2movies")  #moved its own package to reduce the download size of ggplot2.
library(ggplot2movies)

#####################################################################################################


