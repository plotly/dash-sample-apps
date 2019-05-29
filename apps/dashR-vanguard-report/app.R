#install.packages('devtools')
#
library(devtools) # devtools: Tools to Make Developing R Packages Easier
# install_github('plotly/dashR', ref= "0.0.7-debug") # The core dash backend
# install_github('plotly/dash-html-components') # HTML components
# install_github('plotly/dash-core-components') # Supercharged components

# Install a package that is necessary and has not been installed, or load one that is.

usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}


setwd('plotly/dash-sample-apps/apps/dashR-vanguard-report')
#Source assets
source("assets/VanguardFunctions.R")

#setwd("Vanguard")


# Load Necessary Packages

usePackage('dashR')

usePackage('dashCoreComponents')

usePackage('dashHtmlComponents')

usePackage('plyr')

usePackage('plotly')

usePackage('lubridate')

usePackage('dashTable')

usePackage('taRifx')


#################################################################################

external_css1 = list("https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",
                     "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
                     "https://codepen.io/bcd/pen/KQrXdb.css",
                     "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
                     "https://cdn.rawgit.com/plotly/dash-app-stylesheets/5047eb29e4afe01b45b27b1d2f7deda2a942311a/goldman-sachs-report.css")


app <- Dash$new(external_stylesheets = external_css1)





