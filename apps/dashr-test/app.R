library(dashR)
library(dashHtmlComponents)
library(dashCoreComponents)

css_list <- list("https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css",
                 "https://dash.plot.ly/assets/base.css?m=1552356970.0",
                 "https://dash.plot.ly/assets/daq.css?m=1552356970.0",
                 "https://dash.plot.ly/assets/dosis.css?m=1552356970.0",
                 "https://dash.plot.ly/assets/local-css-example.css?m=1552356970.0",
                 "https://dash.plot.ly/assets/open-sans.css?m=1552356970.0",
                 "https://dash.plot.ly/assets/override.css?m=1552356970.0",
                 "https://dash.plot.ly/assets/pdf-docs.css?m=1552356970.0",
                 "https://dash.plot.ly/assets/tabs-styled-with-classes.css?m=1552356970.0")

app <- Dash$new(external_stylesheets=css_list)
