appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)
  Sys.setenv(
    DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
    DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix
  )
  setwd(sprintf("/app/apps/%s", appName))
}

library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)

app <- Dash$new()

app$layout(
  htmlA(list(
    htmlImg(id = "banner-image", src = "assets/image.png")
  ), className = "logo",
  href = "https://dashr.plot.ly/"),
  htmlDiv(list(
    dccGraph(
      id = '3dplot',
      figure = plot_ly(
        x = rnorm(10),
        y = rnorm(10),
        z = rnorm(10),
        type = "scatter3d"
      ),
      clear_on_unhover = TRUE
    ),
    htmlDiv(id = "hoverdataOutput", className = "shiny-container"),
    htmlBr(),
    htmlDiv(id = "clickdataOutput", className = "shiny-container")
  ))
)

app$callback(
  output = list(id = "hoverdataOutput", property = "children"),
  params = list(input(id = "3dplot", property = "hoverData")),
  function(hover_data) {
    if (is.null(hover_data[[1]])) {
      d <- "Hover events appear here (unhover to clear)"
    } else {
      d <- paste("Curve Number:", hover_data[[1]][[1]]$curveNumber,
            "Point Number:", hover_data[[1]][[1]]$pointNumber,
            "x:", hover_data[[1]][[1]]$x,
            "y:", hover_data[[1]][[1]]$y,
            "z:", hover_data[[1]][[1]]$z)
    }
    return(d)
  }
)


app$callback(
  output = list(id = "clickdataOutput", property = "children"),
  params = list(input(id = "3dplot", property = "clickData")),
  function(click_data) {
    if (is.null(click_data[[1]])) {
      d <- "Click events appear here (double-click to clear)"
    } else {
      d <- paste(
        "Curve Number:", click_data$points[[1]]$curveNumber,
        "Point Number:", click_data$points[[1]]$pointNumber,
        "x:", click_data$points[[1]]$x,
        "y:", click_data$points[[1]]$y,
        "z:", click_data$points[[1]]$z
      )
    }
    return(d)
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(showcase = TRUE)
}
