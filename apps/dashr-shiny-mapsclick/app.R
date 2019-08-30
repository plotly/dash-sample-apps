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
      id = 'mapplot',
      figure = plot_ly(
        z = state.area,
        text = state.name,
        locations = state.abb,
        type = 'choropleth',
        locationmode = 'USA-states'
      ) %>%
        layout(geo = list(
          scope = 'usa',
          projection = list(type = 'albers usa'),
          lakecolor = toRGB('white')
        ))
    ),
    htmlDiv(id = "clickdataOutput")
  ))
)

app$callback(
  output = list(id = "clickdataOutput", property = "children"),
  params = list(input(id = "mapplot", property = "clickData")),
  function(click_data) {
    if (is.null(click_data[[1]])) {
      d <- "Click on a state to view event data"
    } else {
      d <- paste(
        "Curve Number:", click_data$points[[1]]$curveNumber,
        "Point Number:", click_data$points[[1]]$pointNumber,
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
