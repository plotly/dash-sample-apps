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

# compute a correlation matrix
correlation <- round(cor(mtcars), 3)
nms <- names(mtcars)

app <- Dash$new()

app$layout(htmlDiv(list(
  htmlDiv(htmlA(list(
    htmlImg(id = "banner-image", src = "assets/image.png")
  ), className = "logo",
  href = "https://dashr.plot.ly/"), className = 'row'),
  htmlDiv(list(
    htmlDiv(list(
      dccGraph(
        id = 'heatmapplot',
        figure = plot_ly(
          x = nms,
          y = nms,
          z = correlation,
          key = correlation,
          type = "heatmap",
          source = "heatplot"
        ) %>%
          layout(xaxis = list(title = ""),
                 yaxis = list(title = ""))
      )
    ), className = 'six columns'),
    
    htmlBr(),
    htmlDiv(list(dccGraph(id = 'scatterplot')), id = "scatterplot-container", 
            className = 'six columns')
  ), className = 'row'),
  htmlDiv(list(htmlDiv(id = "clickoutput", className = "shiny-container")),
          className = 'four columns')
)))

app$callback(
  output = list(id = "clickoutput", property = "children"),
  params = list(input(id = "heatmapplot", property = "clickData")),
  function(click_data) {
    if (length(click_data[[1]][[1]])) {
      paste("Curve Number:", click_data$points[[1]]$curveNumber,
              "x:", click_data$points[[1]]$x,
              "y:", click_data$points[[1]]$y,
              "z:", click_data$points[[1]]$z)
    } else {
      "Click on a cell in the heatmap to display a scatterplot"
    }
  }
)

app$callback(
  output = list(id = "scatterplot", property = "figure"),
  params = list(input(id = "heatmapplot", property = "clickData")),
  function(click_data) {
    if (length(click_data[[1]][[1]])) {
      vars <- c(click_data$points[[1]]$x, click_data$points[[1]]$y)
      d <- setNames(mtcars[vars], c("x", "y"))
      yhat <- fitted(lm(y ~ x, data = d))
      p <- plot_ly(d, x = ~ x) %>%
        add_markers(y = ~ y) %>%
        add_lines(y = ~ yhat) %>%
        layout(
          xaxis = list(title = click_data$points[[1]]$x),
          yaxis = list(title = click_data$points[[1]]$y),
          showlegend = FALSE
        )
    } else {
      p <- plotly_empty()
    }
    return(p)
  }
)

app$callback(
  output(id = "scatterplot-container", property = "style"),
  params = list(input(id = "heatmapplot", property = "clickData")),
  function(click_data) {
    if (length(click_data[[1]][[1]])) {
      style = list("display" = "block")
    }
    else {
      style = list("display" = "none")
    }
    return(style)
  }
)


if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(showcase = TRUE)
}
