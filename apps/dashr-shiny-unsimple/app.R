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


ideal <- read.csv("Data/UN_IdealPoints.csv", stringsAsFactors = F)

app <- Dash$new()

app$layout(
  htmlA(list(
    htmlImg(id = "banner-image", src = "assets/image.png")
  ), className = "logo",
  href = "https://dashr.plot.ly/"),
  htmlDiv(list(htmlH1('Ideal Points'))),
  htmlDiv(list(
    htmlDiv(
      list(
        dccDropdown(
          id = 'country-dropdown',
          className = 'app-controls-block-dropdown',
          options = lapply(as.character(unique(ideal$Name)), function(x) {
            list(label = x, value = x)
          }),
          value = "United States of America",
          multi = TRUE,
          clearable = TRUE,
          searchable = TRUE
        ),
        htmlBr(),
        htmlP(
          "Data: Bailey, Michael, Anton  Strezhnev and Erik Voeten. Forthcoming. 'Estimating Dynamic State Preferences from United Nations Voting Data.' Journal of Conflict Resolution."
        )
      ),
      style = list('marginBottom' = 50),
      className = "four columns"
    ),
    htmlDiv(list(dccGraph(id = 'lineplot')),
            className = 'eight columns')
  ), className = "row")
)

app$callback(
  output = list(id = "lineplot", property = "figure"),
  params = list(input(id = "country-dropdown", property = "value")),

  function(selected_country) {
    filtered_df <- ideal %>% filter(Name %in% selected_country)
    print(selected_country)

    p <- plot_ly(
      filtered_df,
      x = ~Year,
      y = ~Ideal.point,
      color = ~Name,
      hovertemplate = paste('Year: %{x}',
                            '<br>Ideal.point: %{y}',
                            '<br>Name:', filtered_df$Name),
      type = "scatter",
      mode = 'lines'
    ) %>%
      layout(
        title = "Ideal Points for Countries",
        yaxis = list(title = 'Ideology')
      )
    return(p)
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(showcase = TRUE)
}
