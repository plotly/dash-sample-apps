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
  htmlDiv(list(htmlH1('Ideal Points Estimation'))),
  htmlDiv(list(htmlH4('Country Name(s) of Interest'))),
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
        htmlDiv(list(dccGraph(id = 'barplot'))),
        htmlBr(),
        htmlP(
          "Data: Bailey, Michael, Anton  Strezhnev and Erik Voeten. Forthcoming. 'Estimating Dynamic State Preferences from United Nations Voting Data.' Journal of Conflict Resolution."
        )
      ),
      style = list('marginBottom' = 50),
      className = "four columns"
    ),
    htmlDiv(list(dccGraph(id = 'lineplot')),
            className = 'seven columns')
  ), className = "row")
)

app$callback(
  output = list(id = "lineplot", property = "figure"),
  params = list(input(id = "country-dropdown", property = "value")),
  function(selected_country) {
    filtered_df <- ideal %>% filter(Name %in% selected_country)

    p <- plot_ly(
      filtered_df,
      x = ~ Year,
      y = ~ Ideal.point,
      color = ~ Name,
      hovertemplate = paste(
        'Year: %{x}',
        '<br>Ideal.point: %{y}',
        '<br>Name:',
        filtered_df$Name
      ),
      type = "scatter",
      mode = 'lines'
    ) %>%
      layout(
        title = "Ideal Points for Countries",
        yaxis = list(title = 'Ideology'),
        showlegend = FALSE
      )

    #Manually access and modify plot's underlying spec.
    pl <- plotly_build(p)$x
    pl$layout$annotations <- list()

    #Manually adding annotations
    for(i in 1:length(pl$data)){
      pl$layout$annotations[[i]] <- list(
        text = pl$data[[i]]$name,  # The text label of the annotation, e.g. "Canada"
        font = list(color = pl$data[[i]]$line$color), # Match the font color to the line color
        showarrow = FALSE, # Don't show the annotation arrow
        y = pl$data[[i]]$y[[length(pl$data[[i]]$y)]], # set the y position of the annotation to the last point of the line
        yref = "y1", # the "y" coordinates above are with respect to the yaxis
        x = 1, # set the x position of the graph to the right hand side of the graph
        xref = "paper", # the x coordinates are with respect to the "paper", where 1 means the right hand side of the graph and 0 means the left hand side
        xanchor = "left" # position the x coordinate with respect to the left of the text
      )
    }
    # increase the size of the right margin to accommodate more room for the annotation labels
    pl$layout$margin$r <- 170

    return(as_widget(pl))
  }
)

app$callback(
  output = list(id = "barplot", property = "figure"),
  params = list(input(id = "country-dropdown", property = "value")),
  function(selected_country) {
    filtered_df <- ideal %>%
      filter(Name %in% selected_country) %>%
      group_by(Name) %>%
      summarise(terms = n())

    p <- plot_ly(
      filtered_df,
      x = ~ terms,
      y = ~ Name,
      type = "bar",
      orientation = 'horizontal'
    ) %>%
      layout(xaxis = list(title = 'Terms (in years)'),
             yaxis = list(title = ''))
    return(p)
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(showcase = TRUE)
}
