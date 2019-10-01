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


# Load Data ------------------------------

ideal <- read.csv("Data/UN_IdealPoints.csv", stringsAsFactors = F)


# App Start ------------------------------

# Initiate application
app <- Dash$new()


# Create Layout ------------------------------

app$layout(
  htmlA(list(
    htmlImg(id = "banner-image", src = "assets/image.png")
  ), className = "logo",
  href = "https://dashr.plot.ly/"),
  htmlDiv(list(htmlH1("Ideal Points Estimation"))),
  htmlDiv(list(htmlH4("Country Name(s) of Interest"))),
  htmlDiv(list(
    htmlDiv(
      list(
        dccDropdown(
          id = "country-dropdown",
          className = "app-controls-block-dropdown",
          options = lapply(as.character(unique(ideal$Name)), function(x) {
            list(label = x, value = x)
          }),
          value = "United States of America",
          multi = TRUE,
          clearable = TRUE,
          searchable = TRUE
        ),
        htmlBr(),
        htmlDiv(list(dccGraph(id = "barplot",
                              style = list("height" = "200px")))),
        htmlBr(),
        htmlP(
          "Data: Bailey, Michael, Anton  Strezhnev and Erik Voeten. Forthcoming.
          'Estimating Dynamic State Preferences from United Nations Voting Data.'
          Journal of Conflict Resolution."
        )
      ),
      style = list("marginBottom" = 50),
      className = "four columns"
    ),
    htmlDiv(list(dccGraph(id = "lineplot")),
            className = "seven columns")
  ), className = "row")
)


# Callbacks Start ------------------------------

app$callback(
  output = list(id = "lineplot", property = "figure"),
  params = list(input(id = "country-dropdown", property = "value")),

  function(selected_country) {
    filtered_df <- ideal %>% filter(Name %in% selected_country)

    ggideal_point <- ggplot(filtered_df) +
      geom_line(aes(x = Year, y = Ideal.point, color = Name)) +
      labs(x = "Year", y = "Ideology", title = "Ideal Points for Countries") +
      scale_colour_hue("", l = 70, c = 150) + ggthemes::theme_few() +
      theme(legend.direction = "horizontal", legend.position = "bottom")

    # Convert ggplot object to plotly
    gg <- plotly_build(ggideal_point)

    return(gg)
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

    gg <- ggplot(filtered_df, aes(x = reorder(Name, terms), y = terms)) +
      geom_bar(stat = "identity", fill = "#2980b9") + coord_flip() +
      theme_bw() + labs(y = "Terms (in years)", x = "")

    p <- plotly_build(gg)
    return(p)
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv("PORT", 8050))
} else {
  app$run_server(showcase = TRUE)
}
