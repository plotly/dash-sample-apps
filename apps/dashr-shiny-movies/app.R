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
library(ggplot2movies) # contains movies data


# Load & Prep Data ------------------------------

moviesDf <- as.data.frame(ggplot2movies::movies)


# App Start ------------------------------

# Initiate application
app <- Dash$new()


# Create Layout Variables ------------------------------

pageTitle <- htmlH2("Movie Ratings!")

plotlyLogo <- htmlA(
  list(htmlImg( id = "banner-image", src = "assets/image.png")),
  className = "logo",
  href = "https://dashr-docs.herokuapp.com/")


firstP <- htmlDiv(htmlLabel("Number of bins:"))


slider <- dccSlider(
  id = "movies-slider",
  min = 1,
  max = 50,
  marks = setNames(c(as.list(seq(1, 50, 5)), as.list(50)),
                   c(seq(1, 50, 5), 50)),
  value = 10
)



# Create Layout ------------------------------

app$layout(htmlDiv(
  list(
    plotlyLogo,
    pageTitle,
    firstP,
    htmlDiv(list(
      htmlDiv(list(slider),
              style = list("marginBottom" = 50),
              className = "three columns"),
      htmlDiv(list(dccGraph(id = "histogram")),
              className = "six columns")),
      className = "row"
    )
  )
))


# Callbacks Start ------------------------------

app$callback(output = list(id = "histogram", property = "figure"),
             params = list(input(id = "movies-slider", property = "value")),

  function(bins) {

    minx <- min(movies$rating)
    maxx <- max(movies$rating)
    size <- (maxx - minx) / bins

    # Simple histogram of movie ratings
    p <- plot_ly(
      movies,
      x = movies$rating,
      autobinx = FALSE,
      type = "histogram",
      xbins = list(
        start = minx,
        end = maxx,
        size = size
      )
    ) %>% layout(xaxis = list(
      title = "Ratings",
      range = c(minx, maxx),
      autotick = FALSE,
      tick0 = minx,
      dtick = size
      ))
    return(p)
})

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv("PORT", 8050))
} else {
  app$run_server()
}
