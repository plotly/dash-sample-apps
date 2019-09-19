appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("/app/apps/%s", appName))
}

library(plotly)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(mlbench) # to load BreastCancer data
library(dplyr)
library(jsonlite) # to convert eventData to df


# Load & Prep Data ------------------------------

# Load data
data(BreastCancer)

# Remove NAs
BreastCancer <- na.omit(BreastCancer)

# Remeve ID
BreastCancer <- BreastCancer[, -1]

# Store features and actual class in seprate variables
featureList <- colnames(BreastCancer)[-10]
class <- BreastCancer$Class

# Convert to numeric
BreastCancer[, 1:9] <- apply(BreastCancer[, -10], 2, as.numeric)


# App Start ------------------------------

# Initiate application
app <- Dash$new(name = "dashr-shiny-coupledclickevents")


# Create Layout Variables ------------------------------

plotlyLogo <-
  htmlA(list(htmlImg(id = "banner-image", src = "assets/image.png")),
        className = "logo",
        href = "https://dashr.plot.ly")


featureOptions <- lapply(
  names(BreastCancer[1:9]),
  function(x) {
    list(label = x, value = x)
})

firstFeatureDropdown <- dccDropdown(
  id = "first-feature",
  options = featureOptions,
  value = names(BreastCancer)[1]
)

featureOptions <- lapply(
  names(BreastCancer[1:9]),
  function(x) {
    list(label = x, value = x)
  })

secondFeatureDropdown <- dccDropdown(
  id = "second-feature",
  options = featureOptions,
  value = names(BreastCancer)[5]
)


# Create Layout ------------------------------

app$layout(
  htmlDiv(list(
    plotlyLogo,
    htmlH2("Coupled events in plotly charts using Dash"),
    htmlDiv(list(
      htmlH4("This Dash app showcases coupled events using Dash's ",
             style = list("display" = "inline")),
      htmlCode(list("selectedData", " & ", "clickData"),
               style = list("display" = "inline")),
      htmlH6(" properties.", style = list("display" = "inline"))
    )),
    htmlBr(),
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv("1. The first chart showcases",
                style = list("display" = "inline")),
        htmlCode(" selectedData")
      )),
      htmlDiv(list(
        htmlDiv("2. The third chart showcases",
                style = list("display" = "inline")),
        htmlCode(" clickData")
      ), style = list("display" = "inline"))
    )),
    htmlHr(),
    htmlDiv(
      className = "row",
      style = list("width" = "600px"),
      list(
        htmlDiv(list(
        htmlStrong("Select first feature"),
        firstFeatureDropdown
        ), className = "six columns"),
        htmlDiv(list(
        htmlStrong("Select second feature (observed event)"),
        secondFeatureDropdown
        ), className = "six columns")
    )),
    htmlBr(),
    htmlDiv(
      className = "row",
      list(
        htmlDiv(
          dccGraph(id = "plot1"),
          className = "six columns"
          ),
        htmlDiv(
          dccGraph(id = "plot2"),
          className = "six columns"
        )
      )
    ),
    htmlHr(),
    htmlBlockquote("First drag a selection box in the scatter
      plot to populate the barchart.
      Then select one of the bars in the barchat
      to populate the boxplot."),
    htmlDiv(
      className = "row",
      list(
        htmlDiv(
          dccGraph(id = "plot3", style = list("display" = "none")),
          className = "six columns"
        ),
        htmlDiv(
          dccGraph(id = "plot4", style = list("display" = "none")),
          className = "six columns"
        )
      )
    )
  ))
)


# Callbacks Start ------------------------------

app$callback(output = list(
              output(id = "plot1", property = "figure"),
              output(id = "plot2", property = "figure")
             ),
             params = list(
               input(id = "first-feature", property = "value"),
               input(id = "second-feature", property = "value")
             ),

  function(feature1, feature2) {

    # Create a convenience data.frame which can be used for charting
    plotdf <- data.frame(BreastCancer[, feature1],
                       BreastCancer[, feature2],
                       Class = BreastCancer$Class)

    # Add column names
    colnames(plotdf) <- c("x", "y", "Class")

    contour1 <- plot_ly(plotdf, x = ~x, y = ~y,
                        mode = "markers",
                        type = "scatter",
                        color = ~Class,
                        source = "subset",
                        marker = list(size = 30)) %>%
                layout(title = paste(feature1, "vs ", feature2),
                       xaxis = list(title = feature1),
                       yaxis = list(title = feature2),
                       dragmode =  "select",
                       plot_bgcolor = "6A446F")

    contour2 <- plotdf %>%
      group_by(x, y, Class) %>%
      summarize(Count = n()) %>%
      filter(Class == "malignant") %>%
      plot_ly(x = ~x, y = ~y, z = ~Count, type = "contour") %>%
      layout(title = "Contour map of number of malignant cases",
             xaxis = list(title = feature1),
             yaxis = list(title = feature2))

    return(list(
      contour1,
      contour2
    ))
})


app$callback(output = list(
  output(id = "plot3", property = "style"),
  output(id = "plot3", property = "figure"),
  output(id = "plot4", property = "style"),
  output(id = "plot4", property = "figure")
  ),
params = list(
  input(id = "plot1", property = "selectedData"),
  input(id = "plot3", property = "clickData"),
  state(id = "first-feature", property = "value"),
  state(id = "second-feature", property = "value")
  ),

  function(selection, click, feature1, feature2) {

    # Ensure that selection is not empty
    if ( !is.null(selection[[1]]) ) {
      # Create a convenience data.frame which can be used for charting
      plotdf <- data.frame(BreastCancer[, feature1],
                            BreastCancer[, feature2],
                            Class = BreastCancer$Class)
      # Add column names
      colnames(plotdf) <- c("x", "y", "Class")

      selectionList <- fromJSON(toJSON(selection))
      # In case of `undefined columns selected` error can be ignored
      selectedDf <- selectionList[["points"]]
      selectedDf <- as.data.frame(apply(selectedDf, 2, as.integer))

      # Get number of malignant and benign cases from selection
      maligClass <- subset(plotdf, Class == "malignant")[subset(selectedDf, curveNumber == 0)$pointNumber + 1, ]
      benignClass <- subset(plotdf, Class == "benign")[subset(selectedDf, curveNumber == 1)$pointNumber + 1, ]
      # Combine
      plotSubset <- rbind(maligClass, benignClass)
      # Summarize
      plotsumm <- as.data.frame(
        plotSubset %>%
        group_by(x, y, Class) %>%
        summarize(Count = n())
      )

      barplot <- plot_ly(plotsumm,
                         x = ~Class,
                         y = ~Count,
                         type = "bar",
                         source = "select",
                         color = ~Class) %>%
        layout(title = "No. of Malignant and Benign cases <br> in Selection",
               plot_bgcolor = "6A446F",
               yaxis = list(domain = c(0, 0.9)))

      # If click is not null
      if ( !is.null(click[[1]]) ) {

        boxplotstyle <- list("display" = "inline")
        clickList <- fromJSON(toJSON(click))
        clickDf <- clickList[["points"]]

        # If Malignant
        if (clickDf$x[[1]] == "malignant") {

          tab <- subset(plotsumm, Class == "malignant")

          p1 <- plot_ly(tab,
                        x = ~x,
                        y = ~Count,
                        type = "box",
                        showlegend = F) %>%
            layout(yaxis = list(title = "Count"),
                   xaxis = list(title = feature1))

          p2 <- plot_ly(tab,
                        x = ~y,
                        y = ~Count,
                        type = "box",
                        showlegend = F) %>%
            layout(title = "Box plot for Malignant cases",
                   yaxis = list(title = "Count"),
                   xaxis = list(title = feature2))

          boxplot <- subplot(p1, p2)
        } else {
          tab <- subset(plotsumm, Class == "benign")

          p1 <- plot_ly(tab,
                        x = ~x,
                        y = ~Count,
                        type = "box",
                        showlegend = F) %>%
            layout(yaxis = list(title = "Count"),
                   xaxis = list(title = feature1))

          p2 <- plot_ly(tab,
                        x = ~y,
                        y = ~Count,
                        type = "box",
                        showlegend = F) %>%
            layout(title = "Box plot for Benign cases",
                   yaxis = list(title = "Count"),
                   xaxis = list(title = feature2))

          boxplot <- subplot(p1, p2)
        }

      # Click is null
      } else {
        boxplotstyle <- list("display" = "none")
        boxplot <- plotly_empty()

      }

      return(list(
        list("display" = "inline"),
        barplot,
        boxplotstyle,
        boxplot
      ))
    }
})


if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}
