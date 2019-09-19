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
library(jsonlite) # to convert selectedData to df

# Load & Prep Data ------------------------------

# Load data
data(BreastCancer)

# Remove NAs
BreastCancer <- na.omit(BreastCancer)

# Remeve ID
BreastCancer <- BreastCancer[,-1]

# Store features and actual class in seprate variables
featureList <- colnames(BreastCancer)[-10]
class <- BreastCancer$Class

# Convert to numeric
BreastCancer[, 1:9] <- apply(BreastCancer[,-10], 2, as.numeric)

####################################################################################################


# App Start ------------------------------

# Initiate application
app <- Dash$new(name = "dashr-shiny-coupledclickevents")


####################################################################################################


# Create Layout Variables ------------------------------

plotlyLogo <-
  htmlA(list(htmlImg(id = 'banner-image', src = 'assets/image.png')), className = 'logo',
        href = 'https://dashr.plot.ly')


featureOptions <- lapply(
  names(BreastCancer[1:9]),
  function(x) {
    list(label = x, value = x)
})

firstFeatureDropdown <- dccDropdown(
  id = "first-feature",
  options = featureOptions,
  value = names(BreastCancer)[1]
  #className = "six columns"
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
  #className = "six columns"
)


# Create Layout ------------------------------

app$layout(
  htmlDiv(list(
    plotlyLogo,
    htmlH2("Coupled events in plotly charts using Shiny"),
    htmlDiv(list(
      htmlH4("This Shiny app showcases coupled events using Plotly's ",
             style = list("display" = "inline")),
      htmlCode("event_data()", style = list("display" = "inline")),
      htmlH6(" function.", style = list("display" = "inline"))
    )),
    htmlBr(),
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv("1. The first chart showcases", style = list("display" = "inline")),
        htmlCode(" plotly_selected")
      )),
      htmlDiv(list(
        htmlDiv("2. The third chart showcases", style = list("display" = "inline")),
        htmlCode(" plotly_click")
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
    plot.df <- data.frame(BreastCancer[, feature1],
                       BreastCancer[, feature2],
                       Class = BreastCancer$Class)

    # Add column names
    colnames(plot.df) <- c("x", "y", "Class")

    contour1 <- plot_ly(plot.df, x = ~x, y = ~y,
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

    contour2 <- plot.df %>%
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
  output(id = "plot3", property = "figure")
  ),
params = list(
  input(id = "plot1", property = "selectedData"),
  state(id = "first-feature", property = "value"),
  state(id = "second-feature", property = "value")

  ),

  function(selection, feature1, feature2) {

    # Ensure that selection is not empty
    if ( !is.null(selection[[1]]) ) {
      # Create a convenience data.frame which can be used for charting
      plot.df <- data.frame(BreastCancer[, feature1],
                            BreastCancer[, feature2],
                            Class = BreastCancer$Class)
      # Add column names
      colnames(plot.df) <- c("x", "y", "Class")

      selectionList <- fromJSON(toJSON(selection))
      # In case of `undefined columns selected` error can be ignored
      selectedDf <- selectionList[["points"]]
      selectedDf <- as.data.frame(apply(selectedDf, 2, as.integer))

      # Get number of malignant and benign cases from selection
      malig.class <- subset(plot.df, Class == "malignant")[subset(selectedDf, curveNumber == 0)$pointNumber + 1,]
      benign.class <- subset(plot.df, Class == "benign")[subset(selectedDf, curveNumber == 1)$pointNumber + 1,]
      # Combine
      plot.subset <- rbind(malig.class, benign.class)
      # Summarize
      plot.summ <- plot.subset %>%
        group_by(x, y, Class) %>%
        summarize(Count = n())

      barplot <- plot_ly(plot.summ, x = ~Class, y = ~Count, type = "bar", source = "select", color = ~Class) %>%
        layout(title = "No. of Malignant and Benign cases <br> in Selection",
               plot_bgcolor = "6A446F",
               yaxis = list(domain = c(0, 0.9)))
      return(list(
        list("display" = "inline"),
        barplot
      ))
    }
})


####################################################################################################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}
