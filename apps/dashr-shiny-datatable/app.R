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
library(dashTable)
library(jsonlite)
library(dplyr)


# Load & Prep Data ------------------------------

# Load data
data(mtcars)

m <- mtcars %>%
  tibble::rownames_to_column()

# Prevent data point from collapsing on scatter
m$mpg[1] <- m$mpg[1] + 3


# App Start ------------------------------

# Initiate application
app <- Dash$new(name = "dashr-shiny-datatable")


# Create Layout Variables ------------------------------

plotlyLogo <- htmlA(list(
  htmlImg(id = 'banner-image', src = 'assets/image.png')), className = 'logo',
        href = 'https://dashr.plot.ly')


carsDashDT <- dashDataTable(
  id = "cars-table",
  columns = lapply(colnames(m),
                   function(colName){
                     list(
                       id = colName,
                       name = colName
                     )
                   }),
  data = df_to_list(m),
  row_selectable = "multi",
  page_size = 10,
  sort_action="native",
  page_current= 0
)

carsScatter <- plot_ly(data = m,
                       x = ~mpg,
                       y = ~disp,
                       mode = "markers",
                       color = I('black'),
                       name = 'Unfiltered') %>%
               layout(showlegend = T)

downloadData <- htmlA(
  children = htmlButton(
    "Download Filtered Data",
    id = "download-data"
  ),
  href = "data/mtcars-filtered.csv",
  download = "data/mtcars-filtered.csv"
)


# Create Layout ------------------------------

app$layout(
  htmlDiv(list(
    plotlyLogo,
    htmlH1("Plotly & DashTable", style = list("textAlign" = "center")),
    dccGraph(id = "cars-scatter", figure = carsScatter),
    carsDashDT,
    htmlBr(),
    htmlDiv(downloadData, style = list("textAlign" = "center")),
    dccStore(id = "download-store")
  ))
)


# Callbacks Start ------------------------------

# Highlight dashTable rows
app$callback(output = list(
  output(id = "cars-table", property = "style_data_conditional")
),
params = list(
  #  `derived_virtual_selected_rows`
  #  dynamically change row indexes with sorted dashtable,
  # `selected_rows` remain static
  input(id = "cars-table", property = "derived_virtual_selected_rows"),
  input(id = "cars-table", property = "page_current"),
  input(id = "cars-scatter", property = "selectedData")
),

  function(rows, currentpage, selectedScatter) {

    style_data_conditional <- NULL

    # 1 - Handle table row selections
    # Row selection not empty list
    if (length(rows) > 0) {

      # List is not null
      if (!is.null(rows[[1]]) ) {

        for (i in 1:length(rows)) {

          # Deciding which rows to mark since table indexes reset on everypage
          rowIndex <- rows[[i]] - (currentpage) * 10

          style_data_conditional_append = list(list(
            "if" =  list("row_index" = rowIndex),
            "backgroundColor" = "#B0BED9"
          ))
          style_data_conditional <- c(style_data_conditional,
                                      style_data_conditional_append)
        }
      }
    }

    # 2 - Handle figure selections
    # Selection not empty list
    if (length(selectedScatter$points) > 0) {

      # List is not null
      if (!is.null(selectedScatter$points[[1]])) {

        selectedvector <- unlist(fromJSON(toJSON(selectedScatter))$points$pointIndex)
        unlistrows <- unlist(rows)

        # Prioritize highlighting table selections
        if (length(rows) > 0 && !is.null(rows[[1]])) {
          selectedvector <- selectedvector[!(selectedvector %in% unlistrows)]
        }
        for (i in 1:length(selectedvector)) {

          # Deciding which rows to mark since table indexes reset on everypage
          rowIndexScatter <- selectedvector[[i]] - (currentpage) * 10

          style_data_conditional_append_scat = list(list(
            "if" =  list("row_index" = rowIndexScatter),
            "color" = "white",
            "backgroundColor" = "black"
          ))
          style_data_conditional <- c(style_data_conditional,
                                      style_data_conditional_append_scat)
        }
      }
    }

    return(list(
      style_data_conditional
    ))
})


# Mark Selected points
app$callback(output = list(
  output(id = "cars-scatter", property = "figure")
),
params = list(
  input(id = "cars-scatter", property = "selectedData"),
  input(id = "cars-table", property = "selected_rows")
),

  function(selectedScatter, rows) {

    # 1 - Handle table selections
    # Row selection not empty list
    if (length(rows) > 0 ) {
      # List is not null
      if (!is.null(rows[[1]])) {
        rowvector <- unlist(rows) + 1
        m <- mutate(m, groupSelected = ifelse(
          row_number() %in% rowvector, "2", "1"))

        pp <- plot_ly(data = m,
                      x = m$mpg,
                      y = m$disp,
                      mode = "markers",
                      marker = list(color = factor(m$groupSelected,
                                                   labels = c("black","red"))),
                      name = "Filtered") %>%
          layout(showlegend = T)
      } else {
        pp <- carsScatter
      }
    } else {
      pp <- carsScatter
    }

    # 2 - Handle figure selections
    # Selection not empty list
    if (length(selectedScatter$points) > 0 ) {
      # List is not null
      if (!is.null(selectedScatter$points[[1]])) {

        selectedvector <- unlist(fromJSON(toJSON(selectedScatter))$points$pointIndex) + 1

        # Append rows if they exist
        if (exists("rowvector")) {
          selectedvector <- unique(c(rowvector, selectedvector))
        }
        m <- mutate(m, groupSelected = ifelse(
          row_number() %in% selectedvector, "2", "1"))

        pp <- plot_ly(data = m,
                      x = m$mpg,
                      y = m$disp,
                      mode = "markers",
                      marker = list(color = factor(m$groupSelected,
                                                   labels = c("black","red"))),
                      name = 'Filtered') %>%
          layout(showlegend = T)
      } else {
        # Handles the case when single point is unselected
        if (!exists("pp")) {
        pp <- carsScatter
        }
      }
    } else {
      if (!exists("pp")) {
        pp <- carsScatter
      }
    }
    return(list(pp))
})


# Write filtered data to csv
app$callback(
  output(id = "download-store", property = "children"),
params = list(
  input(id = "download-data", property = "n_clicks"),
  state(id = "cars-table", property = "selected_rows"),
  state(id = "cars-scatter", property = "selectedData")
  ),

  function(nclicks, rows, selectedScatter) {

    rowvector <- vector()
    scattervector <- vector()

    if (length(rows) > 0) {
      if (!is.null(rows[[1]]) ) {
        rowvector <- unlist(rows) + 1
      }
    }

    if (length(selectedScatter$points) > 0) {
      if (!is.null(selectedScatter$points[[1]])) {
        scattervector <- unlist(fromJSON(toJSON(selectedScatter))$points$pointIndex) + 1
      }
    }
    combinedvector <- unique(c(rowvector, scattervector))

    write.csv(m[combinedvector, ], "data/mtcars-filtered.csv")
})

app$run_server()
