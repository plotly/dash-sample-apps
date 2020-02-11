appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashTable)
library(dplyr)
source("utils/utils.R")

#################################### LOAD DATA & CREATE GLOBAL OBJECTS #############################

pkData <- read.csv("./data/pkdata.csv", stringsAsFactors = FALSE)
# Read df

####################################################################################################

#################################### APP START #####################################################

app <- Dash$new()
# Initiate application

####################################################################################################

#################################### CREATE LAYOUT VARIABLES #######################################

nSubjects <- length(unique(pkData$subject_index))
nTimes <- length(unique(pkData$time))
# Default values for cols & rows

dashBioLogo <- htmlImg(src = "./assets/dashbio_logo_transparent.png")

githubLogo <- htmlImg(src = "./assets/GitHub-Mark-Light-64px.png",
                      className = "github-logo")

timePointsLabel <- htmlLabel(
  list(htmlDiv(list("Time points")),
       dccInput(
         id = "times-input",
         placeholder = "Enter a value...",
         type = "number",
         value = nTimes,
         min = 3,
         max = 999
       )
  )
)

subjectsLabel <- htmlLabel(
  list(htmlDiv(list("Subjects")),
       dccInput(
         id = "subjects-input",
         placeholder = "Enter a value...",
         type = "number",
         value = nSubjects,
         min = 1,
         max = 48
       )
  )
)

styleConditional <- list(
  list(
    "if" = list("column_id" = "parameter"),
    "textAlign" = "right",
    "paddingRight" = 10
  ),
  list(
    "if" = list("row_index" = "odd"),
    "backgroundColor" = "white"
  )
)

tableHeaderStyle <- list(
  "backgroundColor" = "rgb(2,21,70)",
  "color" = "white",
  "textAlign" = "center"
)

pkDataTable <- dashDataTable(
  id = "data-table",
  columns = GenerateParentList(nSubjects),
  data = PkData2DashT(pkData),
  editable = TRUE,
  style_header = tableHeaderStyle,
  style_cell_conditional = styleConditional,
  active_cell = list(row = 0, column = 0, column_id = "time")
)

resultsTable <- dashDataTable(
  id = "results-table",
  style_header = tableHeaderStyle,
  style_cell_conditional = styleConditional
)

####################################################################################################

#################################### CREATE LAYOUT #################################################

app$layout(htmlDiv(
  className = "",
  children = list(
    htmlDiv(
      className = "pkcalc-banner",
      children = list(
        htmlA(
          id = "dashbio-logo",
          children = list(dashBioLogo),
          href = "/Portal"
          ),
        htmlH2("Noncompartmental Pharmacokinetics Analysis"),
        htmlA(
          id = "gh-link",
          children = list("View on GitHub"),
          href = paste("https://github.com/plotly/",
                  "dash-sample-apps/tree/master/apps/dashr-pk-calc", sep = ""),
          style = list(color = "white", border = "solid 1px white")
        ),
        githubLogo
      )
    ),
    htmlDiv(
      className = "container",
      children = list(
        htmlDiv(
          className = "row",
          style = list(),
          children = list(
            htmlDiv(
              className = "four columns pkcalc-settings",
              children = list(
                htmlP(list("Settings")),
                htmlDiv(list(
                  timePointsLabel,
                  subjectsLabel
                  )
                )
                )
              ),
            htmlDiv(
              className = "eight columns pkcalc-data-table",
              children = list(pkDataTable)
              )
            )
        ),
        htmlDiv(
          className = "row",
          children = list(
            htmlDiv(
              className = "six columns",
              children = list(dccGraph(id = "results-graph"))
            ),
            htmlDiv(
              className = "six columns pkcalc-results-table",
              children = list(resultsTable)
            )
          )
        )
        )
    )
  )
)
)

####################################################################################################

#################################### CALLBACKS START ###############################################

app$callback(output = list(id = "data-table", property = "columns"),
             params = list(
               input(id = "subjects-input", property = "value")
               ),

  function(subjects){
    return(GenerateParentList(subjects))
  }
)

app$callback(output = list(id = "data-table", property = "data"),
             params = list(
               input(id = "times-input", property = "value"),
               input(id = "subjects-input", property = "value"),
               state(id = "data-table", property = "data")
             ),

  function(rows, subjects, records){

   changeRow <- rows - length(records)

   if (changeRow > 0) {
     for (i in (length(records) + 1):(length(records) + changeRow)){
     #iterate through new rows

       records[[i]] <- records[[i - 1]]
       # Copy previous row content to new
       records[[i]][1:length(records[[i]])] <- ""
       # Replace with empty string
     }

   } else if (changeRow < 0) {
     records <- records[1:rows]
   }

   recordDf <- rbindlist(records) %>% select(time, everything())
   # Create df from records time first place

   changeCol <- subjects - (ncol(recordDf) - 1)

   if (changeCol > 0) {

     numLengths <- as.numeric(names(recordDf)[2:length(recordDf)])
     nameAppend <- numLengths[length(numLengths)] + 1
     # Get name for new column
     recordDf <- cbind(recordDf, "")
     # Create new column
     names(recordDf)[ncol(recordDf)] <- nameAppend
     # Set the name
   } else if (changeCol < 0) {
     recordDf <- select(recordDf, 1:(subjects + 1))
     # Remove columns
   }

   records <- df_to_list(recordDf)
   # Convert df back to list
   return(records)
  }
)

app$callback(output = list(id = "results-table", property = "columns"),
             params = list(
               input(id = "data-table", property = "data")
             ),

  function(results){
    return(GenerateResultColumn(results))
  }
)

app$callback(output = list(id = "results-table", property = "data"),
             params = list(
               input(id = "data-table", property = "data")
             ),

  function(records){
    results <- GenerateResultData(records)
  }
)

app$callback(output = list(id = "results-graph", property = "figure"),
             params = list(
               input(id = "data-table", property = "data")
             ),
  function(records){

    df <- rbindlist(records)
    # Convert records to df

    df <- as.data.frame(sapply(df, function(x) as.numeric(as.character(x))))
    # Convert anything non-numeric to NA

    df <- df[, colSums(is.na(df)) < nrow(df)]
    # Remove column if all NA

    dfMelt <- melt(df, id.vars = "time")
    names(dfMelt) <- c("x", "name", "y")
    dfMelt <- dfMelt[, c("x", "y", "name")]
    dfMelt <- mutate(dfMelt, mode = "lines+markers")
    # Reformat df to easily feed into list

    figData <- lapply(
      1:length(unique(dfMelt$name)),
      function(i){
        dfSub <- dfMelt[dfMelt$name == as.character(i - 1), ]
        list(
          x = as.numeric(dfSub$x),
          y = as.numeric(dfSub$y),
          name = paste("Subj", i, sep = ""),
          mode = "lines+markers"
        )
      }
    )
    # Filled the figData list

    figure <- list(
      data = figData,
      layout = list(
        xaxis = list(zeroline = FALSE),
        yaxis = list(
          title = list(
            text = "Conc (uM)",
            font = list(
              family = paste("'Open Sans', 'HelveticaNeue', 'Helvetica Neue',",
              "'Helvetica', 'Arial', 'sans-serif'", sep = ""),
              size = 12)
          ),
          type = "log",
          rangemode = "tozero",
          zeroline = FALSE,
          showticklabels = FALSE
        ),
        margin = list(l = 40, r = 30, b = 50, t = 50),
        showlegend = FALSE,
        height = 294,
        paper_bgcolor = "rgb(245, 247, 249)",
        plot_bgcolor = "rgb(245, 247, 249)"
      )
    )
  return(figure)
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv("PORT", 8050))
} else {
  app$run_server()
}
