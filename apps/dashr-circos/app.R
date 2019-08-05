library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(jsonlite)
library(dashBio)
library(GEOquery)
library(dashTable)


appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf("/app/apps/%s", appName))
}

app <- Dash$new()


circos_dataframe <- jsonlite::fromJSON("assets/circos_graph_data.json")


circos_graph_data <- read_json(path = json_file)

json_file <- "assets/circos_graph_data.json"


source("assets/CircosFunctions.R")

# Load up the preloaded dataset into R structures.



applayout <- circos_layout()


header_colors <- function() {
  return(list('bg_color' = '#262B3D', 'font_color' = '#FFF', 'light_logo' = TRUE))
}
  
  
  
  
 

header <- htmlDiv(
  id = "app-page-header",
  style = list(
    width = "100%",
    background = header_colors()[["bg_color"]],
    color = header_colors()[["font_color"]]
  ),
  children = list(
    htmlA(
      id = "dashbio-logo",
      children = list(
        htmlImg(src='assets/plotly-dash-bio-logo.png', height = '36', width = '190',
                style = list('top' = '10', 'margin-left' = '10px'))
      ),
      href = "/Portal"
    ),
    htmlH2("Circos"),
    htmlA(
      id = "gh-link",
      children = list("View on GitHub"),
      href = "https://github.com/plotly/dash-sample-apps/blob/master/apps/dashr-circos/app.R",
      style = list(color = "white", border = "solid 1px white")
    ),
    htmlImg(
      src = "assets/GitHub-Mark-Light-64px.png"
    )
  )
)



app$layout(htmlDiv(list(
  header,
  htmlBr(),
  htmlBr(),
  applayout
)))




# Callback to update the text for chords graph.

app$callback(
  output(id = "chords-text", property = "children"),
  params = list(
    input(id = 'circos-graph-type', property = 'value')
  ),
  update_desc <- function(option) {
    if (option == 'chords') {
      return(htmlP("Highlight chords on the graph by selecting the 
                   'Table' tab and selecting rows in the 'Chords' dataset."))
    }
    
    else {
      return("    ")
    }
  }
)


# Callback to update the table based on the selected dataset.

app$callback(
  output(id = 'circos-table-container', property = 'children'),
  params = list(
    input(id = 'circos-view-dataset', property = 'value')
  ),

  update_table <- function(selection) {
    df <- circos_dataframe[[selection]]

    df <- flatten(df, recursive = TRUE)
    
    table <- dashDataTable(
      id = 'data-table',
      columns = lapply(colnames(df), 
                       function(colName){
                         list(
                           id = colName,
                           name = colName
                         )
                       }),
      data= df_to_list(df),
      row_selectable='multi',
      selected_rows = character(0),
      # sorting = TRUE,
      # filtering = TRUE,
      # css = list(list(
      #   'selector' = '.dash-cell div.dash-cell-value',
      #   'rule' = 'display: inline;',
      #            'white-space: inherit;',
      #            'overflow: auto;',
      #            'text-overflow: inherit;'
      # )),
      style_cell = list(
        'whiteSpace' = 'no-wrap',
        'overflow' = 'hidden',
        'textOverflow' = 'ellipsis',
        'maxWidth' = 100,
        'fontWeight' = 100,
        'fontSize' = '11pt',
        'fontFamily' = 'Courier New',
        'backgroundColor' = '#1F2132'
      ),
      
      style_header = list(
        'backgroundColor' = '#1F2132',
        'textAlign' = 'center'
      ),
      
      style_table = list(
        'maxHeight' = '310px',
        'width' = '320px',
        'marginTop' = '5px',
        'marginBottom' = '10px'
      )
    )
    

    return(table)
  }
)


# Callback to update the displayed graph based on the selection.

app$callback(output = list( id= 'circos-hold',property = 'children'),
             params = list(input(id = 'circos-graph-type', property = 'value'),
                           input(id = 'circos-size', property = 'value')
             ),
             
             update_graph <- function(key, size, rows) {
               return(get_circos_graph(key = key, size = size))
             }
)


# Callback to show and hide graphs based on the table rows selected or not selected.

app$callback(output = list( id= 'graph-container',property = 'style'),
             params = list(input(id = 'data-table', property = 'selected_rows')
             ),
             hide_graph <- function(selected_rows) {
               if (length(selected_rows) > 0) {
                 style = list('display' = 'none')
                 return(style)
               }
             }
)


 app$callback(output = list( id= 'chords-plot',property = 'style'),
              params = list(input(id = 'data-table', property = 'selected_rows')
              ),
              hide_graph <- function(selected_rows) {
                if (length(selected_rows) <= 0) {
                  style = list('display' = 'none')
                  return(style)
                }
              }
 )



app$callback(output = list( id= 'options', property = 'style'),
             params = list(input(id = 'circos-preloaded-uploaded', property = 'value')
             ),
             show_options <- function(options) {
               if (options == 'preloaded') {
                 style = list('display' = 'none')
                 return(style)
               }
             }
)


app$callback(
  output(id = 'chords-plot', property = 'children'),
  params = list(
    input(id = 'data-table' , property = 'selected_rows'),
    input(id = 'circos-size', property = 'value')
    
  ),

  update_rows <- function(rows, size) {
    if (length(rows) > 0) {
      for (i in 0:length(circos_graph_data$chords)) {
        if (i %in% rows) {
          circos_graph_data$chords[[i + 1]]$color <- '#00cc96'
        }
      }
      return(dashbioCircos(
        id = 'random-circos',
        selectEvent = list('0' = 'both', '1' = 'both'),
        layout = circos_graph_data[['GRCh37']],
        config = list(
          'innerRadius' = size/2 - 80,
          'outerRadius' = size/2 - 30,
          'ticks' = list('display' = FALSE, 'labelDenominator' = 1000000),
          'labels' = list(
            'position' = 'center',
            'display' = TRUE,
            'size' = 11,
            'color' = '#fff',
            'radialOffset' = 75
          )
        ),
  
        tracks = list(
          list(
            'type' = 'HIGHLIGHT',
            'data' = circos_graph_data[['cytobands']],
            'config' = list(
              'innerRadius' = size/2 - 80,
              'outerRadius' = size/2 - 40,
              'opacity' = 0.3,
              'tooltipContent' = list('name' = 'all'),
              'color' = list('name' = 'color')
            )
          ),
          list(
            'type' = 'CHORDS',
            'data' = circos_graph_data[['chords']],
            'config' = list(
              'opacity' = 0.7,
              'color' = list('name' = 'color'),
              'tooltipContent' = list(
                'source' = 'source',
                'sourceID' = 'id',
                'target' = 'target',
                'targetID' = 'id',
                'targetEnd' = 'end'
              )
            )
          )
        ), size = 700
      ))
    }
  }
)




app$callback(
  output(id = 'event-data-select', property = 'children'),
  params = list(
    input(id = 'main-circos', property = 'eventDatum')
  ),
  
  event_data_select <- function(s) {
    contents = list('There is no event data. Hover over the circos graph to access event data.')
    
    output <- lapply(
      1:length(s),
      function (i) {
        return(htmlDiv(list(
          htmlSpan(names(s)[[i]], style = list('fontWeight' = "bold", 'textTransform' = 'capitalize')),
          " - ", 
          s[[i]]
        )))
      }
    )
    
    return(output)
  }
)



if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server(showcase = TRUE)
}


