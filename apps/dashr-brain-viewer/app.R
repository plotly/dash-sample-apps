library(dashColorscales)
library(dash)
library(plotly)
library(rapportools)
library(dashCoreComponents)
library(dashHtmlComponents)
library(rjson)


appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("/app/apps/%s", appName))
}


# Load Data ------------------------------

# The following data are created using source('brainDATA.R')
load('data/alldata.RData', verbose = TRUE)

data <- initialGraphdata


# Helper Variables / Functions  ------------------------------

standard_intensity <- function(x, y, z) {
  return(z)
}

default_colorscale_no_index <- list("#0c3383", "#0a88ba", "#f2d338", "#f28f38", "#d91e1e")

axis_template <- list(
  showbackground = TRUE,
  backgroundcolor = '#0a0a0a',
  gridcolor = '#ffffff',
  zerolinecolor = '#ffffff'
)

plot_layout <- list(
  title = '',
  margin = list(t = 0, b = 0, l = 0,  r = 0  ),
  font = list(size = 12, color = 'white'),
  autoresize = TRUE,
  showlegend = FALSE,
  plot_bgcolor = 'black',
  paper_bgcolor = 'black',
  scene = list(
    xaxis = axis_template,
    yaxis = axis_template,
    zaxis = axis_template,
    aspectratio = list(x = 1, y = 1.2, z = 1),
    camera = list(eye = list( x = 1.25, y = 1.25, z = 1.25)),
    annotations = list()
  )
)

styles <- list(
  'pre' = list(
    'border' = 'thin lightgrey solid',
    'padding' = '10px',
    'marginBottom' = '20px'
  ),
  'graph' = list('userSelect' = 'none',
                 'margin' = 'auto'))


# Create Layout ------------------------------

app = Dash$new()

app$layout(htmlDiv(children = list(

  htmlDiv(list(
    htmlImg(src = "assets/dash-logo.png", "height" = "50px"),
    htmlP(
      'MRI Reconstruction',
      style = list(
        color = '#506784',
        fontSize = '29px',
        fontFamily = 'Script'
      )
    )
  ), style = list("float" = "top", 'display' = "flex")),
  htmlP(
    children = list(
      htmlP(
        'Click on the brain to add an annotation.
        Drag the black corners of the graph to rotate.',
        style = list(fontFamily = 'Trebuchet MS')
      ),
      htmlBr(),
      htmlA(
        children = 'Learn More',
        href = 'https://github.com/plotly/dash-sample-apps',
        style = list(
          'color' = 'white',
          fontSize = '13px',
          'text-decoration' = 'None',
          fontFamily = 'Trebuchet MS',
          border = 'solid',
          'border-width' = '0.1px',
          padding = '15px',
          'margin-right' = '30px'
        )
      ),
      htmlA(
        children = 'View On Github',
        href = 'https://github.com/plotly/dash-brain-surface-viewer',
        style = list(
          'color' = 'white',
          fontSize = '13px',
          'text-decoration' = 'None',
          fontFamily = 'Trebuchet MS',
          border = 'solid',
          'border-width' = '0.1px',
          padding = '15px'
        )
      )
    )
  ),
  htmlBr(),
  htmlDiv(list(
    htmlP(
      children = 'Click colorscale to change:',
      style = list('display' = 'inline-block', 'fontSize' = '12px')
    ),
    htmlDiv(list(
      dashColorscales(id = 'colorscale-picker',
                      colorscale = default_colorscale_no_index)
    ),
    style = list('marginTop' = '-15px', 'marginLeft' = '-30px')),
    htmlDiv(list(
      dccRadioItems(
        options = list(
          list('label' = 'Cortical', 'value' = 'human'),
          list('label' = 'Thickness Mouse Brain', 'value' = 'mouse'),
          list('label' = 'Brain Atlas', 'value' = 'human_atlas')
        ),
        value = 'human_atlas',
        id = 'radio-options',
        labelStyle = list('display' = 'inline-block')
      )
    )),
    htmlBr(),
    dccMarkdown(children = '**Click Data** | Click on points in the graph.'),
    htmlDiv(
      list(
        htmlPre(
          id = 'click-data',
          style = list(
            'border' = 'None',
            'padding' = '10px',
            'marginBottom' = '20%',
            'resize' = 'both',
            'overflow' = 'auto'
          )
      )),
      style = list('background-color' = '#282828')),

    dccMarkdown(children = '**Relayout Data** | Drag the graph corners to rotate it.'),

    htmlDiv(
      list(
        htmlPre(
        id = 'relayout-data',
        style = list(
          'border' = 'None',
          'padding' = '10px',
          'marginBottom' = '20px',
          'resize' = 'both',
          'overflow' = 'auto')
        )
      ),
      style = list('background-color' = '#282828')),
    htmlP(children = list(
      'Dash/R code on ',
      htmlA(
        children = 'GitHub',
        target = '_blank',
        href = 'https://github.com/plotly/dash-brain-surface-viewer',
        style = list('color' = '#F012BE')
      ),
      '. Brain data from Mcgill\'s ACE Lab',
      htmlA(
        children = 'Surface Viewer',
        target = '_blank',
        href = 'https://brainbrowser.cbrain.mcgill.ca/surface-viewer#ct',
        style = list('color' = '#F012BE')
      )
    ))),
    style = list(
      'float' = 'right',
      'background-color' = "#1E1E1E",
      'width' = '29%',
      'height' = '1300px',
      padding = '10px'
    )),

  htmlDiv(
    dccGraph(
      id = 'brain-graph',
      figure = list(data = list(data),
                    layout = plot_layout),
      config = list('editable' = TRUE, 'scrollZoom' = FALSE),
      style = list(
        'userSelect' = 'none',
        'width' = '100%',
        'height' = '900px',
        "flex" = "5",
        "position" = "relative"
      )
    ),
    style = list(
      'float' = 'left',
      'width' = '60%',
      'height' = '50%',
      padding = 'None'
    )
  )
)))


# Callbacks Start ------------------------------

app$callback(output = list(id = 'brain-graph', property = 'figure'),
             params = list(input(id = 'brain-graph', property = 'clickData'),
                           input(id = 'radio-options', property = 'value'),
                           input(id = 'colorscale-picker', 'colorscale'),
                           state(id = 'brain-graph', property = 'figure')),
  function(clickData, val, colorrs, figure) {

    if (data[['name']] != val) {
      if (val == 'human') {
        traces = realct
        figure = list(data = list(traces),
                      layout = plot_layout)
      } else if (val == 'human_atlas') {
        traces = surfreg
        figure = list(data = list(traces),
                      layout = plot_layout)
      } else if (val == 'mouse') {
        traces = mouse_map
        figure = list(data = traces,
                      layout = plot_layout)
      }
    }

    if (!is.null(clickData[[1]])) {
      marker = list(
        x = list(clickData[['points']][[1]][['x']]),
        y = list(clickData[['points']][[1]][['y']]),
        z = list(clickData[['points']][[1]][['z']]),
        mode = 'markers',
        marker = list(size = 15, line = list(width = 3)),
        name = 'Marker',
        type = 'scatter3d',
        text = 'Click point to remove annotation'
      )
      anno = list(
        list(
          x = clickData[['points']][[1]][['x']],
          y = clickData[['points']][[1]][['y']],
          z = clickData[['points']][[1]][['z']],
          font = list(color = 'black'),
          bgcolor = 'white',
          borderpad = 5,
          bordercolor = 'black',
          borderwidth = 1,
          captureevents = TRUE,
          ay = -50,
          arrowcolor = 'white',
          arrowwidth = 2,
          arrowhead = 0,
          text = 'Click here to annotate<br>(Click point to remove)'
        )
      )
      if (length(figure[['data']]) > 0) {
        same_point_found <- FALSE
        for (i in 1:length(figure[['data']])) {
          for (point in figure['data']) {
            if (unlist(point[[i]]['x']) == marker[['x']][[1]] &
                unlist(point[[i]]['y']) == marker[['y']][[1]] &
                unlist(point[[i]]['z']) == marker[['z']][[1]]) {
              anno_trace_index_offset <- 1
              if (val == 'mouse') {
                anno_trace_index_offset <- 2
              }
              figure[['data']][i] <- NULL
              #print(list('DEL. MARKER', i, figure[['layout']][['scene']][['annotations']]))
              if (length(figure[['layout']][['scene']][['annotations']]) >= (i -
                                                                             anno_trace_index_offset)) {
                figure[['layout']][['scene']][['annotations']][i - anno_trace_index_offset] <- NULL
              }
              same_point_found <- TRUE
            }
          }
        }
        if (same_point_found == FALSE) {
          figure[['data']] <- append(figure[['data']], list(marker))
          figure[['layout']][['scene']][['annotations']] <-
            append(figure[['layout']][['scene']][['annotations']], anno)
        }
      } else {
        figure[['data']] <- append(figure[['data']], list(marker))
        figure[['layout']][['scene']][['annotations']] <-
          append(figure[['layout']][['scene']][['annotations']], anno)
      }
    }

    cs <- list()
    for (i in seq_along(colorrs)) {
      cs[[i]] <- list(i / (length(colorrs) - 1) - .20, colorrs[[i]])
    }
    figure[["data"]][[1]][['colorscale']] <- cs
    return(figure)
  }
)


app$callback(output = list(id = 'click-data', property = 'children'),
             params = list(input(id = 'brain-graph', property = 'clickData')),
  function(clickData) {
    return(toJSON(clickData, indent = 4))
  }
)

app$callback(output = list(id = 'relayout-data', property = 'children'),
             params = list(input(id = 'brain-graph', property = 'relayoutData')),
  function(relayoutData){
    return(toJSON(relayoutData, indent = 4))
  }
)


if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv("PORT", 8050))
} else {
  app$run_server()
}
