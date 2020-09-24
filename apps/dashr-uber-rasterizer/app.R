appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}



# Mapbox token for plotly | this one is for plot_mapbox figure
mapboxToken <- ("pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNrOWJqb2F4djBnMjEzbG50amg0dnJieG4ifQ.Zme1-Uzoi75IaFbieBDl3A")

# Setting mapbox token for R environment
Sys.setenv("MAPBOX_TOKEN" = mapboxToken)


library(plotly)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(data.table)
library(viridis)
library(rasterly)

app <- Dash$new()

ridesRaw_1 <- data.table::fread("https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data1.csv", stringsAsFactors = FALSE)
ridesRaw_2 <- data.table::fread("https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data2.csv", stringsAsFactors = FALSE)
ridesRaw_3 <- data.table::fread("https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data3.csv", stringsAsFactors = FALSE)
ridesDf <-  data.table::rbindlist(list(ridesRaw_1, ridesRaw_2, ridesRaw_3))

# subsetting on a matrix would likely be faster
ridesDf <- as.matrix(ridesDf[,-1])

default_colorscale <- lapply(0:256,
                             function(i) {
                               if(i == 0) {
                                 return(list(0, 'black'))
                               } else {
                                 return(list(i/256, viridis::plasma(256)[i]))
                               }
                             }
)

filtered_df_lat <- ridesDf[ridesDf[,"Lat"] > 39.9 & ridesDf[,"Lat"] < 42.1166, ]
filtered_df_lon <- ridesDf[ridesDf[,"Lon"] > -74.929 & ridesDf[,"Lon"] < -72.5, ]

initial_plot <- add_rasterly_heatmap(plot_ly(as.data.frame(filtered_df_lon), x = ~Lon, y = ~Lat,
                                             colorscale = default_colorscale,
                                             colorbar = list(title = "Log(No. of Rides)")))

default_plot <- layout(initial_plot, font = list(color = 'rgb(226, 239, 250)'),
                       margin = list(l = 1, r = 1, b = 1, t = 15, pad = 0, autoexpand = TRUE),
                       paper_bgcolor='rgb(38, 43, 61)',
                       plot_bgcolor='rgb(38, 43, 61)',
                       xaxis = list(title = "Longitude",
                                    constrain = "domain",
                                    automargin = TRUE,
                                    scaleanchor = "y",
                                    scaleratio = cos(40.8*pi/180)),
                       yaxis = list(title = "Latitude",
                                    automargin = TRUE,
                                    constrain = "domain"))


################################################### App Layout #####################################
header <- htmlDiv(
  id = "app-page-header",
  style = list(
    width = "100%",
    background = "#262B3D" ,
    color = "#E2EFFA"
  ),
  children = list(
    htmlA(
      id = "dash-logo",
      children = list(
        htmlImg(src='assets/plotly-dash-logo.png', height = '36', width = '180',
                style = list('top' = '10', 'margin' = '10px'))
      ),
      href = "/Portal"
    ),
    htmlH2("Uber NYC Rasterizer"),
    htmlA(
      id = "gh-link",
      children = list("View on GitHub"),
      href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-uber-rasterizer",
      style = list(color = "white", border = "solid 1px white")
    ),
    htmlImg(
      src = "assets/GitHub-Mark-Light-64px.png"
    )
  )
)


tabs <- htmlDiv(dccTabs(id = 'circos-control-tabs', value = 'what-is', children = list(
  dccTab(
    label = 'About',
    value = 'what-is',
    children = htmlDiv(
      id = 'control-tab', children = list(
        htmlH4("What is Uber NYC Rasterizer?", style = list("font-size" = "24pt", "font-weight" = "200", "letter-spacing" = "1px")),
        dccMarkdown("This Dash app is a simple demonstration of the rasterizing capabilities of the _rasterly_ package.
              The dataset consists of over 4.5 million observations, representing Uber rides taken in New York City in 2014.
              In CSV format, the source data are over 165 MB in size. _rasterly_ is capable of processing datasets an order
              of magnitude larger in similarly brisk fashion. The raster data required to produce the aggregation layers and color
              gradients displayed here are computed efficiently enough to maintain the interactive feel of the application.
              ", style = list("padding" =  "5px")),
        dccMarkdown("The image will eventually transition to a Mapbox map view if you zoom in far enough; to return to the rasterized
              view, just click the Reset Graph button, which is accessible by clicking the Options tab.
                    ", style = list("padding" = "5px")),
        dccMarkdown("Visit the _rasterly_ package repository [here](https://github.com/plotly/rasterly) to learn more.
                    ", style = list("padding" = "5px"))
      )
    )
  ),

  dccTab(
    label = 'Options',
    value = 'data',
    children = htmlDiv(className = 'circos-tab', children = list(
      htmlDiv(className = 'app-controls-block', children= list(
        htmlH4("Colorscale", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
        htmlDiv(dccDropdown(
          id = "cmap",
          value = "plasma",
          options = list(
            list('label' = 'Plasma', 'value' = 'plasma'),
            list('label' = 'Viridis', 'value' = 'viridis'),
            list('label' = 'Blues', 'value' = 'blue'),
            list('label' = 'Magma', 'value' = 'fire')
          )
        ), style = list("color" = "white")),
        htmlH4("Background", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
        htmlDiv(dccDropdown(
          id = "background",
          value = "black",
          options = list(
            list('label' = 'Black', 'value' = 'black'),
            list('label' = 'Grey', 'value' = 'grey'),
            list('label' = 'White', 'value' = 'white')
          )
        ), style = list("color" = "white")),
        htmlH4("Point Scaling", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
        dccDropdown(
          id = 'scaling',
          value = 'log',
          options = list(
            list(label = 'Log', value = 'log'),
            list(label = 'Origin', value = 'origin')
          )
        ),
        htmlH4("Reduction Method", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
        dccDropdown(
          id = "reduc",
          options = list(list(label = "sum", value = "sum"),
                         list(label = "any", value = "any"),
                         list(label = "mean", value = "mean")),

          value = 'sum'
        ),
        htmlH4("Pixel Size", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
        dccSlider(
          id = 'point-size',
          min = 0,
          max = 10,
          step = 1,
          value = 0,
          marks <- as.list(
            setNames(
              as.character(seq(1:10)),
              as.character(seq(1:10))
            )
          )
        ),
        htmlBr(),
        htmlBr(),
        htmlButton(
          id = 'reset-button',
          n_clicks = 0,
          children = 'Reset Graph'
        )
      ))
    ))
  )
)))


options <- htmlDiv(list(
  tabs
), className = 'item-a')


app$layout(htmlDiv(list(
  header,
  htmlDiv(list(
    options,
    htmlDiv(children = dccGraph(id = 'rasterizer-output',
                                figure = default_plot, style = list("height" = "88vh")), className = 'item-b'),
    dccStore(id = 'store')
  ), className = 'container')
)))

################################################### App Callbacks ###############################################

# Callback to reset the reset button.

app$callback(
  output(id = 'reset-button', property = 'n_clicks'),
  params = list(
    input(id = 'rasterizer-output', property = 'relayoutData')
  ),
  reset_clicks <- function(relayout) {
    if (is.null(relayout[1]) == FALSE) {
      return(0)
    }
  }
)

# Callback to subset and return data ranges based on the zoom of the plot.

app$callback(
  output(id = 'store', property = 'data'),
  params = list(
    input(id ='rasterizer-output', property = 'relayoutData'),
    input(id ='reset-button', property = 'n_clicks')
  ),
  update_stored_data <- function(relayout, n_clicks) {
    if (n_clicks > 0) {
      x_range <- c(min(ridesDf[,"Lon"], na.rm=T), -72.5)
      y_range <- c(39.9, max(ridesDf[,"Lat"], na.rm=T))
      differences = c(10,10)
      return(list(x_range, y_range, differences))
    }
    else {
      if (length(relayout) == 4) {
        x_range <- c(relayout$`xaxis.range[0]`, relayout$`xaxis.range[1]`)
        y_range <- c(relayout$`yaxis.range[0]`, relayout$`yaxis.range[1]`)

        x_difference <- relayout$`xaxis.range[1]` - relayout$`xaxis.range[0]`
        y_difference <- relayout$`yaxis.range[1]` - relayout$`yaxis.range[0]`

        differences = c(x_difference, y_difference)
      }

      else {
        x_range <- c(min(ridesDf[,"Lon"], na.rm=T), -72.5)
        y_range <- c(39.9, max(ridesDf[,"Lat"], na.rm=T))
        differences = c(10,10)
      }
    }
    return(list(x_range, y_range, differences))
  }
)

# Callback to generate rasterization plot.

app$callback(
  output(id = 'rasterizer-output', property = 'figure'),
  params = list(
    input(id = 'store', property = 'data'),
    input(id = 'cmap', property = 'value'),
    input(id = 'background', property = 'value'),
    input(id = 'reduc', property = 'value'),
    input(id = 'scaling', property = 'value'),
    input(id = 'point-size', property = 'value')
  ),
  update_graph <- function(data, cmap, background, reduc, scale, point_size) {
    color <- if(cmap == "blue") {
      c("lightblue", "darkblue")
    } else if(cmap =="viridis") {
      color <- viridis(256)
    } else if (cmap =="plasma"){
      color <- plasma(256)
    } else {
      eval(parse(text = cmap))
    }

    if(background != "black") {
      color <- rev(color)
    }

    len_col <- length(color)

    colorscale <- lapply(0:len_col,
                         function(i) {
                           if(i == 0) {
                             return(list(0, background))
                           } else {
                             return(list(i/len_col, color[i]))
                           }
                         }
    )

    x_min <- data[[1]][[1]]
    x_max <- data[[1]][[2]]
    y_min <- data[[2]][[1]]
    y_max <- data[[2]][[2]]


    x_difference = data[[3]][[1]]
    y_difference = data[[3]][[2]]


    filtered_df_lat <- ridesDf[ridesDf[, "Lat"] > y_min & ridesDf[, "Lat"] < y_max, ]
    filtered_df_lon <- filtered_df_lat[filtered_df_lat[,"Lon"] > x_min & filtered_df_lat[,"Lon"] < x_max, ]

    colorbar_title <- ifelse(scale == "log", "Log(No. of Rides)", "No. of Rides")
    # plot_ly requires a data.frame

    if (x_difference < 0.1004274 || y_difference < 0.0409278) {
      return(plot_mapbox(as.data.frame(ridesDf)[sample(nrow(as.data.frame(ridesDf)), 10000, replace = TRUE),], lon = ~Lon, lat = ~Lat) %>%
               layout(mapbox = list(zoom = 12,
                                    accesstoken = mapboxToken,
                                    center = list(lat = mean(as.data.frame(filtered_df_lon)$Lat),
                                                  lon = mean(as.data.frame(filtered_df_lon)$Lon)),
                                    style = "dark"
               ),
               font = list(color = 'rgb(226, 239, 250)'),
               margin = list(l = 5, r = 5, b = 5, t = 15, pad = 0, autoexpand = TRUE),
               paper_bgcolor='rgb(38, 43, 61)',
               xaxis = list(title = "Longitude"),
               yaxis = list(title = "Latitude")
               )
      )
    }
    else{
      return(
        plot_ly(as.data.frame(filtered_df_lon), x = ~Lon, y = ~Lat,
                colorscale = colorscale,
                colorbar = list(title = colorbar_title)) %>%
          add_rasterly_heatmap(reduction_func = reduc, scaling = scale, size = b <- if(point_size == 0) NULL else {point_size}) %>%
          layout(font = list(color = 'rgb(226, 239, 250)'),
                 margin = list(l = 1, r = 1, b = 1, t = 15, pad = 0, autoexpand = TRUE),
                 paper_bgcolor='rgb(38, 43, 61)',
                 plot_bgcolor='rgb(38, 43, 61)',
                 xaxis = list(title = "Longitude",
                              constrain = "domain",
                              scaleanchor = "y",
                              scaleratio = cos(40.8*pi/180)),
                 yaxis = list(title = "Latitude",
                              constrain = "domain"))
      )
    }

  }
)


if(appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(host = "127.0.0.1",
                 port=8050,
                 debug = FALSE)
}
