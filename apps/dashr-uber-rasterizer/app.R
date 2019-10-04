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
library(data.table)
library(rasterizer)
library(viridis)


app <- Dash$new()

ridesRaw_1 <- "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data1.csv" %>%
  data.table::fread(stringsAsFactors = FALSE)
ridesRaw_2 <- "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data2.csv" %>%
  data.table::fread(stringsAsFactors = FALSE)
ridesRaw_3 <- "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data3.csv"  %>%
  data.table::fread(stringsAsFactors = FALSE)
ridesDf <- list(ridesRaw_1, ridesRaw_2, ridesRaw_3) %>%
  data.table::rbindlist()


default_colorscale <- lapply(0:256,
                             function(i) {
                               if(i == 0) {
                                 return(list(0, 'black'))
                               } else {
                                 return(list(i/256, viridis::viridis(256)[i]))
                               }
                             }
)

filtered_df_lat <- ridesDf %>% filter(Lat > 39.6569 & Lat < 42.1166)
filtered_df_lon <- filtered_df_lat %>% filter(Lon > -74.929 & Lon < -72.066)

initial_plot <- plot_ly(filtered_df_lon, x = ~Lat, y = ~Lon, colorscale = default_colorscale) %>%
  add_rasterizer()

default_plot <- layout(initial_plot, font = list(color = 'rgb(226, 239, 250)'),
                       paper_bgcolor='rgb(38, 43, 61)',
                       plot_bgcolor='rgb(38, 43, 61)')


################################################### App Layout ##################################################
header <- htmlDiv(
  id = "app-page-header",
  style = list(
    width = "100%",
    background = "#262B3D" ,
    color = "#E2EFFA"
  ),
  children = list(
    htmlA(
      id = "dashbio-logo",
      children = list(
        htmlImg(src='assets/plotly-dash-bio-logo.png', height = '36', width = '180',
                style = list('top' = '10', 'margin' = '10px'))
      ),
      href = "/Portal"
    ),
    htmlH2("Uber NYC Rasterizer"),
    htmlA(
      id = "gh-link",
      children = list("View on GitHub"),
      href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-clustergram",
      style = list(color = "white", border = "solid 1px white")
    ),
    htmlImg(
      src = "assets/GitHub-Mark-Light-64px.png"
    )
  )
)


options <- htmlDiv(children =htmlDiv(list(
  htmlH2("Chart Options", style = list("font-size" = "23pt", "font-weight" = "200", "letter-spacing" = "1px")),
  htmlH4("Colorscale", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
  htmlDiv(dccDropdown(
    id = "cmap",
    value = "viridis",
    options = list(
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
  htmlBr(),
  htmlH4("Point Scaling", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
  dccDropdown(
    id = 'scaling',
    value = 100
  ),
  htmlH4("Reduction method", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
  dccDropdown(
    id = "reduc",
    options = list(list(label = "sum", value = "sum"),
                   list(label = "any", value = "any"),
                   list(label = "mean", value = "mean")),

    value = 'sum'
  )
)), className = 'item-a')

app$layout(
  header,
  htmlDiv(list(
    options,
    htmlDiv(children = dccGraph(id = 'rasterizer-output',
                                figure = default_plot, style = list("height" = "88vh")), className = 'item-b'),
    dccStore(id = 'store')
  ), className = 'container'))

################################################### App Callbacks ###############################################
app$callback(
  output(id = 'rasterizer-output', property = 'figure'),
  params = list(
    input(id = 'store', property = 'data'),
    input(id = 'cmap', property = 'value'),
    input(id = 'background', property = 'value'),
    input(id = 'reduc', property = 'value')
  ),
  update_graph <- function(data, cmap, background, reduc) {

    color <- if(cmap == "blue") {
      c("lightblue", "darkblue")
    } else if(cmap =="viridis") {
      color <- viridis(256)
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
    
    
    filtered_df_lat <- ridesDf[ridesDf$Lat > x_min & ridesDf$Lat < x_max,]
    filtered_df_lon <- filtered_df_lat[filtered_df_lat$Lon > -y_min & filtered_df_lat$Lon < y_max,]
    
    return(
      plot_ly(filtered_df_lon, x = ~Lat, y = ~Lon, colorscale = colorscale) %>%
        add_rasterizer(reduction_func = reduc)
      %>%
        layout(font = list(color = 'rgb(226, 239, 250)'),
               paper_bgcolor='rgb(38, 43, 61)',
               plot_bgcolor='rgb(38, 43, 61)')
    )
  }
)


app$callback(
  output(id = 'store', property = 'data'),
  params = list(
    input('rasterizer-output', property = 'relayoutData')
  ),
  update_stored_data <- function(relayout) {
    if (length(relayout) == 4) {
      x_range <- c(relayout$`xaxis.range[0]`, relayout$`xaxis.range[1]`)
      y_range <- c(relayout$`yaxis.range[0]`, relayout$`yaxis.range[1]`)

    }

    else {
      x_range <- c(min(ridesDf$Lat), max(ridesDf$Lat))
      y_range <- c(min(ridesDf$Lon), max(ridesDf$Lon))
    }

    return(list(x_range, y_range))
  }
)

if(appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(showcase = TRUE)
}

