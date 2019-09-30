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

app <- Dash$new()

ridesRaw_1 <- "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data1.csv" %>%
  data.table::fread(stringsAsFactors = FALSE)
ridesRaw_2 <- "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data2.csv" %>%
  data.table::fread(stringsAsFactors = FALSE)
ridesRaw_3 <- "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data3.csv"  %>%
  data.table::fread(stringsAsFactors = FALSE)
ridesDf <- list(ridesRaw_1, ridesRaw_2, ridesRaw_3) %>%
  data.table::rbindlist()
test1 <- plot_ly(ridesDf, x = ~Lat, y = ~Lon, colorscale = 'Viridis', reversescale = FALSE) %>%
  add_rasterizer()

test<- layout(test1, font = list(color = 'rgb(226, 239, 250)'),
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


palette <- colorRampPalette(c("darkblue", "blue", "lightblue1",
                              "green","yellow", "red", "darkred"))


options <- htmlDiv(children =htmlDiv(list(
  htmlH2("Chart Options", style = list("font-size" = "23pt", "font-weight" = "200", "letter-spacing" = "1px")),
  htmlH4("Colorscale", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
  dccDropdown(
    id = "colorscale",
    value = "Viridis",
    options = list(
      list('label' = 'Viridis', 'value' = 'Viridis'),
      list('label' = 'Plasma', 'value' = 'Plasma'),
      list('label' = 'Blues', 'value' = 'Blues'),
      list('label' = 'Magma', 'value' = 'Magma'),
      list('label' = 'Greys', 'value' = 'Greys')
    )
  ),
  htmlBr(),
  htmlH4("Point Scaling", style = list("font-size" = "18pt", "font-weight" = "200", "letter-spacing" = "1px")),
  dccDropdown(
    id = 'scaling',
    value = 100
  )
)), className = 'item-a')

app$layout(
  header,
  htmlDiv(list(
    options,
    htmlDiv(children = dccGraph(id = 'rasterizer-output',
                                figure = test, style = list("height" = "88vh")), className = 'item-b'),
    dccStore(id = 'store')
  ), className = 'container'))

################################################### App Callbacks ###############################################
app$callback(
  output(id = 'rasterizer-output', property = 'figure'),
  params = list(
    input(id = 'store', property = 'data'),
    input(id = 'colorscale', property = 'value')
  ),
  update_graph <- function(data, colorscale) {
    
    if(colorscale == 'Blues' || colorscale == 'Plasma') {
      fix_scale = TRUE
    }
    
    else {
      fix_scale = FALSE
    }
    
    x_min <- data[[1]][[1]]
    x_max <- data[[1]][[2]]
    y_min <- data[[2]][[1]]
    y_max <- data[[2]][[2]]
    
    print(c(x_min, x_max, y_min, y_max))
    
    filtered_df_lat <- ridesDf[(ridesDf$Lat > x_min & ridesDf$Lat < x_max),]
    print(str(filtered_df_lat))
    filtered_df_lon <- filtered_df_lat[filtered_df_lat$Lon > y_min & filtered_df_lat$Lon < y_max,]
    print(str(filtered_df_lon))
    return(
      plot_ly(filtered_df_lon, x = ~Lat, y = ~Lon, colorscale = colorscale, reversescale = fix_scale) %>%
        add_rasterizer()
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