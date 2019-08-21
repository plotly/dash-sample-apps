library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashDaq)
library(reticulate)
library(data.table)
library(rasterizer)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf("/app/apps/%s", appName))
}

################################################### helper function #############################################
pandas <- reticulate::import("pandas")
read_parquet <- function(path, columns = NULL) {
  
  path <- path.expand(path)
  path <- normalizePath(path)
  
  if (!is.null(columns)) columns = as.list(columns)
  
  xdf <- pandas$read_parquet(path, columns = columns)
  
  data.table::as.data.table(xdf, stringsAsFactors = FALSE)
}
easting_northing <- function(x) {
  
  longitude <- x[["longitude"]]
  latitude <- x[["latitude"]]
  
  if(is.null(longitude)) longitude <- x[[1]]
  if(is.null(latitude)) latitude <- x[[2]]
  
  origin_shift <- base::pi * 6378137
  easting <- longitude * origin_shift / 180.0
  northing <- log(tan((90 + latitude) * base::pi / 360.0)) * origin_shift / base::pi
  
  list(easting = easting, northing = northing)
}

################################################### load data #############################################
start_time <- Sys.time()
data <- lapply(0:35,
               function(i){
                 filename <- paste0("data/part.", i, ".parquet")
                 read_parquet(filename)
               })
# data <- data.table::rbindlist(data)
data <- do.call(rbind, data)
end_time <- Sys.time()
paste0("Time to load data: ", end_time - start_time)


################################################### some settings #############################################
USA           <- list(longitude = c(-124.72,  -66.95), latitude = c(23.55, 50.06))
LakeMichigan  <- list(longitude = c( -91.68,  -83.97), latitude = c(40.75, 44.08))
Chicago       <- list(longitude = c( -88.29,  -87.30), latitude = c(41.57, 42.00))
Chinatown     <- list(longitude = c( -87.67,  -87.63), latitude = c(41.84, 41.86))
NewYorkCity   <- list(longitude = c( -74.39,  -73.44), latitude = c(40.51, 40.91))
LosAngeles    <- list(longitude = c(-118.53, -117.81), latitude = c(33.63, 33.96))
Houston       <- list(longitude = c( -96.05,  -94.68), latitude = c(29.45, 30.11))
Austin        <- list(longitude = c( -97.91,  -97.52), latitude = c(30.17, 30.37))
NewOrleans    <- list(longitude = c( -90.37,  -89.89), latitude = c(29.82, 30.05))
Atlanta       <- list(longitude = c( -84.88,  -84.04), latitude = c(33.45, 33.84))

plot_width <- 600
plot_height <- 400
################################################### app #############################################
app <- Dash$new()

app$layout(
  htmlDiv(
    list(
      htmlH2("USA census summary"),
      htmlH5("The data has 300 million observations with size 7.36GB (It will take around 3.5 minutes to load the data set, please be patient)"),
      htmlLabel('Region'),
      dccDropdown(
        id = "region",
        options = list(list(label = "USA", value = "USA"),
                       list(label = "LakeMichigan", value = "LakeMichigan"),
                       list(label = "Chinatown", value = "Chinatown"),
                       list(label = "NewYorkCity", value = "NewYorkCity"),
                       list(label = "LosAngeles", value = "LosAngeles"),
                       list(label = "Houston", value = "Houston"),
                       list(label = "Austin", value = "Austin"),
                       list(label = "Atlanta", value = "Atlanta"),
                       list(label = "NewOrleans", value = "NewOrleans"),
                       list(label = "Chicago", value = "Chicago")),
        value = 'USA'
      ),
      htmlLabel('Reduction method'),
      dccDropdown(
        id = "reduc",
        options = list(list(label = "sum", value = "sum"),
                       list(label = "any", value = "any"),
                       list(label = "mean:on x", value = "mpx"),
                       list(label = "mean:on -x", value = "mnx"),
                       list(label = "mean:on y", value = "mpy"),
                       list(label = "mean:on -y", value = "mny")),
        value = 'sum'
      ),
      htmlLabel('Background'),
      dccDropdown(
        id = "background",
        options = list(list(label = "black", value = "black"),
                       list(label = "white", value = "white"),
                       list(label = "gray80", value = "gray80")
                       ),
        value = 'black'
      ),
      htmlLabel('Colour map'),
      dccDropdown(
        id = "cmap",
        options = list(list(label = "fire", value = "fire"),
                       list(label = "viridis", value = "viridis"),
                       list(label = "blue", value = "blue")),
        value = 'fire'
      ),
      htmlLabel('Scaling'),
      dccRadioItems(
        id = "scaling",
        options = list(list(label = "Log", value = "log"),
                       list(label = "Scale", value = "scale"),
                       list(label = "Original", value = "origin")),
        value = 'log'
      ),
      dccStore(id = "store"),
      dccGraph(
        id = 'datashade'
      )
    )
  )
)

app$callback(
  output = list(id = 'store', property = 'data'),
  params = list(input(id = 'region', property = 'value'), 
                input(id = 'datashade', property = 'relayoutData'),
                input(id = 'reduc', property = 'value')),
  function(region, relayout, reduc) {
   
    if(!is.null(unlist(relayout)) && !is.null(relayout$`xaxis.range[0]`)) {
      
      x_range <- c(relayout$`xaxis.range[0]`, relayout$`xaxis.range[1]`)
      y_range <- c(relayout$`yaxis.range[0]`, relayout$`yaxis.range[1]`)
      
    } else {
      
      range <- easting_northing(eval(parse(text = region)))
      x_range <- range$easting
      y_range <- range$northing
    }
    
    switch(reduc,
           "sum" = {
             mapping <- aes(x = easting, y = northing)
             reduction_func <- "sum"
           },
           "any" = {
             mapping <- aes(x = easting, y = northing)
             reduction_func <- "any"
           },
           "mpx" = {
             mapping <- aes(x = easting, y = northing, on = easting)
             reduction_func <- "mean"
           },
           "mnx" = {
             mapping <- aes(x = easting, y = northing, on = -easting)
             reduction_func <- "mean"
           },
           "mpy" = {
             mapping <- aes(x = easting, y = northing, on = northing)
             reduction_func <- "mean"
           },
           "mny" = {
             mapping <- aes(x = easting, y = northing, on = -northing)
             reduction_func <- "mean"
           },
    )
    
    data %>%
      canvas(mapping = mapping,
             plot_height = plot_height,
             plot_width = plot_width,
             x_range = x_range,
             y_range = y_range,
             show_raster = FALSE) %>%
      aggregation_points(reduction_func = reduction_func) %>%
      rasterizer() -> ds
    
    return(
      list(
        M = ds$agg[[1]][[1]],
        xlim = ds$lims[[1]]$xlim,
        ylim = ds$lims[[1]]$ylim
      )
    )
  }
)

app$callback(
  output = list(id = 'datashade', property = 'figure'),
  params = list(input(id = 'store', property = 'data'),
                input(id = 'cmap', property = 'value'),
                input(id = 'scaling', property = 'value'),
                input(id = 'reduc', property = 'value'),
                input(id = 'background', property = 'value')),
  function(data, cmap, scaling, reduc, background) {
    
    color <- if(cmap == "blue") {
      c("lightblue", "darkblue")
    } else {
      eval(parse(text = cmap))
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
    
    agg <- lapply(data[[1]], unlist) %>% 
      as.data.frame() %>% 
      as.matrix() %>%
      t()
    
    # scaling
    switch(
      scaling,
      "log" = {
        dimAgg <- dim(agg)
        agg <- if(reduc %in% c("mny", "mnx", "mpx", "mpy")) {
          not_background <- agg != 0
          if(all(agg[not_background] > 0)) {
            matrix(log(agg + 1), nrow = dimAgg[1])
          } else {
            min_agg <- min(agg)
            max_agg <- max(agg)
            agg[not_background] <- agg[not_background] + abs(min_agg) + abs(max_agg)
            matrix(log(agg + 1), nrow = dimAgg[1])
          }
        } else {
          matrix(log(agg + 1), nrow = dimAgg[1])
        }
      },
      "scale" = {
        dimAgg <- dim(agg)
        agg <- if(reduc %in% c("mny", "mnx", "mpx", "mpy")) {
          not_background <- agg != 0
          if(all(agg[not_background] > 0)) {
            min_agg <- min(agg)
            max_agg <- max(agg)
            matrix((agg - min(agg))/(max(agg) - min(agg)), nrow = dimAgg[1])
          } else {
            min_agg <- min(agg)
            max_agg <- max(agg)
            agg[not_background] <- agg[not_background] + abs(min_agg) + abs(max_agg)
            matrix((agg - min(agg))/(max(agg) - min(agg)), nrow = dimAgg[1])
          }
        } else {
          min_agg <- min(agg)
          max_agg <- max(agg)
          matrix((agg - min(agg))/(max(agg) - min(agg)), nrow = dimAgg[1])
        }

      },
      "origin" = {
        NULL
      }
    )
    
    l <- list(
      data =  list(
        list(
          x = seq(data[[2]][[1]], data[[2]][[2]], length.out = plot_width),
          y = seq(data[[3]][[1]], data[[3]][[2]], length.out = plot_height),
          z = agg,
          type = "heatmap",
          colorscale = colorscale
        )
      ),
      layout = list(
        xaxis = list('title' = 'longitude'),
        yaxis = list('title' = 'latitude'),
        margin = list('l' = 60, 'b' = 60, 't' = 10, 'r' = 10),
        hovermode = 'closest'
      )
    )
    return(l)
  }
)

app$run_server()