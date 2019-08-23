library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
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
  if (!is.null(columns)) columns <- as.list(columns)
  path.expand(path) %>% 
    normalizePath() %>%
    pandas$read_parquet(., columns = columns) %>%
    data.table::as.data.table(., stringsAsFactors = FALSE)
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
data <- data.table::rbindlist(data)
end_time <- Sys.time()
paste0("Time to load data: ", end_time - start_time)

################################################### race #############################################
race_all_region <- data[, .N, by = race]
race_full_names <- c("White", "Hispanics", "Black", "Asian", "Other")

race_table <- function(names, table) {
  if(is.null(table$N)) stop("No counting number")
  if(is.null(table$race)) stop("No race")
  sapply(names, 
         function(name) {
           switch(name,
                  "Hispanics" = {
                    x <- table$N[table$race == "h"]
                    if(length(x) == 0) 0 else x
                  },
                  "White" = {
                    x <- table$N[table$race == "w"]
                    if(length(x) == 0) 0 else x
                  },
                  "Black" = {
                    x <- table$N[table$race == "b"]
                    if(length(x) == 0) 0 else x
                  },
                  "Asian" = {
                    x <- table$N[table$race == "a"]
                    if(length(x) == 0) 0 else x
                  },
                  "Other" = {
                    x <- table$N[table$race == "o"]
                    if(length(x) == 0) 0 else x
                  },
           ) 
         })
}
race_all <- race_table(names = race_full_names, 
                       table = race_all_region)
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
bg <- "#F8F8F8"
################################################### app #############################################
app <- Dash$new()

tab1 <- dccTab(
  label = "Generate Graph",
  children = list(
    htmlDiv(
      className = "control",
      style = list(
        paddingBottom = "0px"
      ),
      children = list(
        htmlH5("Generate a rasterizer object"),
        htmlP(
          "Note that the data has 300 million observations with a file size of 7.36GB and 
          each generation of the graph may take around 5 seconds to update. 
          The initial generation may take more time",
          style = list(color = "gray")
        )
      )
    ),
    htmlDiv(
      className = "control",
      children = list(
        htmlH6("Region"),
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
        )
      )
    ),
    htmlDiv(
      className = "control",
      children = list(
        htmlH6("Reduction method"),
        dccDropdown(
          id = "reduc",
          options = list(list(label = "sum", value = "sum"),
                         list(label = "any", value = "any"),
                         list(label = "mean:on x", value = "mpx"),
                         list(label = "mean:on -x", value = "mnx"),
                         list(label = "mean:on y", value = "mpy"),
                         list(label = "mean:on -y", value = "mny")),
          value = 'sum'
        )
      )
    )
  )
)

tab2 <- dccTab(
  label = "Graph setting",
  children = list(
    htmlDiv(
      className = "control",
      children = list(
        htmlH5("Set Graph aesthetics"),
        htmlP(
          "Changing the aesthetics does not reset the graph. It will take less time than regenerating the graph.",
          style = list(color = "gray")
        )
      ) 
    ),
    htmlDiv(
      className = "control",
      children = list(
        htmlH6('Background'),
        dccDropdown(
          id = "background",
          options = list(
            list(label = "black", value = "black"),
            list(label = "grey", value = bg),
            list(label = "white", value = "white")
          ),
          value = bg
        )
      )
    ),
    htmlDiv(
      className = "control",
      children = list(
        htmlH6('Colour map'),
        dccDropdown(
          id = "cmap",
          options = list(list(label = "fire", value = "fire"),
                         list(label = "viridis", value = "viridis"),
                         list(label = "blue", value = "blue")),
          value = 'fire'
        )
      )
    ),
    htmlDiv(
      className = "control",
      children = list(
        htmlH6('Raster Scaling'),
        dccRadioItems(
          id = "scaling",
          options = list(list(label = "Log transformation", value = "log"),
                         list(label = "Scale to [0,1]", value = "scale"),
                         list(label = "Original", value = "origin")),
          value = 'log'
        )
      )
    )
  )
)

app$layout(
  htmlDiv(
    className = "container scalable",
    children = list(
      htmlDiv(
        id = "banner-div",
        className = "banner",
        children = list(
          htmlH1(id = "appTitle", children = "United States Census"),
          htmlImg(src = "assets/plotly_logo_white.png")
        )
      ),
      htmlBr(),
      htmlDiv(
        id = "app-body",
        children = list(
          htmlDiv(
            id = "left-div",
            children = list(
              htmlDiv(
                id = "controls-div",
                children = list(
                  dccTabs(
                    children = list(
                      tab1,
                      tab2
                    )
                  )
                )
              ),
              dccStore(id = "store"),
              htmlDiv(
                className = "plot-outer",
                id = "race_outer",
                dccGraph(
                  id = "race_barplot"
                )
              )
            )
          ),
          htmlDiv(
            id = "rasterizer-outer",
            dccGraph(
              id = "rasterizer"
            )
          )
        )
      )
    ) 
  )
)

app$callback(
  output = list(id = 'store', property = 'data'),
  params = list(input(id = 'region', property = 'value'), 
                input(id = 'rasterizer', property = 'relayoutData'),
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
      rasterizer::canvas(mapping = mapping,
             plot_height = plot_height,
             plot_width = plot_width,
             x_range = x_range,
             y_range = y_range,
             show_raster = FALSE) %>%
      rasterizer::aggregation_points(reduction_func = reduction_func) %>%
      rasterizer::rasterizer() -> ds

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
  output = list(id = 'rasterizer', property = 'figure'),
  params = list(input(id = 'store', property = 'data'),
                input(id = 'cmap', property = 'value'),
                input(id = 'scaling', property = 'value'),
                input(id = 'reduc', property = 'value'),
                input(id = 'background', property = 'value')),
  function(store, cmap, scaling, reduc, background) {
    
    color <- if(cmap == "blue") {
      c("lightblue", "darkblue")
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
    
    agg <- lapply(store[[1]], unlist) %>% 
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
          x = seq(store[[2]][[1]], store[[2]][[2]], length.out = plot_width),
          y = seq(store[[3]][[1]], store[[3]][[2]], length.out = plot_height),
          z = agg,
          type = "heatmap",
          colorscale = colorscale
        )
      ),
      layout = list(
        xaxis = list('title' = 'longitude'),
        yaxis = list('title' = 'latitude'),
        margin = list('l' = 60, 'b' = 60, 't' = 10, 'r' = 10),
        plot_bgcolor = bg,
        paper_bgcolor = bg,
        hovermode = 'closest'
      )
    )
    return(l)
  }
)


app$callback(
  output = list(id = 'race_barplot', property = 'figure'),
  params = list(input(id = 'store', property = 'data')),
  function(store) {
    xlim <- unlist(store[[2]])
    ylim <- unlist(store[[3]])
    
    gg_color_hue <- function(n, l = 65, c = 100) {
      hues <- seq(15, 375, length = n + 1)
      grDevices::hcl(h = hues, l = l, c = c)[1:n]
    }
    
    range <- easting_northing(USA)
    colours <- gg_color_hue(length(race_full_names))
    
    if(sum(abs(xlim - range$easting)) < 1e-1 && sum(abs(ylim - range$northing)) < 1e-1) {
      
      barplot <- list(            
        data = list(
          list(
            x = race_full_names, 
            y = race_all,
            marker = list(color = colours),
            type = 'bar'
          )
        ),
        layout = list(
          title = "Race bar plot",
          plot_bgcolor = bg,
          paper_bgcolor = bg
        )
      )
      
    } else {
      race_sub_region <- data[easting <= xlim[2] & easting > xlim[1] & northing <= ylim[2] & northing > ylim[1], .N, by = race]
      
      barplot <- list(            
        data = list(
          list(
            x = race_full_names, 
            y = race_all,
            marker = list(color = colours),
            name = "All Regions",
            type = 'bar'
          ),
          list(
            x = race_full_names,
            y = race_table(names = race_full_names, 
                           table = race_sub_region),
            marker = list(color = gg_color_hue(length(race_full_names), l = 95)),
            name = "Sub Regions",
            type = 'bar'
          )
        ),
        layout = list(
          title = "Race barplot",
          plot_bgcolor = bg,
          paper_bgcolor = bg
        )
      )
    }
    
    return(barplot)
  }
)

if(appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
