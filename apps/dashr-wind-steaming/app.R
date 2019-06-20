
library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)
library(purrr)
library(dplyr)
library(DBI)
library(RSQLite)
library(glue)
library(data.table)
library(plotly)
library(readr)
library(VGAM)


appName <- Sys.getenv("DASH_APP_NAME")
pathPrefix <- sprintf("/%s/", appName)

Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
           DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

setwd("/app/apps/dashr-opiod-endemic")
setwd("/app/apps/dashr-wind-streaming")

DB_FILE = "/app/apps/dashr-wind-streaming/db/wind-data.db" 

get_wind_data <- function(start, end) {
  # Query wind data rows between two ranges
  # :params start: start row id
  # :params end: end row id
  # :returns: pandas dataframe object
  con <- dbConnect(SQLite(), DB_FILE)
  res <- dbSendQuery(con, glue("SELECT Speed, SpeedError, Direction FROM Wind WHERE rowid > {start} AND rowid <= {end}"))
  df <- dbFetch(res)
  return(df)
}

get_wind_data_by_id <- function(id){
  # Query a row from the Wind Table
  # :params id: a row id
  # :returns: pandas dataframe object
  con <- dbConnect(SQLite(), DB_FILE)
  res <- dbSendQuery(con, glue("SELECT * FROM Wind WHERE rowid = {id}"))
  df <- dbFetch(res)
  return(df)
}

get_current_time <- function(){
  # Helper function to get the current time in seconds
  total_time <- as.integer(as.numeric(Sys.time()) / 216000)
  return(total_time)
}

app <- Dash$new()

app_color <- list("graph_bg" = "#082255", "graph_line" = "#007ACE")

app$layout(htmlDiv(
  list(
    htmlDiv(
      list(
        htmlDiv(
          list(
            htmlH4("WIND SPEED STREAMING", className = "app__header__title"),
            htmlP(
              "This app continually queries a SQL database and displays live charts of wind speed and wind direction.",
              className = "app__header__title--grey"
            ),
            htmlButton(id = "learn-more", src = "https://plot.ly", 
                       style = list(color = "#082255"))
          ), className = "app__header__desc"
        ),
        htmlDiv(
          list(
            htmlImg(
              src  =  "https://raw.githubusercontent.com/plotly/dash-sample-apps/master/apps/dash-wind-streaming/assets/dash-logo-stripe-inverted.png",
              className = "app__menu__img"
            )
          ), className = "app__header__logo"
        )
      ), className = "app__header"
    ),
    htmlDiv(
      list(
        # wind speed
        htmlDiv(
          list(
            htmlDiv(
              list(htmlH6("WIND SPEED (MPH)", className = "graph__title"))
            ),
            dccGraph(
              id = "wind-speed"
            ),
            dccInterval(
              id = "wind-speed-update", interval = 1000, n_intervals = 1
            )
          ), className = "two-thirds column wind__speed__container"
        ),
        htmlDiv(
          list(
            # histogram
            htmlDiv(
              list(
                htmlDiv(
                  list(
                    htmlH6(
                      "WIND SPEED HISTOGRAM",
                      className = "graph__title"
                    )
                  )
                ),
                htmlDiv(
                  list(
                    dccSlider(
                      id = "bin-slider",
                      min = 1,
                      max = 60,
                      step = 1,
                      value = 20,
                      updatemode = "drag",
                      marks <- list(
                        setNames(
                          seq(from = 20, to = 60, by = 20),
                          seq(from = 20, to = 60, by = 20)
                        )
                      )
                    )
                  ), className = "slider"
                ),
                htmlDiv(
                  list(
                    dccChecklist(
                      id = "bin-auto",
                      options = list(
                        list("label" = "Auto", "value" = "Auto")
                      ),
                      value = list("Auto"),
                      inputClassName = "auto__checkbox",
                      labelClassName = "auto__label"
                    ),
                    htmlP(
                      "# of Bins:  Auto",
                      id = "bin-size",
                      className = "auto__p"
                    )
                  ), className = "auto__container"
                ),
                dccGraph(
                  id = "wind-histogram"
                )
              ), className = "graph__container first"
            ),
            # wind direction
            htmlDiv(
              list(
                htmlDiv(
                  list(
                    htmlH6(
                      "WIND DIRECTION", className = "graph__title"
                    )
                  )
                ),
                dccGraph(
                  id = "wind-direction"
                )
              ), className = "graph__container second"
            )
          ), className = "one-third column histogram__direction"
        )
      ), className = "app__content"
    )
  ), className = "app__container"
))

app$callback(
  output = list(id = "wind-speed", property = "figure"),
  params = list(input(id = "wind-speed-update", property = "n_intervals")),
  
  function(interval){
    # Generate the wind speed graph.
    # params interval :update the graph based on an interval
    
    total_time <- get_current_time()
    df <- get_wind_data((total_time - 200), total_time)
    
    fig <- plot_ly(df, y = df$Speed, 
                   type = "scatter", 
                   mode = "lines",
                   line = list(color = "#42C4F7"), 
                   hoverinfo = "skip",
                   height = 700,
                   error_y = list(
                     type = "data",
                     array = df$SpeedError,
                     thickness = 1.5,
                     width = 2, 
                     color = "#B4E8FC")) %>%
      layout(
        plot_bgcolor = app_color[["graph_bg"]],
        paper_bgcolor = app_color[["graph_bg"]],
        font = list("color" = "#fff"),
        xaxis <- list(
          range = list(0, 200),
          showline = TRUE,
          zeroline = FALSE,
          fixedrange = TRUE,
          tickvals = list(0, 50, 100, 150, 200),
          ticktext = list("200", "150", "100", "50", "0"),
          title = "Time Elapsed (sec)"
        ),
        yaxis <- list(
          range = list(
            min(0, min(df$Speed)),
            max(45, max(df$Speed) + max(df$SpeedError))
          ),
          showgrid = TRUE,
          showline = TRUE,
          fixedrange = TRUE,
          zeroline = FALSE,
          gridcolor = app_color[["graph_line"]],
          nticks = max(6, round(tail(df$Speed, n = 1)) / 10)
        )
      )
    return(fig)
  }
)

app$callback(
  output = list(id = "wind-direction", property = "figure"),
  params = list(input(id = "wind-speed-update", property = "n_intervals")),
  
  function(interval){
    # Generate the wind direction graph.
    # params interval: update the graph based on an interval
    
    total_time <- get_current_time()
    df <- get_wind_data_by_id(total_time)
    val <- tail(df$Speed, n=1)
    direction <- list(0, (df$Direction[1] - 20), (df$Direction[1] + 20), 0)
    
    fig <- plot_ly(
      type = "scatterpolar",
      mode = "lines",
      height = 350
    ) %>%
      add_trace(
        r = list(0, val, val, 0),
        theta = direction,
        fillcolor = "#084E8A",
        fill = "toself",
        line = list(color = "rgba(32, 32, 32, .6)", "width" = 1)
      ) %>%
      add_trace(
        r = list(0, val * 0.65, val * 0.65, 0),
        theta = direction,
        fillcolor = "#B4E1FA",
        fill = "toself",
        line = list(color = "rgba(32, 32, 32, .6)", "width" = 1)
      ) %>%
      add_trace(
        r = list(0, val * 0.3, val * 0.3, 0),
        theta = direction,
        fillcolor = "#EBF5FA",
        fill = "toself",
        line = list(color = "rgba(32, 32, 32, .6)", "width" = 1)
      ) %>% 
      layout(
        plot_bgcolor = app_color[["graph_bg"]],
        paper_bgcolor=app_color[["graph_bg"]],
        font = list("color" = "#fff"),
        autosize = FALSE,
        polar = list(
          "bgcolor" = app_color[["graph_line"]],
          "radialaxis" = list("range" = c(0, 45), "angle" = 45, "dtick" = 10),
          "angularaxis" = list("showline" = FALSE, "tickcolor" = "white")
        ),
        showlegend = FALSE
      )
    return(fig)
  }
)

app$callback(
  output = list(id = "wind-histogram", property = "figure"),
  params = list(
    input(id = "wind-speed-update", property = "n_intervals"),
    list(state(id = "wind-speed", property = "figure"),
         state(id = "bin-slider", property = "value"),
         state(id = "bin-auto", property = "values"))
  ),
  
  function(interval, wind_speed_figure, slider_value, auto_state) {

    try(
      expr = {
        if(!is.null(wind_speed_figure)){
          wind_val <- wind_speed_figure[['x']]$attrs[[1]]$y 
        }
        if("Auto" %in% auto_state){
          bin_val <- hist(
            wind_val,
            breaks = round(max(unlist(wind_val)))
          )
        }else{
          bin_val <- hist(wind_val, breaks = slider_value)
        }
      },
      error = function(e){
        e <- suppressWarnings()
      }
    )
  
    avg_val <- as.numeric(sum(wind_val)) / length(wind_val)
    median_val <- median(unlist(wind_val))
    
    pdf_fitted <- drayleigh(bin_val[[4]], ((tail(bin_val[[1]], n= 1) - bin_val[[1]][1]) / 3))
    y_val <- c(0, (pdf_fitted * max(bin_val[[2]]) * 20))
    y_val_max <- max(y_val)
    bin_val_max <- max(bin_val[[2]])
    
    fig <- plot_ly(
      x = bin_val[[1]][2:length(bin_val[[1]])],
      y = bin_val[[2]],
      marker = list(color = app_color[["graph_line"]]),
      showlegend = FALSE,
      hoverinfo = "x+y",
      type = "bar",
      height = 350
    ) %>%
      add_lines(
        x = bin_val$breaks, #list(bin_val[as.integer(length(bin_val) / 2)]),
        y = list(0),
        mode = "lines",
        line = list(dash = "dash", color = "#2E5266"),
        marker = list(opacity = 0),
        visible = TRUE,
        name = "Average"
      ) %>%
      add_lines(
        x = bin_val$breaks,
        y = list(0),
        mode = "lines",
        line = list(dash = "dot", color = "#BD9391"),
        marker = list(opacity = 0),
        visible = TRUE,
        name = "Median"
      ) %>%
      add_lines(
        x = bin_val$breaks,
        y = y_val,
        mode = "lines",
        line = list(color = "#42C4F7"),
        name = "Rayleigh Fit",
        marker = list(opacity = 0)
      )  %>%
      layout(
        plot_bgcolor = app_color[["graph_bg"]],
        paper_bgcolor = app_color[["graph_bg"]],
        font = list("color"= "#fff"),
        xaxis = list(
          title = "Wind Speed (mph)",
          showgrid = FALSE,
          showline = FALSE,
          fixedrange = TRUE
        ),
        yaxis = list(
          showgrid = FALSE,
          showline = FALSE,
          zeroline = FALSE,
          title = "Number of Samples",
          fixedrange = TRUE
        ),
        autosize = TRUE,
        bargap = 0.01,
        #bargroupgap = 0, 
        hovermode = "closest",
        legend = list(
          orientation = "h",
          yanchor = "bottom",
          xanchor = "center",
          y = 1,
          x = 0.5
        ),
        shapes = list(
          list(
            xref = "x",
            yref = "y",
            y1 = as.integer(max(bin_val_max, y_val_max)) + 0.5,
            y0 = 0,
            x0 = avg_val,
            x1 = avg_val,
            type = "line",
            line = list(dash = "dash", color = "#2E5266", width = 5)
          ),
          list(
            xref = "x",
            yref = "y",
            y1 = as.integer(max(bin_val_max, y_val_max)) + 0.5,
            y0 = 0,
            x0 = median_val,
            x1 = median_val,
            type = "line",
            line = list("dash"= "dot", "color"= "#BD9391", "width"= 5)
          )
        )
      )
    return(fig)
  }
)

app$callback(
  output = list(id = "bin-auto", property = "values"),
  params = list(
    input(id = "wind-speed", property = "figure")),
  
  function(wind_speed_figure){
    
    if(!length(wind_speed_figure[["x"]])){
      return("")
    }
    if(!is.null(wind_speed_figure) && (length(wind_speed_figure[["x"]]$attrs[[1]]$y) > 5)){
      return(list(""))
    }
    return(list("Auto"))
  }
)

app$callback(
  output = list(id = "bin-size", property = "children"),
  params = list(
    input(id = "bin-auto", property = "values"),
    state(id = "bin-slider", property = "value")
  ),
  
  function(autoValue, slider_value){
    if("Auto" %in% autovalue){
      return("# of Bins: Auto")
    }
    return(paste("# of Bins = ", as.character(as.integer(slider_value))))
  }
)

app$run_server(host = "0.0.0.0", port = Sys.getenv("PORT", 8080))