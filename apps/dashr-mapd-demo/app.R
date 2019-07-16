library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashTable)
library(plotly)
library(data.table)

library(fasttime)
# RJDBC requires rJava. See installation notes for Ubuntu:
# https://datawookie.netlify.com/blog/2018/02/installing-rjava-on-ubuntu/
library(RJDBC)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("/app/apps/%s", appName))
}

if (appName != ""){
  host <- Sys.getenv("DB_HOST")
} else {
  host <- "localhost"
}

if (host == "localhost"){
  # insert personal credentials here:
  user <- "admin"
  password <- "HyperInteractive"
  dbname <- "omnisci"
} else {
  # use server credentials:
  user <- "mapd"
  password <- "HyperInteractive"
  dbname <- "mapd"
}

flight_table <- "flights_2008_7M"

db_connect <- function(){
  tryCatch(
    {
      drv <- JDBC(
        "com.omnisci.jdbc.OmniSciDriver",
        # Insert pathname to the mapd / omnisci JDBC driver JAR 
        # file installed on your machine
        "/opt/omnisci/bin/omnisci-jdbc-4.7.1.jar",
        identifier.quote="'"
      )
      conn <- dbConnect(
        drv,
        sprintf("jdbc:%s:%s:6274:%s", dbname, host, dbname),
        user = user,
        password = password,
        host = host,
        dbname = dbname
      )
      if (!flight_table %in% dbListTables(conn)){
        message(
          sprintf(
            "Table %s not found in this database, please load sample data.",
            flight_table
          )
        )
      }
      conn
    },
    error = function(cond){
      message(sprintf("Error connectiong to OmniSci database: %s", cond))
    }
  )
}

conn <- db_connect()

generate_dest_choro <- function(dd_select, start, end){
  start_f <- sprintf("%s 00:00:00", start)
  end_f <- sprintf("%s 00:00:00", end)
  state_col <- ifelse(dd_select == "dep", "origin_state", "dest_state")
  tryCatch(
    {
      dest_df <- dbGetQuery(
        conn,
        sprintf(
          "SELECT AVG(depdelay) AS avg_delay, %s AS state
            FROM %s WHERE dep_timestamp BETWEEN '%s' AND '%s'
            GROUP BY %s",
          state_col, flight_table, start_f, end_f, state_col
        )
      )
      dest_df
    },
    error = function(cond){
      message(sprintf("Error querying for choropleth: ", cond))
      return(list())
    }
  )
  zmin = min(dest_df[, "avg_delay"])
  zmax = max(dest_df[, "avg_delay"])
  origin_or_destination <- ifelse(dd_select == "dep", "Origin", "Destination")
  plot_geo(
    data = dest_df,
    colorbar = list(
      tickvals = c(zmin, zmax),
      tickformat = ".2f",
      ticks = "",
      title = "Average Delay(Min)",
      titlefont = list(family = "Open Sans, sans-serif", color = "#515151"),
      thickness = 15,
      len = 1,
      tickcolor = "#151515",
      tickfont = list(family = "Open Sans, sans-serif", color = "#515151")
    ),
    colorscale = list(list(0, "#71cde4"), list(1, "#ecae50")),
    reversescale = TRUE,
    locations = ~state,
    z = ~avg_delay,
    locationmode = "USA-states",
    marker = list(line = list(color = "rgb(255,255,255)")),
    customdata = ~state
  ) %>%
    layout(
      title = list(
        text = sprintf(
          "Average Departure Delay <b>(Minutes)</b> By %s State",
          origin_or_destination
        ),
        font = list(
          family = "Open Sans, sans-serif",
          size = 15,
          color = "#515151"
        )
      ),
      margin = list(l = 10, r = 10, t = 50, b = 10),
      clickmode = "event+select",
      geo = list(
        scope = "usa", projection = list(type = "albers usa")
      )
    )
}

generate_flights_hm <- function(state, dd_select, start, end, select=FALSE){
  # Total flight count by Day of Week / Hour:
  # Initializing the dataframe with all possible days of week (1:7)
  hm <- data.frame("flight_dayofweek" = 1:7)
  state_query <- ""
  start_f <- sprintf("%s 00:00:00", start)
  end_f <- sprintf("%s 00:00:00", end)
  if (select){
    state_query <- sprintf("origin_state = '%s' AND ", state)
  }
  for (i in 0:23){
    hm_query <- sprintf(
      "SELECT flight_dayofweek, COUNT(*) AS Time_%s
        FROM %s WHERE %s%stime >= %s00 AND %stime < %s00
        AND %s_timestamp BETWEEN '%s' AND '%s'
        GROUP BY flight_dayofweek",
      i, flight_table, state_query,
      dd_select, i, dd_select, i+1,
      dd_select, start_f, end_f
    )
    tryCatch(
      {
        hm <- merge(
          hm, dbGetQuery(conn, hm_query),
          by = "flight_dayofweek", all = TRUE
        )
      },
      error = function(cond){
        message(sprintf("%s", cond))
      }
    )
  }
  hm[is.na(hm)] <- 0
  zmin = min(hm[, 2:ncol(hm)], na.rm = TRUE)
  zmax = max(hm[, 2:ncol(hm)], na.rm = TRUE)
  origin_or_destination <- ifelse(dd_select == "dep", "Origin", "Destination")
  plot_ly(
    data = hm,
    type = "heatmap",
    z = as.matrix(hm[, 2:ncol(hm)]),
    x = lapply(0:23, function(i) sprintf("%s:00", i)),
    y = list("Mon", "Tues", "Wed", "Thu", "Fri", "Sat", "Sun"),
    colorscale = list(list(0, "#71cde4"), list(1, "#ecae50")),
    reversescale = TRUE,
    showscale = TRUE,
    xgap = 2,
    ygap = 2,
    colorbar = list(
      len = 1,
      ticks = "",
      title = "Number of <br>Flights",
      titlefont = list(family = "Gravitas One, sans-serif", color = "#515151"),
      thickness = 15,
      tickcolor = "#515151",
      tickfont = list(family = "Open Sans, sans-serif", color = "#515151"),
      tickval = c(zmin, zmax)
    )
  ) %>%
    layout(
      title = list(
        text = sprintf(
          "%s Flights by days/hours in <b>%s</b>",
          origin_or_destination,
          ifelse(state != "", state, "the USA")
        ),
        font = list(family = "Open Sans, sans-serif", color = "#515151")
      ),
      font = list(family = "Open Sans, sans-serif", size = 13),
      margin = list(l = 100, r = 100, t = 100, b = 100)
    )
}

generate_time_series_chart <- function(state, start, end, dd_select){
  start_f <- sprintf("%s 00:00:00", start)
  end_f <- sprintf("%s 00:00:00", end)
  origin_or_destination <- ifelse(dd_select == "dep", "Departure", "Arrival")
  state_col <- ifelse(dd_select == "dep", "origin_state", "dest_state")
  state_query <- ifelse(
    nchar(state) > 0,
    sprintf("%s = '%s' AND ", state_col, state),
    ""
  )
  ts_query <- sprintf(
    "SELECT %s_timestamp AS ts_timestamp, %stime, %s
      FROM %s WHERE %s%s_timestamp BETWEEN '%s' AND '%s'",
    dd_select, dd_select, state_col, flight_table,
    state_query, dd_select, start_f, end_f
  )
  tryCatch(
    {
      df_ts <- dbGetQuery(conn, ts_query)
    },
    error = function(cond){
      message(sprintf("Error querying for time-series: %s", cond))
    }
  )
  # add 18000 seconds to make up for difference in tz
  df_ts[, "ts_timestamp"] <- fastPOSIXct(df_ts[,"ts_timestamp"]) + 18000
  min_date <- as.POSIXct(start, tz = "UTC")
  max_date <- as.POSIXct(end, tz = "UTC")
  x <- sort(cut(df_ts[, "ts_timestamp"], breaks = "hours"))
  y <- table(x)
  plot_ly(
    x = fastPOSIXct(names(y)) + 18000,
    y = y,
    type = "scatter",
    mode = "lines",
    line = list(color = "#71cde4")
  ) %>%
    layout(
      title = list(
        text = sprintf(
          "Flights by %s Time in <b>%s</b> <br> %s - %s",
          origin_or_destination,
          ifelse(state != "", state, "the USA"),
          start, end
        ),
        font = list(
          family = "Open sans, sans-serif",
          size = 15,
          color = "#515151"
        )
      ),
      font = list(family = "Open sans, sans-serif", size = 13),
      hovermode = "closest",
      margin = list(l = 75, r = 75, t = 100, b = 30),
      xaxis = list(
        rangeslider = list(visible = TRUE),
        range = c(min_date, max_date)
      ),
      yaxis = list(title = "Records")
    )
}

generate_count_chart <- function(state, dd_select, start, end){
  # Initializing the dataframe with all possible days of week (1:7)
  df_count <- data.frame("flight_dayofweek" = 1:7)
  select_f <- ifelse(dd_select == "dep", "origin_state", "dest_state")
  start_f <- sprintf("%s 00:00:00", start)
  end_f <- sprintf("%s 00:00:00", end)
  origin_or_destination <- ifelse(
    dd_select == "dep",
    "Originating in",
    "Arriving in"
  )
  state_query <- ifelse(
    nchar(state) > 0,
    sprintf("%s = '%s' AND ", select_f, state),
    ""
  )
  count_query <- sprintf(
    "SELECT flight_dayofweek, count(*) AS total_count
      FROM %s WHERE %s%s_timestamp BETWEEN '%s' AND '%s'
      GROUP BY flight_dayofweek",
    flight_table, state_query, dd_select, start_f, end_f
  )
  tryCatch(
    {
      df_count <- merge(
        df_count, dbGetQuery(conn, count_query),
        by = "flight_dayofweek", all = TRUE
      )
    },
    error = function(cond){
      message(sprintf("Error querying for count_chart: %s", cond))
      return(plot_ly())
    }
  )
  df_count[is.na(df_count)] <- 0
  plot_ly(
    data = df_count,
    x = ~total_count,
    y = list("Mon", "Tues", "Wed", "Thu", "Fri", "Sat", "Sun"),
    type = "bar",
    orientation = "h",
    marker = list(color = "#71cde4")
  ) %>%
    layout(
      title = list(
        text = sprintf(
          "Flight Counts by Days <br> %s <b>%s<b>",
          origin_or_destination,
          ifelse(state != "", state, "the USA")
        ),
        font = list(size = 15)
      ),
      font = list(family="Open Sans, sans-serif", size=13),
      xaxis = list(
        title = list(
          text = "Total Flight Counts",
          font = list(size = 12)
        ),
        zerolinecolor = "#999999"
      ),
      margin = list(l = 50, r = 50, t = 75, b = 75),
      clickmode = "event+select"
    )
}

generate_city_graph <- function(state_select, dd_select, start, end){
  start_f <- sprintf("%s 00:00:00", start)
  end_f <- sprintf("%s 00:00:00", end)
  origin_or_destination <- ifelse(dd_select == "dep", "Departure", "Arrival")
  days = list("Mon", "Tues", "Wed", "Thu", "Fri", "Sat", "Sun")
  count_df <- data.frame(matrix(nrow = 0, ncol = 1))
  colnames(count_df) <- c("city")
  state_f <- ifelse(nchar(state_select) > 0, state_select, "NY")
  state_col <- ifelse(dd_select == "arr", "origin_state", "dest_state")
  city_col <- ifelse(dd_select == "arr", "origin_city", "dest_city")
  for (i in 1:length(days)){
    count_query <- sprintf(
      "SELECT AVG(arrdelay) AS %s, %s AS city FROM %s
        WHERE %s = '%s' AND flight_dayofweek = %s
        AND %s_timestamp BETWEEN '%s' AND '%s' GROUP BY %s",
      days[[i]], city_col, flight_table, state_col,
      state_f, i, dd_select, start_f, end_f, city_col
    )
    tryCatch(
      {
        count_df <- merge(
          count_df, dbGetQuery(conn, count_query),
          by = "city", all = TRUE
        )
      },
      error = function(cond){
        message(sprintf("%s", cond))
      }
    )
  }
  count_df <- melt(
    count_df,
    id.vars = "city",
    variable.name = "Day",
    value.name = "Minutes"
  )
  plot_ly(
    count_df,
    x = ~Minutes,
    y = ~Day,
    color = ~city,
    type = "scatter",
    mode = "markers",
    customdata = ~city
  ) %>%
    layout(
      title = list(
        text = sprintf(
          "Avg %s Delays by City in <b>%s</b>",
          origin_or_destination, state_f
        ),
        font = list(
          family = "Open Sans, sans-serif",
          size = 15,
          color = "#515151"
        )
      ),
      font = list(family = "Open Sans, sans-serif", size = 13),
      xaxis = list(
        title = list(text = "Minutes", font = list(size = 12)),
        zerolinecolor = "#999999"
      ),
      yaxis = list(
        title = list(text = "Day", font = list(size = 12))
      ),
      margin = list(l = 75, r = 100, t = 70, b = 70),
      hovermode = "closest",
      clickmode = "event+select",
      dragmode = "select",
      showlegend = TRUE
    )
}

query_helper <- function(state_query, dd_select, start, end, weekday_query){
  #conn <- db_connect()
  origin_or_destination <- ifelse(dd_select == "dep", "origin", "dest")
  add_and <- ifelse(
    state_query!="",
    sprintf("%s_state = '%s' AND", origin_or_destination, state_query),
    ""
  )
  query <- sprintf(
    "SELECT uniquecarrier AS carrier,
    flightnum, dep_timestamp, arr_timestamp, origin_city, dest_city
    FROM %s WHERE %s %s_timestamp BETWEEN '%s 00:00:00' AND '%s 00:00:00' %s 
    LIMIT 100",
    flight_table, add_and, dd_select, start, end, weekday_query
  )
  tryCatch(
    {
      dff <- dbGetQuery(conn, query)
      dff[, "dep_timestamp"] <- fastPOSIXct(dff[, "dep_timestamp"]) + 18000
      dff[, "arr_timestamp"] <- fastPOSIXct(dff[, "arr_timestamp"]) + 18000
      dff[, "flightnum"] <- paste0(
        dff[, "carrier"],
        as.character(dff[, "flightnum"])
      )
      dff <- dff[, !names(dff) == "carrier"]
      return(dff)
    },
    error = function(cond){
      message(sprintf("Error Querying %s: %s", query, cond))
    }
  )
}

generate_flight_table <- function(...){
  # ... takes arguments from query_helper
  dashDataTable(
    id = "flights-table",
    columns = lapply(
      list(
        "flightnum",
        "dep_timestamp",
        "arr_timestamp",
        "origin_city",
        "dest_city"
      ),
      function(col){
        list(name = col, id = col)
      }
    ),
    filter_action = "native",
    data = df_to_list(
      query_helper(...)
    ),
    style_as_list_view = TRUE,
    style_header = list(
      textTransform = "Uppercase",
      fontWeight = "bold",
      backgroundColor = "#ffffff",
      padding = "10px 0px"
    ),
    style_cell_conditional = list(
      list(
        "if" = list(row_index = "odd"),
        backgroundColor = "#f5f6f7"
      )
    ),
    # TODO: Fix widths. See https://github.com/plotly/dash-table/issues/489
    style_cell = list(width = "240px"),
    fixed_rows = list(headers = TRUE, data = 0),
    style_table = list(maxHeight = "36rem", overflowY = "scroll")
  )
}

app <- Dash$new(name = "DashR MapD Demo")

app$layout(
  htmlDiv(
    className = "container scalable",
    children = list(
      htmlDiv(
        id = "banner",
        className = "banner",
        children = list(
          htmlH6("US Flights"),
          htmlImg(src = "assets/plotly_logo.png")
        )
      ),
      htmlDiv(
        className = "app_main_content",
        children = list(
          htmlDiv(
            id = "dropdown-select-outer",
            children = list(
              htmlDiv(
                className = "selector",
                children = list(
                  htmlP("Select Departure/Arrival"),
                  dccDropdown(
                    id = "dropdown-select",
                    options = list(
                      list(label = "Departure", value = "dep"),
                      list(label = "Arrival", value = "arr")
                    ),
                    value = "dep"
                  )
                )
              ),
              htmlDiv(
                id = "date-picker-outer",
                className = "selector",
                children = list(
                  htmlP("Select Date Range"),
                  dccDatePickerRange(
                    id = "date-picker-range",
                    min_date_allowed = as.POSIXct("2008-01-01"),
                    max_date_allowed = as.POSIXct("2008-12-31"),
                    initial_visible_month = as.POSIXct("2008-01-01"),
                    display_format = "MMM Do, YY",
                    start_date = as.POSIXct("2008-01-01"),
                    end_date = as.POSIXct("2008-01-08")
                  )
                )
              )
            )
          ),
          htmlDiv(
            id = "top-row",
            className = "row",
            children = list(
              htmlDiv(
                id = "map_geo_outer",
                className = "seven columns",
                children = dccLoading(
                  children = dccGraph(id = "choropleth")
                )
              ),
              htmlDiv(
                id = "flights_by_day_hm_outer",
                className = "five columns",
                children = dccLoading(
                  children = dccGraph(id = "flights_hm")
                )
              )
            )
          ),
          htmlDiv(
            id = "middle-row",
            className = "row",
            children = list(
              htmlDiv(
                id = "Flights-by-city-outer",
                className = "six columns",
                children = dccLoading(
                  children = dccGraph(id = "value_by_city_graph")
                )
              ),
              htmlDiv(
                id = "time-series-outer",
                className = "six columns",
                children = dccLoading(
                  children = dccGraph(
                    id = "flights_time_series",
                    figure = generate_time_series_chart(
                      "",
                      "2008-01-01",
                      "2008-01-08",
                      "dep"
                    )
                  )
                )
              )
            )
          ),
          htmlDiv(
            id = "bottom-row",
            className = "row",
            children = list(
              htmlDiv(
                id = "Count_by_days_outer",
                className = "four columns",
                children = dccLoading(
                  children = dccGraph(id = "count_by_day_graph")
                )
              ),
              htmlDiv(
                id = "flight_info_table_outer",
                className = "eight columns",
                children = dccLoading(
                  id = "table-loading",
                  children = generate_flight_table(
                    "", "dep", "2008-01-01", "2008-01-08", ""
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

app$callback(
  output("choropleth", "figure"),
  list(
    input("dropdown-select", "value"),
    input("date-picker-range", "start_date"),
    input("date-picker-range", "end_date")
  ),
  function(dd_select, start, end){
    generate_dest_choro(dd_select, start, end)
  }
)

app$callback(
  output("flights-table", "data"),
  list(
    input("flights_time_series", "relayoutData"),
    input("count_by_day_graph", "clickData"),
    input("value_by_city_graph", "selectedData"),
    input("choropleth", "figure"),
    state("dropdown-select", "value"),
    state("date-picker-range", "start_date"),
    state("date-picker-range", "end_date"),
    state("choropleth", "clickData"),
    state("flights-table", "data")
  ),
  function(
    ts_select,
    count_click,
    city_select,
    choro_fig,
    dd_select,
    start,
    end,
    choro_click,
    old_table_data
  ){
    ctx <- app$callback_context()
    inputs <- ctx$inputs
    states <- ctx$states
    prop_id <- strsplit(
      as.character(ctx[["triggered"]]["prop_id"]), "\\."
    )[[1]][1]
    state_query <- ""
    wk_map <- list(
      Mon = 1, Tues = 2, Wed = 3, Thu = 4, Fri = 5, Sat = 6, Sun = 7
    )
    if (!is.null(unlist(choro_click))){
      state <- choro_click[["points"]][[1]][["location"]]
      state_query <- state
    }
    if (prop_id == "choropleth"){
      return(df_to_list(query_helper(state_query, dd_select, start, end, "")))
    } else if (prop_id == "flights_time_series") {
      if (!"xaxis.range[0]" %in% names(
        inputs[["flights_time_series.relayoutData"]])
      ){
        return(old_table_data)
      }
      range_min <- inputs[[
        "flights_time_series.relayoutData"
        ]]["xaxis.range[0]"]
      range_max <- inputs[[
        "flights_time_series.relayoutData"
        ]]["xaxis.range[1]"]
      return(
        df_to_list(
          query_helper(state_query, dd_select, range_min, range_max, "")
        )
      )
    } else if (prop_id == "count_by_day_graph"){
      wk_day <- wk_map[[inputs[["count_by_day_graph.clickData"]]
                      [["points"]][[1]][["y"]]]]
      wk_day_query = sprintf(" AND flight_dayofweek = %s", wk_day)
      return(
        df_to_list(
          query_helper(state_query, dd_select, start, end, wk_day_query)
        )
      )
    } else if (prop_id == "value_by_city_graph"){
      cities <- lapply(city_select[["points"]], function(i) i[["customdata"]])
      wk_days <- lapply(city_select[["points"]], function(i) i[["y"]])
      q_template <- "SELECT uniquecarrier AS carrier,
                      flightnum, dep_timestamp,
                      arr_timestamp, origin_city, dest_city FROM %s
                      WHERE %s_timestamp BETWEEN '%s 00:00:00' and '%s 00:00:00'
                      AND flight_dayofweek = %s
                      AND %s = '%s' LIMIT 25"
      city_col <- ifelse(dd_select == "dep", "origin_city", "dest_city")

      DF <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(DF) <- c(
        "carrier",
        "flightnum",
        "dep_timestamp",
        "arr_timestamp",
        "origin_city",
        "dest_city"
      )
      for (wk_day in wk_days){
        for (city in cities){
          q <- sprintf(
            q_template,
            flight_table,
            dd_select,
            start,
            end,
            wk_map[[wk_day]],
            city_col,
            city
          )
          tryCatch(
            {
              dff <- dbGetQuery(conn, q)
              dff[, "flightnum"] <- paste0(
                dff[, "carrier"],
                as.character(dff[, "flightnum"])
              )
              dff <- dff[, !names(dff) == "carrier"]
              DF <- rbind(DF, dff)
            },
            error = function(cond){
              message(sprintf("Error Querying %s: %s", q, cond))
            }
          )
        }
      }
      if (nrow(DF) < 1){
        return(old_table_data)
      } else {
        return(df_to_list(DF))
      }
    }
  }
)

app$callback(
  output("flights_hm", "figure"),
  list(
    input("choropleth", "clickData"),
    input("choropleth", "figure"),
    state("dropdown-select", "value"),
    state("date-picker-range", "end_date"),
    state("date-picker-range", "start_date")
  ),
  function(choro_click, choro_figure, dd_select, end, start){
    if (!is.null(unlist(choro_click))){
      state <- lapply(
        choro_click[["points"]],
        function(point){point[["location"]]}
      )
      return(
        generate_flights_hm(state[[1]], dd_select, start, end, select = TRUE)
      )
    }
    generate_flights_hm("", dd_select, start, end, select = FALSE)
  }
)

app$callback(
  output("flights_time_series", "figure"),
  list(
    input("choropleth", "clickData"),
    input("choropleth", "figure"),
    state("dropdown-select", "value"),
    state("date-picker-range", "end_date"),
    state("date-picker-range", "start_date")
  ),
  function(choro_click, choro_figure, dd_select, end, start){
    if (!is.null(unlist(choro_click))){
      state <- lapply(
        choro_click[["points"]],
        function(point){point[["location"]]}
      )
      return(generate_time_series_chart(state[[1]], start, end, dd_select))
    }
    generate_time_series_chart("", start, end, dd_select)
  }
)

app$callback(
  output("count_by_day_graph", "figure"),
  list(
    input("choropleth", "clickData"),
    input("choropleth", "figure"),
    state("dropdown-select", "value"),
    state("date-picker-range", "end_date"),
    state("date-picker-range", "start_date")
  ),
  function(choro_click, choro_figure, dd_select, end, start){
    if (!is.null(unlist(choro_click))){
      state <- lapply(
        choro_click[["points"]],
        function(point){point[["location"]]}
      )
      return(generate_count_chart(state[[1]], dd_select, start, end))
    }
    generate_count_chart("", dd_select, start, end)
  }
)

app$callback(
  output("value_by_city_graph", "figure"),
  list(
    input("choropleth", "clickData"),
    input("choropleth", "figure"),
    state("dropdown-select", "value"),
    state("date-picker-range", "end_date"),
    state("date-picker-range", "start_date")
  ),
  function(choro_click, choro_figure, dd_select, end, start){
    if (!is.null(unlist(choro_click))){
      state <- lapply(
        choro_click[["points"]],
        function(point){point[["location"]]}
      )
      return(generate_city_graph(state[[1]], dd_select, start, end))
    }
    generate_city_graph("", dd_select, start, end)
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}
