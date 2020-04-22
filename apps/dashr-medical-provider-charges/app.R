library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashTable)
library(plotly)
library(data.table)

appName <- Sys.getenv("DASH_APP_NAME")

if (!appName == "") {
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

# Plotly mapbox token
mapbox_access_token <- "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNrOWJqb2F4djBnMjEzbG50amg0dnJieG4ifQ.Zme1-Uzoi75IaFbieBDl3A"
Sys.setenv("MAPBOX_TOKEN" = mapbox_access_token)

state_map <- list(
  "AK" = "Alaska",
  "AL" = "Alabama",
  "AR" = "Arkansas",
  "AZ" = "Arizona",
  "CA" = "California",
  "CO" = "Colorado",
  "CT" = "Connecticut",
  "DC" = "District of Columbia",
  "DE" = "Delaware",
  "FL" = "Florida",
  "GA" = "Georgia",
  "HI" = "Hawaii",
  "IA" = "Iowa",
  "ID" = "Idaho",
  "IL" = "Illinois",
  "IN" = "Indiana",
  "KS" = "Kansas",
  "KY" = "Kentucky",
  "LA" = "Louisiana",
  "MA" = "Massachusetts",
  "MD" = "Maryland",
  "ME" = "Maine",
  "MI" = "Michigan",
  "MN" = "Minnesota",
  "MO" = "Missouri",
  "MS" = "Mississippi",
  "MT" = "Montana",
  "NC" = "North Carolina",
  "ND" = "North Dakota",
  "NE" = "Nebraska",
  "NH" = "New Hampshire",
  "NJ" = "New Jersey",
  "NM" = "New Mexico",
  "NV" = "Nevada",
  "NY" = "New York",
  "OH" = "Ohio",
  "OK" = "Oklahoma",
  "OR" = "Oregon",
  "PA" = "Pennsylvania",
  "RI" = "Rhode Island",
  "SC" = "South Carolina",
  "SD" = "South Dakota",
  "TN" = "Tennessee",
  "TX" = "Texas",
  "UT" = "Utah",
  "VA" = "Virginia",
  "VT" = "Vermont",
  "WA" = "Washington",
  "WI" = "Wisconsin",
  "WV" = "West Virginia",
  "WY" = "Wyoming"
)

state_list <- names(state_map)

dataList <- list()
for (st in state_list){
  csv_path <- sprintf("data/processed/df_%s_lat_lon.csv", st)
  dataList[[st]] <- fread(csv_path)
}

cost_metric <- list(
  "Average Covered Charges",
  "Average Total Payments",
  "Average Medicare Payments"
)

init_region <- as.list(
  unique(
    dataList[[state_list[2]]][,
      "Hospital Referral Region (HRR) Description"
      ]
  )
)[[1]]

generateAggregation <- function(df, metric){
  agg <- function(x) list(min = min(x), mean = mean(x), max = max(x))
  aggRes <- df[,
    j = as.list(unlist(lapply(.SD, agg))),
    by = c("Hospital Referral Region (HRR) Description", "Provider Name"),
    .SDcols = unlist(metric)
    ]
  lat_lon_add <- df[,
    j = lapply(.SD, first),
    by = "Provider Name",
    .SDcols = c("lat", "lon", "Provider Street Address")
    ]
  merge(aggRes, lat_lon_add, by = "Provider Name", all.x = TRUE)
}

buildUpperLeftPanel <- function(){
  htmlDiv(
    id = "upper-left",
    className = "six columns",
    children = list(
      htmlH5(
        className = "section-title",
        paste(
          "Choose hospitals on the map or",
          "procedures from the list below to see costs"
        )
      ),
      htmlDiv(
        id = "select-state-metric",
        children = list(
          htmlDiv(
            id = "state-select-outer",
            style = list(width = "45%"),
            children = list(
              htmlLabel("Select a State"),
              dccDropdown(
                id = "state-select",
                options = lapply(
                  state_list,
                  function(x){
                    list(label = x, value = x)
                  }
                ),
                value = state_list[2]
              )
            )
          ),
          htmlDiv(
            id = "metric-select-outer",
            style = list(width = "45%"),
            children = list(
              htmlLabel("Choose a Cost Metric:"),
              dccDropdown(
                id = "metric-select",
                options = lapply(
                  cost_metric,
                  function(x){
                    list(label = x, value = x)
                  }
                ),
                value = cost_metric[[1]]
              )
            )
          )
        )
      ),
      htmlDiv(
        id = "region-select-outer",
        children = list(
          htmlLabel("Pick a Region:"),
          htmlDiv(
            id = "checklist-container",
            children = dccChecklist(
              id = "region-select-all",
              options = list(
                list(label = "Select All Regions", value = "All")
              ),
              value = list()
            )
          ),
          htmlDiv(
            id = "region-select-dropdown-outer",
            children = dccDropdown(
              id = "region-select",
              options = lapply(
                init_region,
                function(x){
                  list(label = x, value = x)
                }
              ),
              value = init_region[1:4],
              multi = TRUE,
              searchable = TRUE
            )
          )
        )
      ),
      htmlDiv(
        children = list(
          htmlDiv(
            id = "cost-stats-outer-container",
            children = list(
              htmlDiv(
                id = "table-left",
                className = "twelve columns",
                children = list(
                  htmlH5(
                    className = "section-title",
                    "Hospital Charges Summary",
                    style = list(marginTop = "1em")
                  ),
                  dccLoading(
                    children = htmlDiv(
                      id = "cost-stats-container",
                      children = generateDataTable(data.table(), "cost")
                    )
                  )
                )
              ),
              htmlDiv(
                id = "table-right",
                className = "twelve columns",
                children = list(
                  htmlH5(
                    className = "section-title",
                    "Procedure Charges Summary",
                    style = list(marginTop = "1em")

                  ),
                  dccLoading(
                    children = htmlDiv(
                      id = "procedure-stats-container",
                      children = generateDataTable(data.table(), "procedure")
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
}

generateGeoMap <- function(geo_data, selected_metric,
                           region_select, procedure_select){
  filtered_data <- geo_data[
    geo_data[["Hospital Referral Region (HRR) Description"]] %in% region_select
    ]
  co <- paste0(selected_metric, ".mean")
  dMin <- filtered_data[, lapply(.SD, min), .SDcols = co]
  dMax <- filtered_data[, lapply(.SD, max), .SDcols = co]
  dMid <- (dMin + dMax) / 2
  dLowMid <- (dMin + dMid) / 2
  dHighMid <- (dMid + dMax) / 2
  cuts <- unlist(c(dMin, dLowMid, dMid, dHighMid, dMax))
  filtered_data[, f := cut(
    filtered_data[[co]],
    cuts, labels = FALSE, include.lowest = TRUE
    )
  ]
 selected_indices <- which(
    filtered_data[["Provider Name"]] %in% procedure_select$hospital
  )
  plot_mapbox(
    data  = filtered_data,
    x = ~lon, y = ~lat,
    type = "scatter", mode = "markers",
    # Only selects single hospital if the length(selected_indices) > 1
    selectedpoints = c(-1, (selected_indices - 1)),
    selected = list(marker = list(color = "#FFFF00")),
    customdata = filtered_data[["Provider Name"]],
    marker = list(
      color = ~f,
      showscale = TRUE,
      colorscale = list(
        c(0, "#21c7ef"),
        list(0.33, "#76f2ff"),
        list(0.66, "#ff6969"),
        list(1, "#ff1717")
      ),
      # Different scaling method from python app (min max scaling)
      size = 15 *
        (
          1 + (filtered_data[[co]] - dMin[[1]]) / (dMax[[1]] - dMin[[1]])
        ),
      colorbar = list(
        title = list(
          text = "Average Cost",
          font = list(color = "#737a8d", family = "Arial", size = 17)
        ),
        tickmode = "array",
        tickvals = c(~min(f), ~max(f)),
        ticktext = sapply(
          c(dMin[[1]], dMax[[1]]),
          function(x){
            paste("$", round(x, 2))
          }
        ),
        ticks = "outside",
        tickfont = list(family = "Arial", color = "#737a8d"),
        thickness = 15,
        x = 0.8,
        len = 0.7
      )
    ),
    opacity = 0.8,
    hoverinfo = "text",
    text = paste(
      filtered_data[["Provider Name"]],
      filtered_data[["Hospital Referral Region (HRR) Description"]],
      sprintf("Average Procedure Cost: $%s", round(filtered_data[[co]], 2)),
      sep = "<br>"
    )
  ) %>%
    layout(
      margin  = list(l = 10, r = 10, b = 10, t = 10, pad = 5),
      plot_bgcolor = "#171b26",
      paper_bgcolor = "#171b26",
      clickmode = "event+select",
      hovermode = "closest",
      showlegend = FALSE,
      autosize = TRUE,
      mapbox = list(
        accesstoken = mapbox_access_token,
        bearing = 10,
        center = list(lat = ~mean(lat), lon = ~mean(lon)),
        pitch = 5,
        zoom = 5,
        style = "mapbox://styles/plotlymapbox/cjvppq1jl1ips1co3j12b9hex"
      )
    )
}

generateProcedurePlot <- function(raw_data, cost_select,
                                  region_select, provider_select){
  procedure_data <- raw_data[
    raw_data[["Hospital Referral Region (HRR) Description"]] %in% region_select
    ]
  selected_indices <- which(
    procedure_data[["Provider Name"]] %in% provider_select
  )
  plot_ly(
    data = procedure_data,
    y = procedure_data[["DRG Definition"]],
    x = procedure_data[[cost_select]],
    type = "box", boxpoints = "all", jitter = 0,
    pointpos = 0,
    hoveron = "points",
    fillcolor = "rgba(0,0,0,0)",
    line = list(color = "rgba(0,0,0,0)"),
    hoverinfo = "text",
    hovertext = paste(
      procedure_data[["Provider Name"]],
      paste("<b>", procedure_data[["DRG Definition"]], "</b>"),
      paste0("Average Procedure Cost: $", procedure_data[[cost_select]]),
      sep = "<br>"
    ),
    name = "",
    customdata = procedure_data[["Provider Name"]],
    # Same issue as generateGeoMap for the selected point indices
    selectedpoints = selected_indices - 1,
    selected = list(marker = list(color = "#FFFF00", size = 13)),
    unselected = list(marker = list(opacity = 0.2)),
    marker = list(
      line = list(width = 1, color = "#000000"),
      color = "#21c7ef",
      opacity = 0.7,
      symbol = "square",
      size = 12
    )
  ) %>%
    layout(
      showlegend = FALSE,
      hovermode = "closest",
      dragmode = "select",
      clickmode = "event+select",
      xaxis = list(
        zeroline = FALSE,
        automargin = TRUE,
        showticklabels = TRUE,
        title = list(text = "Procedure Cost", font = list(color = "#737a8d")),
        linecolor = "#737a8d",
        tickfont = list(color = "#737a8d"),
        type = "log"
      ),
      yaxis = list(
        automargin = TRUE,
        showticklabels = TRUE,
        tickfont = list(color = "#737a8d"),
        gridcolor = "#171b26"
      ),
      plot_bgcolor = "#171b26",
      paper_bgcolor = "#171b26"
    )
}

generateDataTable <- function(DT, type = c("procedure", "cost")){
  # Create datatable w/ no rows if DT is empty
  if (nrow(DT) > 0){
    d <- df_to_list(DT)
  } else {
    d <- list()
  }
  dashDataTable(
    id = ifelse(
      type == "cost",
      "cost-stats-table",
      "procedure-stats-table"
    ),
    columns = lapply(
      colnames(DT),
      function(x){
        list(name = x, id = x)
      }
    ),
    data = d,
    filter_action = "native",
    sort_action = ifelse(type == "cost", "none", "native"),
    sort_mode = "multi",
    page_action = "native",
    page_size = 5, page_current= 0,
    style_cell = list(
      backgroundColor = "#171b26",
      color = "#7b7d8d",
      overflow = "hidden",
      textOverflow = "ellipsis",
      whiteSpace = "normal"
    ),
    style_header = list(
      backgroundColor = "#1f2536"
    )
  )
}

app <- Dash$new(name = "DashR Medical Provider Charges")

app$layout(
  htmlDiv(
    list(
      htmlDiv(
        id = "banner",
        className = "banner",
        children = list(
          htmlImg(src = "assets/plotly_logo.png"),
          htmlH6("Dash Clinical Analytics")
        )
      ),
      htmlDiv(
        id = "upper-container",
        className = "row",
        children = list(
          buildUpperLeftPanel(),
          htmlDiv(
            id = "geo-map-outer",
            className = "six columns",
            children = list(
              htmlP(
                id = "map-title",
                children =
                  htmlH5(
                    className = "section-title",
                    sprintf(
                    "Medical Provider Charges in the State of %s",
                    state_map[state_list[2]]
                  )
                )
              ),
              htmlDiv(
                id = "geo-map-loading-outer",
                children = list(
                  dccLoading(
                    id = "loading",
                    children = dccGraph(
                      id = "geo-map",
                      figure = plot_ly() %>%
                        layout(
                          plot_bgcolor = "#171b26",
                          paper_bgcolor = "#171b26"
                        )
                    )
                  )
                )
              )
            )
          )
        )
      ),
      htmlDiv(
        id = "lower-container",
        children = dccGraph(
          id = "procedure-plot",
          figure = generateProcedurePlot(
            dataList[[2]], cost_metric[[1]], init_region, list()
          )
        )
      )
    )
  )
)

app$callback(
  output("region-select-dropdown-outer", "children"),
  list(input("state-select", "value")),
  function(state_select){
    state_raw_data <- dataList[[state_select]]
    regions <- unique(
      state_raw_data[["Hospital Referral Region (HRR) Description"]]
    )
    dccDropdown(
      id = "region-select",
      options = lapply(regions, function(x) list(label = x, value = x)),
      value = regions[1:4],
      multi = TRUE,
      searchable = TRUE
    )
  }
)

app$callback(
  output("map-title", "children"),
  list(input("state-select", "value")),
  function(state_select){
    htmlH5(
      className = "section-title",
      sprintf(
        "Medical Provider Charges in the State of %s",
        state_map[state_select]
      )
    )
  }
)

app$callback(
  output("region-select", "value"),
  list(
    input("region-select-all", "values"),
    state("region-select", "value"),
    state("region-select", "options")
  ),
  function(select_all, selected_options, options){
    if (length(select_all) > 0){
      return(lapply(
        options,
        function(x){
          x[["value"]]
        })
      )
    }
    return(selected_options)
  }
)

app$callback(
  output("cost-stats-container", "children"),
  list(
    input("geo-map", "selectedData"),
    input("procedure-plot", "selectedData"),
    input("metric-select", "value"),
    input("state-select", "value")
  ),
  function(geo_select, procedure_select, cost_select, state_select){
    state_agg <- generateAggregation(
      dataList[[state_select]], cost_metric
    )
    geoDataList <- list(
      "Provider Name" = list(),
      "City" = list(),
      "Street Address" = list(),
      "Maximum Cost ($)" = list(),
      "Minimum Cost ($)" = list()
    )
    # Create empty data table as default
    DT <- as.data.table(geoDataList)
    costmin <- paste0(cost_select, ".min")
    costmax <- paste0(cost_select, ".max")
    ctx <- app$callback_context()
    prop_id <- strsplit(
      as.character(ctx[["triggered"]]["prop_id"]), "\\."
      )[[1]][1]
    if ( (prop_id == "procedure-plot") &
        (!is.null(unlist(procedure_select)))){
      point_select <- procedure_select
    } else if ( (prop_id == "geo-map") &
               (!is.null(unlist(geo_select)))){
      point_select <- geo_select
    } else {
      return(generateDataTable(DT, "cost"))
    }
    for (i in seq_along(point_select[["points"]])){
      point <- point_select[["points"]][i]
      dfrow <- state_agg[
        state_agg[["Provider Name"]] == point[[1]]["customdata"]
        ]
      hrr <- dfrow[["Hospital Referral Region (HRR) Description"]]
      city <- strsplit(
        as.character(hrr), " - ", fixed = TRUE
      )[[1]][2]
      adrs <- dfrow[["Provider Street Address"]]
      geoDataList[["Provider Name"]][[i]] <- dfrow[["Provider Name"]]
      geoDataList[["City"]][[i]] <- city
      geoDataList[["Street Address"]][[i]] <- adrs
      geoDataList[["Maximum Cost ($)"]][[i]] <- dfrow[[costmax]]
      geoDataList[["Minimum Cost ($)"]][[i]] <- dfrow[[costmin]]
    }
    DT <- as.data.table(geoDataList)
    DT <- DT[!duplicated(DT[["Provider Name"]])]
    generateDataTable(DT, "cost")
  }
)

app$callback(
  output("procedure-stats-container", "children"),
  list(
    input("procedure-plot", "selectedData"),
    input("geo-map", "selectedData"),
    input("metric-select", "value"),
    state("state-select", "value")
  ),
  function(procedure_select, geo_select, cost_select, state_select){
    procedureList <- list(
      "DRG" = list(),
      "Procedure" = list(),
      "Provider Name" = list(),
      "Cost Summary" = list()
    )
    ctx <- app$callback_context()
    prop_id <- strsplit(
      as.character(ctx[["triggered"]]["prop_id"]),
      "\\."
      )[[1]][1]
    # Displays all the procedures offered at selected hospital
    if ( (prop_id == "geo-map") &
        (!is.null(unlist(geo_select)))){
      if (is.null(unlist(geo_select$points))){
        return(
          generateDataTable(
            as.data.table(procedureList), "procedure")
        )
      }
      provider_select <- list()
      i <- 1
      for (point in geo_select[["points"]]){
        provider_select[[i]] <- point[["customdata"]]
        i <- i + 1
      }
      state_raw_data <- dataList[[state_select]]
      filtered <- state_raw_data[
        state_raw_data[["Provider Name"]] %in% unlist(provider_select)
        ]
      for (i in 1:nrow(filtered)){
        fullProcedureName <- filtered[i, `DRG Definition`]
        splitProcedureName <- strsplit(
          as.character(fullProcedureName), " - ", fixed = TRUE
        )[[1]]
        drg <- splitProcedureName[1]
        proc <- splitProcedureName[2]
        procedureList[["DRG"]][[i]] <- drg
        procedureList[["Procedure"]][[i]] <- proc
        procedureList[["Provider Name"]][[i]] <- filtered[i, `Provider Name`]
        procedureList[["Cost Summary"]][[i]] <- filtered[i, ..cost_select][[1]]
      }
    } else if (!is.null(unlist(procedure_select))){
      for (i in seq_along(procedure_select[["points"]])){
        point <- procedure_select[["points"]][[i]]
        fullProcedureName <- point[["y"]]
        splitProcedureName <- strsplit(
          as.character(fullProcedureName), " - ", fixed = TRUE
        )[[1]]
        drg <- splitProcedureName[1]
        proc <- splitProcedureName[2]
        procedureList[["DRG"]][[i]] <- drg
        procedureList[["Procedure"]][[i]] <- proc
        procedureList[["Provider Name"]][[i]] <- point[["customdata"]]
        procedureList[["Cost Summary"]][[i]] <- point[["x"]]
      }
    }
    DT <- as.data.table(procedureList)
    generateDataTable(DT, "procedure")
  }
)

app$callback(
  output("geo-map", "figure"),
  list(
    input("metric-select", "value"),
    input("region-select", "value"),
    input("procedure-plot", "selectedData"),
    input("state-select", "value")
  ),
  function(cost_select, region_select, procedure_select, state_select){
    state_agg_data <- generateAggregation(
      dataList[[state_select]], cost_metric
    )
    provider_data <- list(procedure = list(), hospital = list())
    if (!is.null(unlist(procedure_select))){
      i <- 1
      for (point in procedure_select[["points"]]){
        provider_data[["procedure"]][[i]] <- point[["y"]]
        provider_data[["hospital"]][[i]] <- point[["customdata"]]
        i <- i + 1
      }
    }
    generateGeoMap(state_agg_data, cost_select, region_select, provider_data)
  }
)

app$callback(
  output("procedure-plot", "figure"),
  list(
    input("metric-select", "value"),
    input("region-select", "value"),
    input("geo-map", "selectedData"),
    input("state-select", "value")
  ),
  function(cost_select, region_select, geo_select, state_select){
    state_raw_data <- dataList[[state_select]]
    provider_select <- list()
    if (!is.null(unlist(geo_select))){
      i <- 1
      for (point in geo_select[["points"]]){
        provider_select[[i]] <- point[["customdata"]]
        i <- i + 1
      }
    }
    generateProcedurePlot(
      state_raw_data, cost_select, region_select, provider_select
    )
  }
)

if (!appName == ""){
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(debug = TRUE)
}

