library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(data.table)

# Plotly mapbox token
mapbox_access_token <- "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNqdnBvNDMyaTAxYzkzeW5ubWdpZ2VjbmMifQ.TXcBE-xg9BFdV2ocecc_7g"
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
  if (basename(getwd()) != "dashr-medical-provider-charges"){
    csv_path <- paste0("apps/dashr-medical-provider-charges", csv_path)
  }
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
    .SDcols = unlist(metric)]
  lat_lon_add <- df[,
    j = lapply(.SD, first),
    by = "Provider Name",
    .SDcols = c("lat", "lon", "Provider Street Address")]
  merge(aggRes, lat_lon_add, by = "Provider Name", all.x = TRUE)
}

buildUpperLeftPanel <- function(){
  htmlDiv(
    id = "upper-left",
    className = "six columns",
    children = list(
      htmlDiv(
        id = "state-select-outer",
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
        id = "select-metric-outer",
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
              values = list("All")
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
  plot_mapbox(
    data = filtered_data, x = ~lon, y = ~lat, 
    type = "scatter",
    mode = "markers",
    marker = list(
      color = ~f,
      showscale = TRUE,
      colorscale = list(
        c(0, "#21c7ef"),
        list(0.33, "#76f2ff"),
        list(0.66, "#ff6969"),
        list(1, "#ff1717")
      ),
      cmin = dMin,
      cmax = dMax,
      # Different scaling method used here (min max scaling)
      size = 10 * 
        (
          1 + (filtered_data[[co]] - dMin[[1]]) / (dMax[[1]] - dMin[[1]])
        ),
      colorbar = list(
        title = list(
          text = "Average Cost",
          font = list(color = "#737a8d", family = "Open Sans")
        ),
        titleside = "top",
        tickmode = "array",
        tickvals = c(dMin, dMax),
        ticktext = list(round(dMin, 2), round(dMax,2)),
        ticks = "outside",
        tickfont = list(family = "Open Sans", color = "#737a8d")
      )
    ),
    opacity = 0.8,
    #selectedpoints = selected_index,
    #selected = list(marker = list(color = "#ffff00")),
    #customdata = list(list(~"Provider Name", ~"Hospital Referral Region (HRR) Description")),
    hoverinfo = "text",
    text = paste(
      filtered_data[["Provider Name"]],
      filtered_data[["Hospital Referral Region (HRR) Description"]],
      sprintf("Average Procedure Cost: $%s", round(filtered_data[[co]], 2)),
      sep = "<br>"
    )
  ) %>%
    layout(
      margin = list(l = 10, r = 10, b = 10, t = 10, pad = 5),
      plot_bgcolor = "#171b26",
      paper_bgcolor = "#171b26",
      clickmode = "event+select",
      hovermode = "closest",
      showlegend = FALSE,
      mapbox = list(
        accesstoken = mapbox_access_token,
        bearing = 10,
        center = list(lat = ~mean(lat), lon = ~mean(lon)),
        pitch = 5,
        zoom = 7,
        style = "mapbox://styles/plotlymapbox/cjvppq1jl1ips1co3j12b9hex"
      )
    )
}

##############################################
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
          htmlP(
            className = "section title",
            children = paste(
              "Choose hospital on the map or procedures",
              "from the list below to see costs"
            )
          ),
          buildUpperLeftPanel(),
          htmlDiv(
            id = "geo-map-outer",
            className = "six columns",
            children = list(
              htmlP(
                id = "map-title",
                children = sprintf(
                  "Medical Provider Charges in the State of %s",
                  state_map[state_list[2]]
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
        id = "middle-container",
        className = "row",
        children = list(
          htmlDiv(
            id = "cost-stats-outer-container",
            children = list(
              htmlDiv(
                id = "table-left",
                className = "six columns",
                children = list(
                  htmlP("Hospital Charges Summary"),
                  dccLoading(
                    children = htmlDiv(id = "cost-stats-container")
                  )
                )
              ),
              htmlDiv(
                id = "table-right",
                className = "six columns",
                children = list(
                  htmlP("Procedure Charges Summary"),
                  dccLoading(
                    children = htmlDiv(id = "procedure-stats-container")
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
          id = "procedure plot",
          figure = plot_ly()
        )
      )
    )
  )
)

app$callback(
  output("geo-map", "figure"),
  list(
    input("metric-select", "value"),
    input("region-select", "value"),
    # input("procedure-plot", "selectedData"),
    input("state-select", "value")
  ),
  function(cost_select, region_select, state_select){
    state_agg_data <- generateAggregation(dataList[[state_select]], cost_metric)
    generateGeoMap(state_agg_data, cost_select, region_select)
  }
)

#state_agg_data <- generateAggregation(dataList[["AL"]], cost_metric)
#generateGeoMap(state_agg_data, "Average Covered Charges", list("AL - Birmingham", "AL - Dothan", "AL - Montgomery", "AL - Huntsville"))


app$run_server()
