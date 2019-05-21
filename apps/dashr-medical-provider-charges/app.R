library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(data.table)

# Plotly mapbox token
mapbox_access_token <- "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNqdnBvNDMyaTAxYzkzeW5ubWdpZ2VjbmMifQ.TXcBE-xg9BFdV2ocecc_7g"

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

cost_metric = list(
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
            children = "Choose hospital on the map or procedures from the list below to see costs"
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
                  state_map[state_list[1]]
                )
              ),
              htmlDiv(
                id = "geo-map-loading-outer",
                children = list(
                  dccLoading(
                    id = "loading",
                    children = dccGraph(
                      id = "geo_map",
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

app$run_server()
