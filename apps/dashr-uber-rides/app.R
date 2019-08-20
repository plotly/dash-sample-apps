library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(data.table)
library(dplyr)
library(lubridate)
library(fasttime)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf("/app/apps/%s", appName))
}
#################################### LOAD DATA & CREATE GLOBAL OBJECTS #############################

ridesPart1 <- "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data1.csv"
ridesPart2 <- "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data2.csv"
ridesPart3 <- "https://raw.githubusercontent.com/plotly/datasets/master/uber-rides-data3.csv" 
# Links for datasets

ridesRaw_1 <- data.table::fread(ridesPart1, stringsAsFactors = FALSE)
ridesRaw_2 <- data.table::fread(ridesPart2, stringsAsFactors = FALSE)
ridesRaw_3 <- data.table::fread(ridesPart3, stringsAsFactors = FALSE)
# Read data

ridesDf <- rbind(ridesRaw_1, ridesRaw_2, ridesRaw_3)
# Combine partitions of data

names(ridesDf) <- c("dateTime", "Lat", "Lon")
# Rename columns


ridesDf <- ridesDf %>% 
  dplyr::mutate(
    Date.Time = fasttime::fastPOSIXct(dateTime) + 
      (as.POSIXct(dateTime[1], format = "%Y-%m-%d %H:%M:%OS") - 
         fasttime::fastPOSIXct(dateTime[1])
      )
  ) %>% 
  dplyr::as_tibble()

locationCoordinates <- stats::setNames(
  list(
    list(40.7272, -73.991251), list(40.7505, -73.9934),
    list(40.8296, -73.9262), list(40.7484, -73.9857),
    list(40.7069, -74.0113), list(40.644987, -73.785607),
    list(40.7527, -73.9772), list(40.7589, -73.9851)
  ),
  c("default", "Madison Square Garden", "Yankee Stadium",
    "Empire State Building", "New York Stock Exchange",
    "JFK Airport", "Grand Central Station", "Times Square")
)

# Location coordinates for map-graph

colorMapOrj <- c("#F4EC15", "#DAF017", "#BBEC19", "#9DE81B",
                 "#80E41D", "#66E01F", "#4CDC20", "#34D822",
                 "#24D249", "#25D042", "#26CC58", "#28C86D",
                 "#29C481", "#2AC093", "#2BBCA4", "#2BB5B8",
                 "#2C99B4", "#2D7EB0", "#2D65AC", "#2E4EA4",
                 "#2E38A4", "#3B2FA0", "#4E2F9C", "#603099")
# Color codes for bars

csSeq <- seq(0, 1, length.out = 23)
# Create 24 evenly spaced floats from 0 to 1:

colorScale <- list(
  list(0, "#F4EC15"),
  list(csSeq[1], "#DAF017"), list(csSeq[2], "#BBEC19"),
  list(csSeq[3], "#9DE81B"), list(csSeq[4], "#80E41D"),
  list(csSeq[5], "#66E01F"), list(csSeq[6], "#4CDC20"),
  list(csSeq[7], "#34D822"), list(csSeq[8], "#24D249"),
  list(csSeq[9], "#25D042"), list(csSeq[10], "#26CC58"),
  list(csSeq[11], "#28C86D"), list(csSeq[12], "#29C481"),
  list(csSeq[13], "#2AC093"), list(csSeq[14], "#2BBCA4"),
  list(csSeq[15], "#2BB5B8"), list(csSeq[16], "#2C99B4"),
  list(csSeq[17], "#2D7EB0"), list(csSeq[18], "#2D65AC"),
  list(csSeq[19], "#2E4EA4"), list(csSeq[20],  "#2E38A4"),
  list(csSeq[21], "#3B2FA0"), list(csSeq[22], "#4E2F9C"),
  list(csSeq[23], "#603099")
)                
# Create custom colorScale list 
ridesDf <- ridesDf %>%
  dplyr::mutate(rideHour = as.numeric(as.POSIXlt(Date.Time)$hour))
# Extract hour from dates and append corresponding color codes as a column for coloring later
# Factor works too

mapboxToken <- paste0("pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNqdnBvNDMyaT",
                      "AxYzkzeW5ubWdpZ2VjbmMifQ.TXcBE-xg9BFdV2ocecc_7g")
# Mapbox token for plotly | this one is for plot_mapbox figure

Sys.setenv("MAPBOX_TOKEN" = mapboxToken)
# Setting mapbox token for R environment 

locationOptions <- list(
  list(label = "New York Main", value = "default"),
  list(label = "Madison Square Garden", value = "Madison Square Garden"),
  list(label = "Yankee Stadium", value = "Yankee Stadium"),
  list(label = "Empire State Building", value = "Empire State Building"),
  list(label = "New York Stock Exchange", value = "New York Stock Exchange"),
  list(label = "JFK Airport", value = "JFK Airport"),
  list(label = "Grand Central Station", value = "Grand Central Station"),
  list(label = "Times Square", value = "Times Square"))

allHours <- c(
  sapply(
    0:23,
    function(h){
      hr <- formatC(h, width = 2, flag = "0")
      paste0(" ", hr, ":00:00")
    }
  )
)

hourOptions <- c(
  list(list(label = "Entire Day", value = "entire_day")),
  lapply(
    0:23, 
    function(h){
      hr <- formatC(h, width = 2, flag = "0")
      list(label = paste0(hr , ":00"), value = paste0(" ", hr, ":00:00"))
    }
  )
)

####################################################################################################

#################################### APP START #####################################################

app <- Dash$new(name = "dashr-uber-rides")
# Initiate application

####################################################################################################

#################################### CREATE LAYOUT VARIABLES #######################################

pageTitle <- htmlH2("DASH - UBER DATA APP", className = "title-left")
firstP <- htmlP(paste("Select different dates, locations or ",
                      "hours using dropdowns below or by selecting ",
                      "different time frames on the histogram", sep = ""))
plotlyLogo <- htmlA(
  list(
    htmlImg(src = "assets/dash-logo-new.png",
            className = "logo")),
  href = "https://dashr-docs.herokuapp.com/")
totalRidesDay <- htmlDiv(id = "total-rides")
locationDropdown <- dccDropdown(id = "location-dropdown",
                                options = locationOptions,
                                value = "default",
                                className = "location-dropdown")

hourDropdown <- dccDropdown(id = "hour-dropdown",
                            options = hourOptions,
                            multi = TRUE,
                            value = "entire_day",
                            className = "hour-dropdown")

datePicker <- dccDatePickerSingle(
  id = "date-picker",
  date = lubridate::make_date(year = 2014L, month = 4L, day = 1L),
  display_format = "MMMM D, YYYY",
  min_date_allowed = lubridate::make_date(year = 2014L, month = 4L, day = 1L),
  max_date_allowed = lubridate::make_date(year = 2014L, month = 9L, day = 30L),
  className = "date-picker", style = list(display = "block"))

sourceObj <- dccMarkdown(
  paste(
    "Source: [FiveThirtyEight](https://github.com/fivethirtyeight/uber-",
    "tlc-foil-response/tree/master/uber-trip-data)", sep = ""),
  className = "source-obj")

####################################################################################################

#################################### CREATE LAYOUT #################################################

app$layout(
  htmlDiv(
    className = "row",
    list(
      htmlDiv(list(
        plotlyLogo,
        pageTitle,
        firstP,
        htmlDiv(datePicker, className = "div-for-dropdown"),
        htmlDiv(id = "hidden-date-picker", hidden = TRUE),
        htmlDiv(locationDropdown, className = "div-for-dropdown"),
        htmlDiv(hourDropdown, className = "div-for-dropdown"),
        totalRidesDay,
        sourceObj
      ), className = "four columns div-user-controls"
      ),
      
      htmlDiv(list(
        dccGraph(id = "map-graph", className = "map-graph"),
        dccGraph(id = "histogram", className = "histogram",
                 selectedData = list(points = list(), range = NULL),
                 clickData = list(points = list(), range = NULL))
      ), className = "div-for-charts eight columns"
      )
    )
  )
)

####################################################################################################

#################################### HELPER FUNCS FOR CALLBACKS ####################################

# Filters out the date selected by user from dataset
filterDay <- function(date) {
  
  datePaste <- paste(date, " 00:00:00", sep = "")
  # Combine month-day-year to %Y-%m-%d %H:%M:%OS format
  
  dayCurrent <- as.POSIXct(datePaste)
  dayNext <- as.POSIXct(datePaste) + 86400
  #Created today and tomorrows posix object for filtering
  
  dfDay <- ridesDf %>%
    dplyr::filter(Date.Time >= dayCurrent & Date.Time < dayNext) %>% 
    dplyr::arrange(`Date.Time`)
  # Filter out data
  
  return(dfDay)
}

# Returns hour list for callbacks
entireDayConverter <- function(hour) {
  
  hourVec <- unlist(hour)
  
  if ("entire_day" %in% hourVec) {
    hourVec <- allHours
    # Assign allHours vector to hour if "entire_day" selected
  }
  return(hourVec)
}

# Filters out the month-day-hour selected by user from dataset
dfToMap <- function(date, hour) {
  
  hourVec <- entireDayConverter(hour)
  if (length(hourVec) >= 24) {
    return(filterDay(date))
    # Avoid loop if entire day selected for faster filtering
  } else {
    
    dfEmpty <- data.table::data.table()
    
    for (h in hourVec) {
      
      hourPaste <- paste(date, h, sep = "")
      # Combine month-day-year to %Y-%m-%d %H:%M:%OS format
      
      hourStart <- as.POSIXct(hourPaste)
      hourEnd <- hourStart + 3600
      # Define df slicing interval
      
      dfEmpty <- ridesDf %>% 
        dplyr::filter(Date.Time >= hourStart & Date.Time < hourEnd) %>% 
        rbind(dfEmpty, .)
      # Append hours to get df to plot
    }
    return (dfEmpty)
  }
}

####################################################################################################

#################################### CALLBACKS START ###############################################

app$callback(output = list(id = "hidden-date-picker", property = "children"),
             params = list(input(id = "date-picker", property = "date")),
             
             function(date) {
               return (as.character(date))
             }
)

app$callback(output = list(id = "map-graph", property = "figure"),
             params = list(
               input(id = "hidden-date-picker", property = "children"),
               input(id = "hour-dropdown", property = "value"),
               input(id = "location-dropdown", property = "value"),
               state(id = "map-graph", property = "relayoutData")),
             
             function(date, hour, location, relayout) {
               subsetIndex <- which(names(locationCoordinates) %in% location)
               
               latInitial <- locationCoordinates[[subsetIndex]][[1]]
               lonInitial <- locationCoordinates[[subsetIndex]][[2]]
               zoom <- 15
               bearing <- 0
               # Set layout variables for locations
               
               if (location == "default") {
                 zoom <- 12
               }
               return (
                 plotly::plot_mapbox(
                   dfToMap(date, hour), lat = ~Lat, lon = ~Lon, split = ~rideHour,
                   mode = "markers", hoverinfo = "lat+lon+text",
                   marker = list(color = ~rideHour, cmin = 0, cmax = 23,
                                 colorscale = colorScale, size = 4)) %>% 
                   plotly::layout(title = "Uber Rides in NYC",
                                  autosize = TRUE,
                                  margin = list(l = 0, r = 0, t = 0, b = 0),
                                  showlegend = TRUE,
                                  legend = list(orientation = "v",
                                                bgcolor = "#1e1e1e", yanchor = "top"),
                                  font = list(color = "white"),
                                  plot_bgcolor = "#1e1e1e", paper_bgcolor = "#1e1e1e",
                                  mapbox = list(zoom = zoom,
                                                accesstoken = mapboxToken,
                                                # The plot DOES NOT render without providing token!
                                                style = "dark",
                                                center = list(lat = latInitial,
                                                              lon = lonInitial),
                                                bearing = bearing),
                                  updatemenus = list(
                                    list(
                                      buttons = list(
                                        list(
                                          args = list(list(
                                            "mapbox.zoom" = 12,
                                            "mapbox.center.lon" = "-73.991251",
                                            "mapbox.center.lat" = "40.7272",
                                            "mapbox.bearing" = 0,
                                            "mapbox.style" = "dark"
                                          )),
                                          label = "Reset Zoom",
                                          method = "relayout"
                                        )
                                      ),
                                      direction = "left",
                                      pad = list("r" = 0, "t" = 0, "b" = 0, "l" = 0),
                                      showactive = FALSE,
                                      type = "buttons",
                                      x = 0.45,
                                      xanchor = "left",
                                      yanchor = "bottom",
                                      bgcolor = "#323130",
                                      borderwidth = 1,
                                      bordercolor = "#6d6d6d",
                                      font = list(color = "#FFFFFF"),
                                      y = 0.02
                                    )
                                  )
                   )
               )
             }
)

app$callback(output = list(id = "histogram", property = "figure"),
             params = list(
               input(id = "hidden-date-picker", property = "children"),
               input(id = "hour-dropdown", property = "value")),
             
             function(date, hour) {
               
               # Aggregate number of hourly rides to df
               dfAgg <- filterDay(date) %>% 
                 dplyr::mutate(Date_Hours = format(Date.Time, "%Y-%m-%d %H")) %>% 
                 dplyr::group_by(Date_Hours) %>%
                 dplyr::summarise(Ride_Counts = n()) %>%
                 dplyr::transmute(Date_Hours = seq(0, 23), Ride_Counts = Ride_Counts)
               
               # Reformatting for getting the x-axis ticks correct
               numericHours <- as.numeric(substr(entireDayConverter(hour), start = 2, stop = 3))
               
               if (length(numericHours) > 0 && length(numericHours) < 24) {
                 colorMapOrj[numericHours + 1] <- "FFFFFF"
               }
               return (
                 plotly::plot_ly(
                   x = dfAgg$Date_Hours,
                   y = dfAgg$Ride_Counts,
                   name = "histogram",
                   marker = list(
                     color = colorMapOrj),
                   type = "bar") %>% 
                   plotly::layout(
                     title = list(
                       text = paste("Hold click and create rectangle ",
                                    "to section data by time", sep = ""),
                       x = 0.02,
                       y = 0.9,
                       font = list(color = "white", size = 12)),
                     plot_bgcolor = "#2d2c2b",
                     bargap = 0.01,
                     barmode = "group",
                     paper_bgcolor = "#2d2c2b",
                     margin = list(l = 10, r = 0, t = 30, b = 30),
                     dragmode = "select",
                     xaxis = list(range = list(-0.5, 23.5),
                                  nticks = 25,
                                  ticksuffix = ":00",
                                  showgrid = FALSE,
                                  fixedrange = TRUE,
                                  color = "white",
                                  type = "category"
                     ),
                     yaxis = list(showticklabels = FALSE,
                                  showgrid = FALSE,
                                  fixedrange = TRUE,
                                  rangemode = "nonnegative",
                                  zeroline = FALSE)) %>%
                   plotly::add_annotations(text = dfAgg$Ride_Counts,
                                           showarrow = FALSE,
                                           font = list(color = "white"),
                                           xanchor = "center",
                                           yanchor = "bottom")
               )
             }
)

app$callback(output = list(id = "hour-dropdown", property = "value"),
             params = list(
               input(id = "histogram", property = "selectedData"),
               input(id = "histogram", property = "clickData")),
             
             function(select, click) {
               
               dfSelectedHour <- data.table::rbindlist(select$points)
               # Convert select list [points] to df
               
               if (nrow(dfSelectedHour) > 0) {
                 dfClick <- data.table::data.table()
                 # Reset dfClick when selectedData
               } else {
                 dfClick <- data.table::rbindlist(click$points)
                 # Add clicked data to df
               }
               
               dfSelectedHour <- rbind(dfSelectedHour, dfClick)
               # Append selection and click df"s together
               
               if (nrow(dfSelectedHour) > 0){
                 hoursToReturn <- with(
                   dfSelectedHour,
                   ifelse(
                     pointNumber < 10,
                     paste(" 0", pointNumber, ":00:00", sep = ""),
                     paste(" ", pointNumber, ":00:00", sep = "")
                   )
                 )
                 # Create a new column with the right hour formatting
               } else {
                 hoursToReturn <- "entire_day"
               }
               return(hoursToReturn)
             }
)

app$callback(output = list(id = "total-rides", property = "children"),
             params = list(
               input(id = "hidden-date-picker", property = "children"),
               input(id = "hour-dropdown", property = "value")),
             
             function(date, h) {
               
               dfDay <- filterDay(date) %>% 
                 dplyr::mutate(rideHourNum = as.numeric(as.POSIXlt(Date.Time)$hour)) 
               
               # Filter df
               dfDayHours <- dfDay %>% 
                 dplyr::filter(rideHourNum %in% as.numeric(substr(entireDayConverter(h),
                                                                  start = 2, stop = 3)))
               
               return (
                 htmlDiv(children = list(
                   htmlP(paste("Total Rides: ",
                               as.character(nrow(dfDay)), sep = ""),
                         className = "totalsP"),
                   htmlP(paste("Total Selected Rides: ",
                               as.character(nrow(dfDayHours)), sep = ""),
                         className = "totalsP"),
                   htmlP(paste("Selected Intervals: ",
                               as.character(length(unique(dfDayHours$rideHourNum)))
                               , sep = ""),
                         className = "totalsP")
                 )
                 )
               )
             }
)

####################################################################################################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}