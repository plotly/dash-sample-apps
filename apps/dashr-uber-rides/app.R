appName <- Sys.getenv("DASH_APP_NAME")
pathPrefix <- sprintf("/%s/", appName)

Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
           DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

setwd("/app/apps/dashr-uber-rides")

library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(data.table)
library(dplyr)
library(lubridate)
library(Hmisc)

#################################### LOAD DATA & CREATE GLOBAL OBJECTS #############################
ridesRaw <- read.csv("data/output.csv",
                     stringsAsFactors = FALSE)
# Read data

ridesDf <- ridesRaw[, c("Date.Time","Lat","Lon")]
# Remove extra columns

ridesDf$Date.Time <- as.POSIXct(ridesDf$Date.Time, format = "%Y-%m-%d %H:%M:%OS")
# Convert Date.Time column to posix type || THIS IS ABOUT A MINUTE WITH THE ACTUAL DF 

locationCoordinates <- list(
  list(40.7272, -73.991251), list(40.7505, -73.9934),
  list(40.8296, -73.9262), list(40.7484, -73.9857),
  list(40.7069, -74.0113), list(40.644987, -73.785607),
  list(40.7527, -73.9772), list(40.7589, -73.9851))
names(locationCoordinates) <-
             c("default", "Madison Square Garden", "Yankee Stadium",
               "Empire State Building", "New York Stock Exchange",
               "JFK Airport", "Grand Central Station", "Times Square")
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
    c(0, "#F4EC15"),
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
    list(csSeq[23], "#603099"))
                
# Create custom colorScale list 

ridesDf$rideHour <- as.numeric(as.POSIXlt(ridesDf$Date.Time)$hour)
# Extract hour from dates and append corresponding color codes as a column for coloring later
# Factor works too

mapboxToken <- paste("pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNqdnBvNDMyaT",
                "AxYzkzeW5ubWdpZ2VjbmMifQ.TXcBE-xg9BFdV2ocecc_7g", sep = "")
#mapbox token for plotly | this one is for plot_mapbox figure

Sys.setenv("MAPBOX_TOKEN" = mapboxToken)
#setting mapbox token for R environment 

locationOptions <- list(
  list(label = "New York Main", value = "default"),
  list(label = "Madison Square Garden", value = "Madison Square Garden"),
  list(label = "Yankee Stadium", value = "Yankee Stadium"),
  list(label = "Empire State Building", value = "Empire State Building"),
  list(label = "New York Stock Exchange", value = "New York Stock Exchange"),
  list(label = "JFK Airport", value = "JFK Airport"),
  list(label = "Grand Central Station", value = "Grand Central Station"),
  list(label = "Times Square", value = "Times Square"))
                  
allHours <- c(" 00:00:00", " 01:00:00", " 02:00:00", " 03:00:00",
               " 04:00:00", " 05:00:00", " 06:00:00", " 07:00:00",
               " 08:00:00", " 09:00:00", " 10:00:00", " 11:00:00",
               " 12:00:00", " 13:00:00", " 14:00:00", " 15:00:00",
               " 16:00:00", " 17:00:00", " 18:00:00", " 19:00:00",
               " 20:00:00", " 21:00:00", " 22:00:00", " 23:00:00")
               

hourOptions <- list(
                  list(label = "Entire Day", value = "entire_day"),
                  list(label = "00:00", value = " 00:00:00"),
                  list(label = "01:00", value = " 01:00:00"),
                  list(label = "02:00", value = " 02:00:00"),
                  list(label = "03:00", value = " 03:00:00"),
                  list(label = "04:00", value = " 04:00:00"),
                  list(label = "05:00", value = " 05:00:00"),
                  list(label = "06:00", value = " 06:00:00"),
                  list(label = "07:00", value = " 07:00:00"),
                  list(label = "08:00", value = " 08:00:00"),
                  list(label = "09:00", value = " 09:00:00"),
                  list(label = "10:00", value = " 10:00:00"),
                  list(label = "11:00", value = " 11:00:00"),
                  list(label = "12:00", value = " 12:00:00"),
                  list(label = "13:00", value = " 13:00:00"),
                  list(label = "14:00", value = " 14:00:00"),
                  list(label = "15:00", value = " 15:00:00"),
                  list(label = "16:00", value = " 16:00:00"),
                  list(label = "17:00", value = " 17:00:00"),
                  list(label = "18:00", value = " 18:00:00"),
                  list(label = "19:00", value = " 19:00:00"),
                  list(label = "20:00", value = " 20:00:00"),
                  list(label = "21:00", value = " 21:00:00"),
                  list(label = "22:00", value = " 22:00:00"),
                  list(label = "23:00", value = " 23:00:00"))

externalSheets <- list("https://codepen.io/chriddyp/pen/bWLwgP.css")

####################################################################################################

#################################### APP START #####################################################

app <- Dash$new(name = "dashr-uber-rides", external_stylesheets = externalSheets)
# Initiate application

#, meta_tags=list(list("name"= "viewport", "content"= "width=device-width, initial-scale=1")))
# meta_tags=[
#   {"name": "viewport", "content": "width=device-width, initial-scale=1"}
#   ]
  
####################################################################################################

#################################### CREATE LAYOUT VARIABLES #######################################

pageTitle <- htmlH2("DASH - UBER DATA APP", className = "title-left")
firstP <- htmlP(paste("Select different dates, locations or ",
                "hours using dropdowns below or by selecting ",
                "different time frames on the histogram", sep = ""))
plotlyLogo <- htmlA(
                  list(
                    htmlImg(src = "assets/dashR-logo-stripe.png",
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
  date = make_date(year = 2014L, month = 4L, day = 1L),
  display_format = "MMMM Y, DD",
  min_date_allowed = make_date(year = 2014L, month = 4L, day = 1L),
  max_date_allowed = make_date(year = 2014L, month = 9L, day = 30L),
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
                ), className = "three columns div-user-controls"
        ),

        htmlDiv(list(
                  dccGraph(id = "map-graph", className = "map-graph"),
                  dccGraph(id = "histogram", className = "histogram",
                     selectedData = list(points = list(), range = NULL),
                     clickData = list(points = list(), range = NULL))
                ), className = "div-for-charts nine columns"
        )
      )
    )
)

####################################################################################################

#################################### HELPER FUNCS FOR CALLBACKS ####################################

#Filters out the date selected by user from dataset
FilterDay <- function(date){

  datePaste <- paste(date, " 00:00:00", sep = "")
  #Combine month-day-year to %Y-%m-%d %H:%M:%OS format

  dayCurrent <- as.POSIXct(datePaste)
  dayNext <- as.POSIXct(datePaste) + 86400
  #Created today and tomorrows posix object for filtering

  dfDay <- ridesDf[ridesDf$Date.Time >= dayCurrent &
                       ridesDf$Date.Time < dayNext, ]
  #Filter out data

  return(dfDay %>% arrange(`Date.Time`))
}

# returns hour list for callbacks
EntireDayConverter <- function(hour){

  hourVec <- unlist(hour)

  if ("entire_day" %in% hourVec) {
    hourVec <- allHours
    # assign allHours vector to hour if "entire_day" selected
  }
  return(hourVec)
}

# Filters out the month-day-hour selected by user from dataset
DfToMap <- function(date, hour){

  hourVec <- EntireDayConverter(hour)

  if (length(hourVec) >= 24){
    return(FilterDay(date))
  # Avoid loop if entire day selected for faster filtering
  }else{

    dfEmpty <- ridesDf[0, ]

    for (h in hourVec){

      hourPaste <- paste(date, h, sep = "")
      #Combine month-day-year to %Y-%m-%d %H:%M:%OS format

      hourStart <- as.POSIXct(hourPaste)
      hourEnd <- hourStart + 3600
      #Define df slicing interval

      ridesHour <- ridesDf[ridesDf$Date.Time >= hourStart &
                           ridesDf$Date.Time < hourEnd, ]
      #Filter out data

      dfEmpty <- rbind(dfEmpty, ridesHour)
      #Append hours to get df to plot
    }
    return(dfEmpty)
  }
}

####################################################################################################

#################################### CALLBACKS START ###############################################

app$callback(output = list(id = "hidden-date-picker", property = "children"),
             params = list(input(id = "date-picker", property = "date")),

  function(date){
    return(as.character(date))
  }
)

app$callback(output = list(id = "map-graph", property = "figure"),
             params = list(
                       input(id = "hidden-date-picker", property = "children"),
                       input(id = "hour-dropdown", property = "value"),
                       input(id = "location-dropdown", property = "value"),
                       state(id = "map-graph", property = "relayoutData")),

  function(date, hour, location, relayout){

    dfToMap <- DfToMap(date, hour)

      subsetIndex <- which(names(locationCoordinates) %in% location)

      latInitial <- locationCoordinates[[subsetIndex]][[1]]
      lonInitial <- locationCoordinates[[subsetIndex]][[2]]
      zoom <- 15
      bearing <- 0
      # Set layout variables for locations

      if (location == "default"){
      zoom <- 12
      }
      return(
        plot_mapbox(dfToMap, lat = ~Lat, lon = ~Lon, split = ~rideHour,
                     mode = "markers", hoverinfo = "lat+lon+text",
                     marker = list(color = ~rideHour, cmin = 0, cmax = 23,
                                   colorscale = colorScale, size = 4)
        ) %>% layout( title = "Uber Rides in NYC",
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

  function(date, hour){

    dfDay <- FilterDay(date)
    #Filter out data

    dfAgg <- aggregate(Date.Time ~
              format(dfDay$Date.Time, "%Y-%m-%d %H"),
              data = dfDay, length)
    # aggregate number of hourly rides to df

    names(dfAgg) <- c("Date-Hours", "Ride_Counts")
    # name aggregate df appropriately

    dfAgg$`Date-Hours` <- seq(0, 23)
    #reformatting for getting the x-axis ticks correct

    hourVec <- EntireDayConverter(hour)

    colorMap <- colorMapOrj

    numericHours <- as.numeric(substr(hourVec, start = 2, stop = 3))

    if (length(numericHours) > 0 && length(numericHours) < 24){
      colorMap[numericHours + 1] <- "FFFFFF"
    }
    return(
      plot_ly(
               x = dfAgg$`Date-Hours`,
               y = dfAgg$Ride_Counts,
               name = "histogram",
               marker = list(
                   color = colorMap),
               type = "bar") %>% layout(
                  title = list(
                    text = paste("Select any of the bars",
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
                    add_annotations(text = dfAgg$Ride_Counts,
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

  function(select, click){

    dfSelectedHour <- rbindlist(select$points)
    # Convert select list [points] to df

    if (nrow(dfSelectedHour) > 0){
      dfClick <- dfSelectedHour[0, ]
     # Reset dfClick when selectedData
    }else{
      dfClick <- rbindlist(click$points)
     # Add clicked data to df
    }

   dfSelectedHour <- rbind(dfSelectedHour, dfClick)
   # Append selection and click df"s together

   if (nrow(dfSelectedHour) > 0){
        dfSelectedHour <- dfSelectedHour %>% mutate(
          hoursToReturn = ifelse(pointNumber < 10,
              paste(" 0", pointNumber, ":00:00", sep = ""),
              paste(" ", pointNumber, ":00:00", sep = "")))
   # Create a new column with the right hour formatting
   }else{
     dfSelectedHour$hoursToReturn <- "entire_day"
   }
   return(dfSelectedHour$hoursToReturn)
  }
)

app$callback(output = list(id = "total-rides", property = "children"),
             params = list(
                 input(id = "hidden-date-picker", property = "children"),
                 input(id = "hour-dropdown", property = "value")),

  function(date, h){

    dfDay <- FilterDay(date)

    dfDay <- dfDay %>% mutate(
      rideHourNum = as.numeric(as.POSIXlt(dfDay$Date.Time)$hour))
    # Created rideHourNum

    hoursToFilter <- as.numeric(substr(EntireDayConverter(h),
                                         start = 2, stop = 3))
    # Format hours to filter

    dfDayHours <- dfDay[dfDay$rideHourNum %in% hoursToFilter, ]
    # Filter df

    return(
      htmlDiv(children = list(
      htmlP(paste("# Total Rides: ",
                  as.character(nrow(dfDay)), sep = ""),
            className = "totalsP"),
      htmlP(paste("# Total Selected Rides: ",
                  as.character(nrow(dfDayHours)), sep = ""),
            className = "totalsP"),
      htmlP(paste("# Selected Intervals: ",
                  as.character(length(unique(dfDayHours$rideHourNum)))
                  , sep = ""),
            className = "totalsP")
              )
      )
    )
  }
)
####################################################################################################



app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))

