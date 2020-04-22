library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dplyr)
library(stringr)
library(glue)
library(data.table)
library(plotly)
library(tidyr)


appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

app <- Dash$new()

# get lat long data and transform
df_lat_lon <- fread("lat_lon_counties.csv")
df_lat_lon$FIPS <- str_pad(df_lat_lon$FIPS, 5, pad = "0")

# get age adjusted  data and transform
df_full_data <- fread("age_adjusted_death_rate_no_quotes.csv")
df_full_data$`County Code` <- str_pad(df_full_data$`County Code`, 5, pad = "0")
df_full_data$County <- paste0(df_full_data$County, ", ", df_full_data$State)
df_full_data <- df_full_data[, -c(1,3)]
df_full_data$`Age Adjusted Rate` <- as.numeric(gsub("[^0-9.]", "", df_full_data$`Age Adjusted Rate`))
df_full_data$Deaths <- as.numeric(gsub("[^0-9.]", "", df_full_data$Deaths))

YEARS = c(2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015)

BINS <- list("0-2", "2.1-4", "4.1-6", "6.1-8", "8.1-10", "10.1-12", "12.1-14",
             "14.1-16", "16.1-18", "18.1-20", "20.1-22", "22.1-24", "24.1-26",
             "26.1-28", "28.1-30", ">30")

DEFAULT_COLORSCALE <- list("#f2fffb", "#bbffeb", "#98ffe0", "#79ffd6", "#6df0c8",
                           "#69e7c0", "#59dab2", "#45d0a5", "#31c194", "#2bb489",
                           "#25a27b", "#1e906d", "#188463", "#157658", "#11684d",
                           "#10523e")

DEFAULT_OPACITY <- 0.8

mapbox_access_token <- "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNrOWJqb2F4djBnMjEzbG50amg0dnJieG4ifQ.Zme1-Uzoi75IaFbieBDl3A"
mapbox_style <- "mapbox://styles/plotlymapbox/cjvprkf3t1kns1cqjxuxmwixz"

#####################
##### APP LAYOUT ####
#####################

app$layout(htmlDiv(
  id ="root",
  children = list(
    htmlDiv(list(
      htmlH4(children = "Rate of US Poison-Induced Deaths"),

      htmlImg(
        src  =  "assets/dash-logo-new.png",
        style  = list('float' = "right", 'height' = "50px")
      )
    )),
    htmlP(id = "description",
          children = "† Deaths are classified using the International Classification of Diseases, \
          Tenth Revision (ICD–10). Drug-poisoning deaths are defined as having ICD–10 underlying \
          cause-of-death codes X40–X44 (unintentional), X60–X64 (suicide), X85 (homicide), or Y10–Y14 \
          (undetermined intent)."
    ),

    htmlDiv(
      id = "app-container",
      children = list(
        htmlDiv(
          id = "left-column",
          className = "seven columns",
          children = list(
            htmlDiv(
              id = "slider-container",
              children = list(
                htmlP(id = "slider-text",
                      children = "Drag the slider to change the year:"),
                dccSlider(
                  id = "years-slider",
                  min = min(YEARS),
                  max = max(YEARS),
                  value = min(YEARS),
                  marks <- as.list(
                    setNames(
                      as.character(seq(from = min(YEARS), to = max(YEARS), by = 1)),
                      as.character(seq(from = min(YEARS), to = max(YEARS), by = 1))
                    )
                  )
                )
              )
            ),
            htmlDiv(
              id="heatmap-container",
              children = list(
                htmlP(glue("Heatmap of age adjusted mortality rates from poisonings in year {min(YEARS)}"),
                      id = "heatmap-title"
                ),
                dccGraph(
                  id = "county-choropleth",
                  figure = list(
                    data = list(list(
                      lat = df_lat_lon[["Latitude"]],
                      lon = df_lat_lon[["Longitude"]],
                      text = df_lat_lon[["Hover"]],
                      type = "scattermapbox"
                    )),
                    layout = list(
                      mapbox = list(
                        layers = list(),
                        accesstoken = mapbox_access_token,
                        style = mapbox_style,
                        center = list(
                          lat = 38.72490,
                          lon = -95.61446
                        ),
                        pitch = 0,
                        zoom = 3.5
                      ),
                      autosize = TRUE
                    )
                  )
                )
              )
            )
          )),
        htmlDiv(
          id = "graph-container",
          children = list(
            htmlP(id="chart-selector",
                  children="Select chart:"
            ),
            dccDropdown(
              options = list(list("label"= "Histogram of total number of deaths (single year)",
                                  "value"= "show_absolute_deaths_single_year"),
                             list("label"= "Histogram of total number of deaths (1999-2016)",
                                  "value"= "absolute_deaths_all_time"),
                             list("label"= "Age-adjusted death rate (single year)",
                                  "value"= "show_death_rate_single_year"),
                             list("label"= "Trends in age-adjusted death rate (1999-2016)",
                                  "value"= "death_rate_all_time")),
              value="show_death_rate_single_year",
              id="chart-dropdown"
            ),
            dccGraph(
              id = "selected-data", #style = list('padding' = "30px"),
              figure = list(
                data = list(list(x = 0, y = 0)),
                layout = list(
                  paper_bgcolor = "#F4F4F8",
                  plot_bgcolor = "#F4F4F8",
                  autofill = TRUE

                )
              )
            )
          ), className = "five columns")
      ))
  )
))
app$callback(
  output = list("id" = "county-choropleth", "property" = "figure"),
  params = list(
    input("id" = "years-slider", "property" = "value"),
    state("id" = "county-choropleth", "property" = "figure")
  ),

  function(year, figure){
    cm <- mapply(c, BINS, DEFAULT_COLORSCALE, SIMPLIFY = FALSE)
    data <- list(list(
      lat = df_lat_lon[["Latitude"]],
      lon = df_lat_lon[["Longitude"]],
      text = df_lat_lon[["Hover"]],
      type = "scattermapbox",
      hoverinfo = "text",
      marker = list(size = 5, color = "white", opacity = 0)
    ))
    annotations <- list(list(
      showarrow = FALSE,
      align = "right",
      text = "<b>Age-adjusted death rate<br>per county per year</b>",
      font = list(color = "#2cfec1"),
      bgcolor = "#1f2630",
      x = 0.95,
      y = 0.95
    ))

    Bin_Rev <- rev(BINS)
    cm_name <- sapply(cm, function(n) n[1])
    annotations_1 <- lapply(1:length(Bin_Rev),
                            function(i) {
                              index <- which(cm_name %in% Bin_Rev[[i]])
                              color <- cm[[index]][2]
                              bin <- cm[[index]][1]
                              list(
                                arrowcolor = color,
                                text = bin,
                                x = 0.95,
                                y = 0.85 - ((i-1) / 20),
                                ax = -60,
                                ay = 0,
                                arrowwidth = 5,
                                arrowhead = 0,
                                bgcolor = "#1f2630",
                                font = list(color = "#2cfec1")
                              )
                            })

    annotations <- append(annotations, annotations_1)

    if("layout" %in% figure){
      lat = figure[["layout"]][["mapbox"]][["center"]]["lat"]
      lon = figure[["layout"]][["mapbox"]][["center"]]["lon"]
      zoom = figure[["layout"]][["mapbox"]]["zoom"]
    } else{
      lat = 38.72490
      lon = -95.61446
      zoom = 3.5
    }

    layout <- list(
      mapbox = list(
        layers = list(),
        accesstoken = mapbox_access_token,
        style = mapbox_style,
        center = list(lat = lat, lon = lon),
        zoom = zoom
      ),
      hovermode = "closet",
      margin = list(r = 0, l = 0, t = 0, b = 0),
      annotations = annotations,
      dragmode = "lasso"
    )

    base_url = "https://raw.githubusercontent.com/jackparmer/mapbox-counties/master/"

    cm_name <- sapply(cm, function(n) n[1])
    geo_layer <- lapply(1:length(BINS),
                        function(i){
                          index <- which(cm_name %in% BINS[[i]])
                          color <- cm[[index]][2]
                          list(
                            sourcetype = "geojson",
                            source = paste0(base_url, as.character(year), "/", BINS[[i]], ".geojson"),
                            type = "fill",
                            color = color,
                            opacity = DEFAULT_OPACITY,
                            fill = list(outlinecolor = "#afafaf")
                          )
                        })
    layout$mapbox$layers = geo_layer
    fig <- list("data" = data, "layout" = layout)
    return(fig)
  }
)

app$callback(
  output = list("id" = "heatmap-title", property = "children"),
  params = list(input("id" = "years-slider", property = "value")),

  function(year){
    # Function prints year selected from range slider
    return(glue("Heatmap of age adjusted mortality rates from poisonings in year {year}"))
  }
)

app$callback(
  output = list(id = "selected-data", property = "figure"),
  params = list(input(id = "county-choropleth", property = "selectedData"),
                input(id = "chart-dropdown", property = "value"),
                input(id = "years-slider", property = "value")),

  function(selectedData, chart_dropdown, year){
    # function returns blank figure if no data is selected
    if(is.null(selectedData[[1]])){
      return(
        list(
          data = list(list(x = 0, y = 0)),
          layout = list(
            title = "Click-drag on the map to select counties",
            paper_bgcolor = "#1f2630",
            plot_bgcolor = "#1f2630",
            font = list(
              color="#2cfec1"
            )
          )
        )
      )
    }else{
      pts <- selectedData[["points"]]
      fips <- lapply(1:length(pts),
                     function(i){
                       gsub("[^0-9.]", "", pts[[i]]["text"])
                     })

      fips <- sapply(1:length(fips),
                     function(i){

                       if(is.na(fips[[i]])){
                         "0" #fips <- fips[-fips[[i]]]
                       }else{
                         if(nchar(fips[[i]]) == 4){
                           fips[[i]] = paste0("0", fips[[i]])
                         }else{fips[[i]] = fips[[i]]}
                       }
                     })

      fips <- fips[fips %in% "0" == FALSE]
      dff <- filter(df_full_data, df_full_data$`County Code` %in% fips)
      dff <- dff[with(dff, order(Year)), ]
      if(chart_dropdown != "death_rate_all_time"){
        title <- "Absolute deaths per county, <b>1999-2016</b>"
        AGGREGATE_BY <- "Deaths"
        if("show_absolute_deaths_single_year" == chart_dropdown){
          dff <- filter(dff, dff$Year == year) # need to check year
          title <- glue("Absolute deaths per county, <b>{year}</b>")
        }else if("show_death_rate_single_year" == chart_dropdown){
          dff <- filter(dff, dff$Year == year)
          title <- glue("Absolute deaths per county, <b>{year}</b>")
          AGGREGATE_BY <- "Age Adjusted Rate"
        }

        deaths_or_rate_by_fips <- aggregate(x = dff[[AGGREGATE_BY]], by = list(County=dff$County),
                                            FUN = function(x){sum((x[!is.na(x)]))})

        # Only look at non-zero rows:
        deaths_or_rate_by_fips <- filter(deaths_or_rate_by_fips, deaths_or_rate_by_fips[[2]] > 0)
        deaths_or_rate_by_fips <- deaths_or_rate_by_fips[with(deaths_or_rate_by_fips, order(x)), ]
        deaths_or_rate_by_fips <- filter(deaths_or_rate_by_fips, !is.null(deaths_or_rate_by_fips[2]))

        # plot for histograms (Death by year)
        x <- deaths_or_rate_by_fips[,1]
        y <- deaths_or_rate_by_fips[,2]

        fig <- plot_ly(x = x, y = y, type = "bar",
                       marker = list(color = "#2cfec1", size = 4, line = list(width = 0)),
                       textposition = "outside", opacity = 1) %>%
          layout(title = title,
                 autosize = TRUE,
                 paper_bgcolor = "#1f2630",
                 plot_bgcolor = "#1f2630",
                 font = list(color = "#2cfec1"),
                 xaxis = list(title = "",
                              tickfont = list(color = "#2cfec1"),
                              gridcolor = "#5b5b5b",
                              categoryorder = "array",
                              categoryarray = c(sort(deaths_or_rate_by_fips[,2]))
                 ),
                 yaxis = list(title = "",
                              tickfont = list(color = "#2cfec1"),
                              gridcolor = "#5b5b5b"),
                 margin = list(
                   r = 30,
                   t = 40,
                   b = 0,
                   l = 30
                 )
          )
        return(fig)
      }

      fig <- plot_ly(data = dff , x = ~Year, y = ~Deaths, color = ~County, type = "scatter",
                     colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
                                "#e6ab02", "#a6761d", "#666666", "#1b9e77"),
                     textposition = "outside", opacity = 1,  mode = "markers+lines",
                     marker = list(size = 4, line = list(width = 1))) %>%
        layout(title = list(glue("<b>{length(fips)} counties selected")),
               paper_bgcolor = "#1f2630",
               plot_bgcolor = "#1f2630",
               font = list(color = "#2cfec1"),
               autosize = TRUE,
               hovermode = "closest",
               xaxis = list(title = "",
                            tickfont = list(color = "#2cfec1"),
                            gridcolor = "#5b5b5b",
                            fixedrange = FALSE),
               yaxis = list(title = "Age-adjusted death rate per county per year",
                            tickfont = list(color = "#2cfec1"),
                            gridcolor = "#5b5b5b",
                            fixedrange = TRUE),
               margin = list(
                 r = 30,
                 t = 40,
                 b = 0,
                 l = 30
               )
        )
      return(fig)
    }
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}



