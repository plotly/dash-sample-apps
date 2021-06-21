library(dash)
library(dashHtmlComponents)
library(dashCoreComponents)
library(rapportools)
library(devtools)
library(gsubfn)
library(lubridate)
library(data.table)
library(rapportools)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

app <- Dash$new()
source("points.R") #Ensure you are able to source these and run.
source("controls.R")

# Create global chart template
mapbox_access_token = 'pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNrOWJqb2F4djBnMjEzbG50amg0dnJieG4ifQ.Zme1-Uzoi75IaFbieBDl3A'
county_options <- lapply(1:length(COUNTIES),
                         function(x){list(label=unname(COUNTIES[x]),
                                          value=labels(COUNTIES[x]))})

well_status_options = lapply(1:length(WELL_STATUSES),
                             function(x){list(label=unname(WELL_STATUSES[x]),
                                              value=labels(WELL_STATUSES[x]))})

well_type_options = lapply(1:length(WELL_TYPES),
                           function(x){list(label=unname(WELL_TYPES[x]),
                                            value=labels(WELL_TYPES[x]))})

# Load data
df <- read.csv('data/wellspublic.csv.gz')
df$Date_Well_Completed <- as.Date(df$Date_Well_Completed, format = "%Y-%m-%d %H:%M:%S")
df <- df[df$Date_Well_Completed > as.Date('1960-01-01'),]
rownames(df) <- 1:nrow(df)
trim <- na.omit(df[,c('API_WellNo', 'Well_Type', 'Well_Name')])

rownames(trim) <- seq.int(nrow(trim))

app$layout(htmlDiv(
  list(
    dccStore(id="aggregate_data"),
    htmlDiv(
      list(
        htmlDiv(
          list(
            htmlImg(
              src="/assets/dash-logo.png",
              id="plotly-image",
              style=list(
                "height"= "60px",
                "width"= "auto",
                "margin-bottom"= "25px"
              )
            )
          ),
          className="one-third column"
        ),
        htmlDiv(
          list(
            htmlDiv(
              list(
                htmlH3(
                  "New York Oil and Gas",
                  style=list("margin-bottom"= "0px")
                ),
                htmlH5(
                  "Production Overview", style=list("margin-top"= "0px")
                )
              )
            )
          ),
          className="one-half column",
          id="title"
        ),
        htmlDiv(
          list(
            htmlA(
              htmlButton("Learn More", id="learn-more-button"),
              href="https=//plot.ly/dash/pricing/"
            )
          ),
          className="one-third column",
          id="button"
        )
      ),
      id="header",
      className="row flex-display",
      style=list("margin-bottom"= "25px")
    ),
    htmlDiv(
      list(
        htmlDiv(
          list(
            htmlP(
              "Filter by construction date (or select range in histogram):",
              className="control_label"
            ),
            dccRangeSlider(
              id="year_slider",
              min=1960,
              max=2017,
              value=list(1990, 2010),
              className="dcc_control"
            ),
            htmlP("Filter by well status=", className="control_label"),
            dccRadioItems(
              id="well_status_selector",
              options=list(
                list("label"= "All ", "value"= "all"),
                list("label"= "Active only ", "value"= "active"),
                list("label"= "Customize ", "value"= "custom")
              ),
              value="active",
              labelStyle=list("display"= "inline-block"),
              className="dcc_control"
            ),
            dccDropdown(
              id="well_statuses",
              options=well_status_options,
              multi=TRUE,
              value=list(labels(WELL_STATUSES)),
              className="dcc_control"
            ),
            dccChecklist(
              id="lock_selector",
              options=list(list("label"= "Lock camera", "value"= "locked")),
              value=list(),
              className="dcc_control"
            ),
            htmlP("Filter by well type=", className="control_label"),
            dccRadioItems(
              id="well_type_selector",
              options=list(
                list("label"= "All ", "value"= "all"),
                list("label"= "Productive only ", "value"= "productive"),
                list("label"= "Customize ", "value"= "custom")
              ),
              value="productive",
              labelStyle=list("display"= "inline-block"),
              className="dcc_control"
            ),
            dccDropdown(
              id="well_types",
              options=well_type_options,
              multi=TRUE,
              value=list(labels(WELL_TYPES)),
              className="dcc_control"
            )
          ),
          className="pretty_container four columns",
          id="cross-filter-options"
        ),
        htmlDiv(
          list(
            htmlDiv(
              list(
                htmlDiv(
                  list(htmlH6(id="well_text"), htmlP("No. of Wells")),
                  id="wells",
                  className="mini_container"
                ),
                htmlDiv(
                  list(htmlH6(id="gasText"), htmlP("Gas")),
                  id="gas",
                  className="mini_container"
                ),
                htmlDiv(
                  list(htmlH6(id="oilText"), htmlP("Oil")),
                  id="oil",
                  className="mini_container"
                ),
                htmlDiv(
                  list(htmlH6(id="waterText"), htmlP("Water")),
                  id="water",
                  className="mini_container"
                )
              ),
              id="info-container",
              className="row container-display"
            ),
            htmlDiv(
              list(dccGraph(id="count_graph")),
              id="countGraphContainer",
              className="pretty_container"
            )
          ),
          id="right-column",
          className="eight columns"
        )
      ),
      className="row flex-display"
    ),
    htmlDiv(
      list(
        htmlDiv(
          list(dccGraph(id="main_graph")),
          className="pretty_container seven columns"
        ),
        htmlDiv(
          list(dccGraph(id="individual_graph")),
          className="pretty_container five columns"
        )
      ),
      className="row flex-display"
    ),
    htmlDiv(
      list(
        htmlDiv(
          list(dccGraph(id="pie_graph")),
          className="pretty_container seven columns"
        ),
        htmlDiv(
          list(dccGraph(id="aggregate_graph")),
          className="pretty_container five columns"
        )
      ),
      className="row flex-display"
    )
  ),
  id="mainContainer",
  style=list("display"= "flex", "flex-direction"= "column")
))

filter_Dataframe <- function(df,
                             well_statuses,
                             well_types,
                             year_slider){
  df <- as.data.table(df)
  min_year <- as.Date(paste0(year_slider[[1]], "-01-01"))
  max_year <- as.Date(paste0(year_slider[[2]], "-01-01"))
  df <- df[Well_Status %in% well_statuses & Well_Type %in% well_types,]
  df <- df[Date_Well_Completed > min_year & Date_Well_Completed < max_year, ]
  return(df)
}

fetch_individual <- function(api){

  index <- as.list(c(min(dfIn[dfIn$API_WellNo == api,'Reporting.Year']):
                       max(dfIn[dfIn$API_WellNo == api,'Reporting.Year'])))
  gas <- list()
  oil <- list()
  water <- list()

  for (year in index) {

    gas[[year]]  <- dfIn[dfIn$API_WellNo == api & dfIn$Reporting.Year == year,'Gas.Produced..MCF']
    oil[[year]]  <- dfIn[dfIn$API_WellNo == api & dfIn$Reporting.Year == year,'Oil.Produced..bbl']
    water[[year]] <- dfIn[dfIn$API_WellNo == api & dfIn$Reporting.Year == year,'Water.Produced..bbl']
  }
  return(list(index,water[unlist(index)], oil[unlist(index)], gas[unlist(index)]))
}

fetch_aggregate <- function(dff,year_slider){
  gas1 = 0
  oil1 = 0
  water1 = 0
  min_year <- paste0(year_slider[[1]], "-01-01")

  index = max(year(as.Date(min_year)):1985):2017

  gas <- list()
  oil <- list()
  water <- list()
  inn <- merge(dff,dfIn, by = 'API_WellNo')
  for (year in index) {
    count_gas = 0
    count_oil = 0
    count_water = 0
    count_gas <- sum(inn$Gas.Produced..MCF[inn$Reporting.Year %in% year])
    count_water <- sum(inn$Water.Produced..bbl[inn$Reporting.Year %in% year])
    count_oil <- sum(inn$Oil.Produced..bbl[inn$Reporting.Year %in% year])
    gas[[year]] = count_gas
    water[[year]] = count_water
    oil[[year]] = count_oil
  }
  gas1 <- Reduce("+", gas[index])
  oil1 <- Reduce("+", oil[index])
  water1 <- Reduce("+", water[index])
  return(list(gas1,oil1,water1))
}

fetch_abs <- function(selected, year_slider){
  gas1 = 0
  oil1 = 0
  water1 = 0
  min_year <- paste0(year_slider[[1]], "-01-01")

  index = max(year(as.Date(min_year)):1985):2017

  gas = list()
  oil = list()
  water = list()
  inn <- merge(selected,dfIn, by = 'API_WellNo')
  for (year in index) {
    count_gas = 0
    count_oil = 0
    count_water = 0
    count_gas <- sum(inn$Gas.Produced..MCF[inn$Reporting.Year %in% year])
    count_water <- sum(inn$Water.Produced..bbl[inn$Reporting.Year %in% year])
    count_oil <- sum(inn$Oil.Produced..bbl[inn$Reporting.Year %in% year])
    gas[[year]] = count_gas
    water[[year]] = count_water
    oil[[year]] = count_oil
  }
  return(list(index,gas,oil,water))
}
# Radio -> multi
app$callback(output=list(id='well_statuses', property='value'),
             params=list(
               input(id='well_status_selector', property='value')),
             function(selector){
               if (selector == 'all') {
                 return(labels(WELL_STATUSES))
               } else if (selector == 'active') {
                 return(list('AC'))
               } else
                 return(list())
             })

app$callback(output=list(id='well_types', property='value'),
             params=list(
               input(id='well_type_selector', property='value')),
             function(selector){
               if (selector == 'all') {
                 return(labels(WELL_TYPES))
               } else if (selector == 'productive') {
                 return(list('GD', 'GE', 'GW', 'IG', 'IW', 'OD', 'OE', 'OW'))
               } else
                 return(list())
             })
app$callback(output=list(id='well_text', property='children'),
             params=list(
               input(id='well_statuses', property='value'),
               input(id='well_types', property='value'),
               input(id='year_slider', property='value')),
             function(well_statuses, well_types, year_slider){
               dff <- filter_Dataframe(df, well_statuses, well_types, year_slider)
               sprintf('%g ⑆', nrow(dff))})


app$callback(output = list(id ='gasText', property ='children'),
             params=list(
               input(id='well_statuses', property='value'),
               input(id='well_types', property='value'),
               input(id='year_slider', property='value')),
             function(well_statuses, well_types, year_slider){
               dff <- filter_Dataframe(df, well_statuses, well_types, year_slider)
               sprintf("%g M mcg ⑆", round(as.numeric(fetch_aggregate(dff,year_slider)[1])/1000000,2))
             })
app$callback(output = list(id ='oilText', property ='children'),
             params=list(
               input(id='well_statuses', property='value'),
               input(id='well_types', property='value'),
               input(id='year_slider', property='value')),
             function(well_statuses, well_types, year_slider){
               dff <- filter_Dataframe(df, well_statuses, well_types, year_slider)
               sprintf("%g M bbl ⑆", round(as.numeric(fetch_aggregate(dff,year_slider)[2])/1000000,2))
             })
app$callback(output = list(id ='waterText', property ='children'),
             params=list(
               input(id='well_statuses', property='value'),
               input(id='well_types', property='value'),
               input(id='year_slider', property='value')),
             function(well_statuses, well_types, year_slider){
               dff <- filter_Dataframe(df, well_statuses, well_types, year_slider)
               sprintf("%g M bbl ⑆", round(as.numeric(fetch_aggregate(dff,year_slider)[3])/1000000,2))
             })

app$callback(output=list(id='year_text', property='children'),
             params=list(
               input(id='year_slider', property='value')),
             function(year_slider){
               sprintf("%g ➮ %g", year_slider[1], year_slider[2])
             })

# Slider -> count graph
app$callback(output=list(id='year_slider', property='value'),
             params=list(
               input(id='count_graph', property='selectedData')),
             function(selectedData){
               if (length(selectedData[['points']]) == 0){
                 return(list(1990,2010))
               } else{
                 nums = list()
                 nums = unlist(lapply(selectedData[["points"]], "[[", "pointNumber"))
                 return(list(min(unlist(nums)) + 1960, max(unlist(nums)) + 1961))
               }
             }
)

#MAIN
app$callback(output=list(id='main_graph', property='figure'),
             params=list(
               input(id='well_statuses', property='value'),
               input(id='well_types', property='value'),
               input(id='year_slider', property='value'),
               state('lock_selector', 'values'),
               state('main_graph', 'relayoutData')),
             function(well_statuses, well_types, year_slider,selector, main_graph_layout){
               filtered_dff <- filter_Dataframe(df,well_statuses, well_types, year_slider)
               traces <- list()
               type <- unique(df$Well_Type)
               for (i in 1:length(type)){
                 df_by_type <- split(filtered_dff, as.factor(filtered_dff$Well_Type==type[i]))$`TRUE`
                 traces[[i]] <- list(
                   type='scattermapbox',
                   lon = df_by_type$Surface_Longitude,
                   lat = df_by_type$Surface_latitude,
                   opacity=0.6,
                   customdata = df_by_type$API_WellNo,
                   text = df_by_type$Well_Name,
                   mode = 'markers',
                   marker = list(
                     'size'= 4,
                     'line' = list('width' = 0.5, 'color' = 'black')
                   ),
                   name = as.character(WELL_TYPES[type[i]])
                 )}

               if (!is.null(main_graph_layout) & 'locked' %in% selector){
                 lon = -78.05
                 lat = 42.54
                 zoom = 7
               } else{
                 lon = -78.05
                 lat = 42.54
                 zoom = 7
               }

               return (list(
                 'data' = traces,
                 'layout'= list(
                   autosize=TRUE,
                   height=500,
                   font=list(color='#777777'),
                   titlefont=list(color='#777777', size='14'),
                   margin=list(
                     'l'=35,
                     'r'=35,
                     'b'=35,
                     't'=45
                   ),
                   hovermode="closest",
                   plot_bgcolor="#F9F9F9",
                   paper_bgcolor="#F9F9F9",
                   legend=list(font=list(size=10), orientation='h'),
                   title='Satellite Overview',
                   mapbox=list(
                     accesstoken=mapbox_access_token,
                     style="light",
                     center=list(
                       lon=-78.05,
                       lat=42.54
                     ),
                     zoom=7
                   ))
               ))
             })

#INDIVIDUAL
app$callback(output=list(id='individual_graph', property='figure'),
             params=list(
               input(id='main_graph', property='hoverData')),
             function(main_graph_hover){
               if (length(main_graph_hover[['points']]['customdata']) == 0){
                 #chosen <- list(31003002470000) #Default value,
                 index <- list(1995:2012)
                 water <- list(0,0,0,260,192,0,0,0,0,93,0,1342,0,26829,8952,17,0)
                 oil <- list(replicate(18,0))
                 gas <- list(replicate(18,0))
               } else{

                 chosen <- unlist(lapply(main_graph_hover[["points"]], "[[", "customdata"))
                 index <- fetch_individual(chosen[[1]])[1] #fetch ind is slow.
                 gas <- fetch_individual(chosen[[1]])[2]
                 oil <- fetch_individual(chosen[[1]])[3]
                 water <- fetch_individual(chosen[[1]])[4]
               }

               if (is.null(index)){
                 annotation = list(
                   text='No data available',
                   x=0.5,
                   y=0.5,
                   align="center",
                   showarrow=FALSE,
                   xref="paper",
                   yref="paper"
                 )
                 data = list()
               }

               else{
                 data = list(
                   list(
                     type='scatter',
                     mode='lines+markers',
                     name='Gas Produced (mcf)',
                     x=c(unlist(index)),
                     y=c(unlist(gas)),
                     line=list(
                       shape="spline",
                       smoothing=2,
                       width=1,
                       color='#fac1b7'
                     ),
                     marker=list(symbol='diamond-open')
                   ),
                   list(
                     type='scatter',
                     mode='lines+markers',
                     name='Oil Produced (bbl)',
                     x=c(unlist(index)),
                     y=c(unlist(oil)),
                     line=list(
                       shape="spline",
                       smoothing=2,
                       width=1,
                       color='#a9bb95'
                     ),
                     marker=list(symbol='diamond-open')
                   ),
                   list(
                     type='scatter',
                     mode='lines+markers',
                     name='Water Produced (bbl)',
                     x=c(unlist(index)),
                     y=c(unlist(water)),
                     line=list(
                       shape="spline",
                       smoothing=2,
                       width=1,
                       color='#92d8d8'
                     ),
                     marker=list(symbol='diamond-open')
                   ))
               }

               return(list(
                 'data' = data,
                 'layout' = list(
                   autosize=TRUE,
                   height=500,
                   font=list(color='#777777'),
                   titlefont=list(color='#777777', size='14'),
                   margin=list(
                     'l'=35,
                     'r'=35,
                     'b'=35,
                     't'=45
                   ),
                   hovermode="closest",
                   plot_bgcolor="#F9F9F9",
                   paper_bgcolor="#F9F9F9",
                   legend=list(font=list(size=10), orientation='h'),
                   title=sprintf('Individual Production: %s' ,unlist(lapply(main_graph_hover[["points"]], "[[", "text"))[1]),
                   mapbox=list(
                     accesstoken=mapbox_access_token,
                     style="light",
                     center=list(
                       lon=-78.05,
                       lat=42.54
                     ),
                     zoom=7
                   ))
               ))})

#COUNT
app$callback(output=list(id='count_graph', property='figure'),
             params=list(
               input(id='well_statuses', property='value'),
               input(id='well_types', property='value'),
               input(id='year_slider', property='value')),
             function(well_statuses, well_types, year_slider){
               year_slider = list(1960,2017)
               dff = filter_Dataframe(df, well_statuses, well_types, year_slider)
               g = dff[,c('API_WellNo', 'Date_Well_Completed')]
               gyeartotals <- data.frame(table(year(g$'Date_Well_Completed')))

               data = list(
                 list(
                   type='bar',
                   x=as.list(gyeartotals$Var1),
                   y=as.list(gyeartotals$Freq),
                   name='All Wells',
                   marker=list(
                     color='#87CEEB'
                   )
                 )
               )
               return(list(
                 'data' = data,
                 'layout' = list(
                   autosize=TRUE,
                   height=500,
                   font=list(color='#777777'),
                   titlefont=list(color='#777777', size='14'),
                   margin=list(
                     'l'=35,
                     'r'=35,
                     'b'=35,
                     't'=45
                   ),
                   hovermode="closest",
                   plot_bgcolor="#F9F9F9",
                   paper_bgcolor="#F9F9F9",
                   legend=list(font=list(size=10), orientation='h'),
                   title='Completed Wells / Year',
                   mapbox=list(
                     accesstoken=mapbox_access_token,
                     style="light",
                     center=list(
                       lon=-78.05,
                       lat=42.54
                     ),
                     zoom=7
                   ))
               ))


             })

#PIE
app$callback(output=list(id='pie_graph', property='figure'),
             params=list(
               input(id='well_statuses', property='value'),
               input(id='well_types', property='value'),
               input(id='year_slider', property='value')),

             function(well_statuses, well_types, year_slider){
               dff <- filter_Dataframe(df, well_statuses, well_types, year_slider)
               gas <- fetch_aggregate(dff,year_slider)[1]
               oil <- fetch_aggregate(dff,year_slider)[2]
               water <- fetch_aggregate(dff,year_slider)[3]
               aggregate = data.frame(table(dff$'Well_Type'))

               i <- 1
               label <- vector()
               for (j in aggregate$Var1) {
                 label[i] <- WELL_TYPES[j]
                 i <- i + 1
               }
               i <- 1
               colors <- list()
               for (j in aggregate$Var1) {

                 colors[[i]] <- WELL_COLORS[j]
                 i <- i + 1
               }

               data = list(
                 list(
                   type='pie',
                   labels=list('Gas', 'Oil', 'Water'),
                   values=list(sum(unlist(gas)), sum(unlist(oil)), sum(unlist(water))),
                   name='Production Breakdown',
                   text=list('Total Gas Produced (mcf)', 'Total Oil Produced (bbl)', 'Total Water Produced (bbl)'),  # noqa: E501
                   hoverinfo="text+value+percent",
                   textinfo="label+percent+name",
                   hole=0.5,
                   marker=list(
                     colors=list('#fac1b7', '#a9bb95', '#92d8d8')
                   ),
                   domain=list("x" = list(0, .45), 'y' = list(0.2, 0.8))
                 ),
                 list(
                   type='pie',
                   labels=label,
                   values=as.list(aggregate$Freq),
                   name='Well Type Breakdown',
                   hoverinfo="label+text+value+percent",
                   textinfo="label+percent+name",
                   hole=0.5,
                   marker=list(
                     colors=colors
                   ),
                   domain=list("x" = list(0.55, 1), 'y' = list(0.2, 0.8))
                 ))
               return (list(
                 'data' = data,
                 layout = list(
                   autosize=TRUE,
                   height=500,
                   font=list(color='#777777', size='7.5'),
                   titlefont=list(color='#777777', size='14'),
                   margin=list(
                     'l'=35,
                     'r'=35,
                     'b'=35,
                     't'=45
                   ),
                   hovermode="closest",
                   plot_bgcolor="#F9F9F9",
                   paper_bgcolor="#F9F9F9",
                   legend=list(font=list(size=10), orientation='h'),
                   title= sprintf('Production Summary: %s to %s', year_slider[1], year_slider[2]),
                   mapbox=list(
                     accesstoken=mapbox_access_token,
                     style="light",
                     center=list(
                       lon=-78.05,
                       lat=42.54
                     ),
                     zoom=7
                   ))))
             })

#AGGREGATE
app$callback(output=list(id='aggregate_graph', property='figure'),
             params=list(
               input(id='well_statuses', property='value'),
               input(id='well_types', property='value'),
               input(id='year_slider', property='value'),
               input(id='main_graph', property='hoverData')),

             function(well_statuses, well_types, year_slider,
                      main_graph_hover){
               chosen <- unlist(lapply(main_graph_hover[["points"]], "[[", "customdata"))
               vc <- chosen[1]
               welltype <- df[df$API_WellNo %in% vc,'Well_Type']
               dff <- filter_Dataframe(df, well_statuses, well_types, year_slider)
               selected <- dff[dff$Well_Type == welltype,]
               index <- fetch_abs(selected,year_slider)[1]
               gas <- fetch_abs(selected,year_slider)[2]
               gas1 <- gas[[1]][c(unlist(index))]
               oil <- fetch_abs(selected,year_slider)[3]
               oil1 <- oil[[1]][c(unlist(index))]
               water <- fetch_abs(selected,year_slider)[4]
               water1 <- water[[1]][c(unlist(index))]

               data = list(
                 list(
                   type='scatter',
                   mode='lines',
                   name='Gas Produced (mcf)',
                   x=c(unlist(index)),
                   y=c(unlist(gas1)),
                   line=list(
                     shape="spline",
                     smoothing="2",
                     color='#87CEEB'
                   )
                 ),
                 list(
                   type='scatter',
                   mode='lines',
                   name='Oil Produced (bbl)',
                   x=c(unlist(index)),
                   y=c(unlist(oil1)),
                   line=list(
                     shape="spline",
                     smoothing="2",
                     color='#849E68'
                   )
                 ),
                 list(
                   type='scatter',
                   mode='lines',
                   name='Water Produced (bbl)',
                   x=c(unlist(index)),
                   y=c(unlist(water1)),
                   line=list(
                     shape="spline",
                     smoothing="2",
                     color='#59C3C3'
                   )
                 )
               )
               return(list(
                 'data' = data,
                 'layout' = list(
                   autosize=TRUE,
                   height=500,
                   font=list(color='#777777'),
                   titlefont=list(color='#777777', size='14'),
                   margin=list(
                     'l'=35,
                     'r'=35,
                     'b'=35,
                     't'=45
                   ),
                   hovermode="closest",
                   plot_bgcolor="#F9F9F9",
                   paper_bgcolor="#F9F9F9",
                   title='Aggregate Production',
                   mapbox=list(
                     accesstoken=mapbox_access_token,
                     #style="dark",
                     center=list(
                       lon=-78.05,
                       lat=42.54
                     ),
                     zoom=7
                   ))
               ))
             })

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}

