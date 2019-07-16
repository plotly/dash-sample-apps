library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(chron)
library(dplyr)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("./app/apps/%s", appName))
}

app <- Dash$new()

setwd(sprintf("./app/apps/%s", appName))

df <- read.csv('./data/clinical_analytics.csv', header = TRUE, sep = ",")


#list of clinics
clinics = unique(df$Clinic.Name)
clinic_list <- list()
for (i in 1:length(clinics)){
  clinic_list[[i]] <- list(label = clinics[i], value = clinics[i])
}

#list of check in times
times = unique(df$Check.In.Time)
times_list <- list()
for (i in 1:length(clinics)){
  times_list[[i]] <- list(label = times[i], value = times[i])
}

#list of all departments 
deps = unique(df$Department)
all_departments = list()
for(i in 1:length(deps)){
  all_departments[[i]] = toString(deps[[i]])
}

df$Admit.Source[df$Admit.Source==""] <- "Not Identified"

df$Days.of.Wk = df$Check.In.Hour = df$Check.In.Time

#add Days of Wk column to df
#this is so slow sos
df$Days.of.Wk = weekdays(as.POSIXlt(df$Days.of.Wk))

#add Check In Hour column to df
printHour <- function(date_time_string) {
  s1 = strsplit(date_time_string, split=' ', fixed=TRUE)
  #s1 = "Y-M-D" "H:M:S" "AM/PM"
  
  am_pm = unlist(s1)[[3]]
  
  hour = strsplit(unlist(s1)[[2]], split=":", fixed = TRUE)
  
  hour = as.numeric(unlist(hour)[[1]])
  
  #want the form 01 PM or 01 AM
  if(hour<10){
    hour = paste(0,hour, sep="")
  }
  
  time = paste(hour, am_pm)
}

#also very slow :(
df$Check.In.Hour = lapply(df$Check.In.Hour, printHour)

df$Check.In.Time = as.Date(df$Check.In.Time)

admits = unique(df$Admit.Source)
admit_list <- list()
for (i in 1:length(admits)){
  admit_list[[i]] <- toString(admits[i])
}

#used for annotations for heat map
admit_list1 <- list()
for (i in 1:length(admits)){
  admit_list1[[i]] <- list(label = admits[i], value = admits[i])
}

y_axis = list('Sunday', 'Saturday', 'Friday', 'Thursday','Wednesday','Tuesday', 'Monday')

x_axis = list('12 AM', '01 AM', '02 AM', '03 AM', '04 AM', '05 AM', '06 AM', '07 AM', '08 AM', '09 AM', '10 AM', '11 AM', 
              '12 PM', '01 PM', '02 PM', '03 PM', '04 PM', '05 PM', '06 PM', '07 PM', '08 PM', '09 PM', '10 PM', '11 PM') 

### FUNCTIONS ###

#bold labels for the wait time graph
bold_labels <- function(department){
  dep <- unique(department)
  text <- list()
  for(i in 1:length(dep)){
    text[[i]] <- paste('<b>', toString(dep[[i]]),'</b>')
  }
  return(text)
}

#generate heat map
gen_heatmap1 = function(clinic_name,hm_click, start_date, end_date, admit_source){
  
  #filter df by clinic name, start date, end date, and admit list
  df_clinic = filter(df, df$Clinic.Name==clinic_name)
  date_sequence = seq(as.Date(start_date), as.Date(end_date), "days")
  df_date = filter(df_clinic,is.element(df_clinic$Check.In.Time, date_sequence))
  df_date = df_date[order(df_date$Check.In.Time),]
  df_date = filter(df_date,is.element(df_date$Admit.Source, admit_source))
  
  # z_axis holds the number of patient records for that given time and day
  z_axis = matrix(0, 7,24)
  annotations = list()
  k=1
  
  for(i in 1:length(y_axis)){
    df_day = filter(df_date, df_date$Days.of.Wk == y_axis[[i]])
    for(j in 1:length(x_axis)){
      df_hour = filter(df_day, df_day$Check.In.Hour==x_axis[[j]])
      sum = sum(df_hour$Number.of.Records)
      annotation = list(
        showarrow= FALSE,
        text = paste("<b>",toString(sum),"</b>"),
        xref = "x",
        x=x_axis[[j]],
        y = y_axis[[i]]
      )
      annotations[[k]]=annotation
      k = k+1
      z_axis[i,j] = sum
    }
  }
  
  # creates an orange box around the rectangle a user clicks on in the heatmap
  shapes = list()
  if(length(hm_click$points)==1){
    hour_of_day = hm_click$points[[1]]$x
    weekday = hm_click$points[[1]]$y
    x0 = (match(hour_of_day, x_axis)-1)/24
    x1 = x0+1 /24
    y0 = (match(weekday, y_axis)-1)/7
    y1 = y0+1 /7
    
    shapes = list(type = 'rect',
                  x0=x0,
                  x1=x1,
                  y0=y0,
                  y1=y1,
                  xref = 'paper',
                  yref = 'paper',
                  line = list(color = "#ff6347"))
  }
  
  #hover text for heatmap
  text1 = paste("<br>","%{y}" ,"%{x}", "<br><br>",
                "%{z}", " Patient Records")

  heat_map = plot_ly(x=x_axis, y = y_axis, z = z_axis,
                     type = 'heatmap', name = "",
                     hovertemplate = text1,
                     colorscale=list(list(0, "#caf3ff"), list(1, "#2c82ff")),
                     showscale = FALSE) %>%
    layout(
      yaxis = list(tickcolor = "rgb(250,250,250)"),
      xaxis = list(tickcolor = "rgb(250,250,250)", side = "top"),
      hovermode = "closest",
      showlegend = FALSE,
      annotations = annotations,
      modebar = list(orientation = "v"),
      shapes = shapes
    )
  return(heat_map)
}

#tick vals for y axis for wait time graph
make_tick_vals = function(dep_length){
  vals <- list()
  for(i in 1:dep_length){
    vals[[i]] = i-1
  }
  return(vals)
}

#generate wait time graph
generate_wait_time_graph = function(df1,clinic_name, admit_source, start_date, end_date){

  #filter by clinic, admit source, start date, end date
  df_filter = filter(df1, df1$Clinic.Name==clinic_name)
  df_filter = df_filter[is.element(df_filter$Admit.Source, admit_source),]
  date_sequence = seq(as.Date(start_date), as.Date(end_date), "days")
  df_filter = filter(df_filter,is.element(df_filter$Check.In.Time, date_sequence))

  #get average wait time and care score
  x1 <- df_filter %>%
    group_by(Department, Encounter.Number, Care.Score, Check.In.Time) %>%
    summarize(mean_wait = mean(Wait.Time.Min),
              mean_care_score = mean(Care.Score))
  
  #hover text
  text_wait = paste("Patient # : ", x1$Encounter.Number, 
                    "<br>Check in Time: ", x1$Check.In.Time,
                    "<br>Wait Time: ", x1$mean_wait, " Minutes <br>Care Score: ",
                    x1$mean_care_score)

  p1 = plot_ly(data = x1, x = ~mean_wait, 
               y = ~Department, 
               type = 'scatter',
               mode = "marker",
               hoverinfo = "text", text = text_wait,
               marker = list(size = 14,
                             line = list(width=1, color="#ffffff"),
                             color = "#2678b2")) %>%
    layout(
      title="",
      autosize= TRUE,
      yaxis = list(zeroline = FALSE,
                   autotick = FALSE,
                   showgrid = FALSE,
                   automargin = TRUE,
                   showline =  FALSE,
                   tickmode = 'list',
                   tickvals = make_tick_vals(length(unique(x1$Department))),
                   ticktext = bold_labels(unique(x1$Department))
                   
      ),
      xaxis = list(zeroline = FALSE,
                   showgrid = FALSE,
                   showline =  TRUE,
                   showticklabels = TRUE,
                   titlefont = list(color = "rgb(0,0,0)")
      ),
      
      paper_bgcolor="rgba(0,0,0,0)",
      plot_bgcolor="rgba(0,0,0,0)")
  
  return(p1)
}

#generate care score graph
generate_care_score_graph = function(df2,clinic_name, admit_source, start_date, end_date){
  #filter by clinic, admit source, start date, end date
  df_filter = filter(df2, df2$Clinic.Name==clinic_name)
  df_filter = df_filter[is.element(df_filter$Admit.Source, admit_source),]
  date_sequence = seq(as.Date(start_date), as.Date(end_date), "days")
  df_filter = filter(df_filter,is.element(df_filter$Check.In.Time, date_sequence))
  
  #get average care score and wait time
  x1 <- df_filter %>%
    group_by(Department, Encounter.Number, Check.In.Time) %>%
    summarize(mean_care_score = mean(Care.Score),
              mean_wait = mean(Wait.Time.Min))
  
  #hover text
  text_care = paste("Patient # : ", x1$Encounter.Number, 
                    "<br>Check in Time: ", x1$Check.In.Time,
                    "<br>Wait Time: ", x1$mean_wait, " Minutes <br>Care Score: ",
                    x1$mean_care_score)

  p2 = plot_ly(data = x1, x = ~mean_care_score, y = ~Department, type = 'scatter',
               hoverinfo = "text", text = text_care, 
               mode = "marker",
               marker = list(size = 14,
                             line = list(width=1, color="#ffffff"), 
                             color = "#2678b2")) %>%
    layout( title = '',
            autosize = TRUE,
      yaxis = list(zeroline = FALSE,
                   showgrid = FALSE,
                   showline =  FALSE,
                   showtitle = FALSE,
                   showticklabels = FALSE
      ),
      xaxis = list(zeroline = FALSE,
                   showgrid = FALSE,
                   showline =  TRUE,
                   showticklabels = TRUE,
                   titlefont = list(color = "rgb(0,0,0)"),
                   showtitle = FALSE,
                   range = list(1.5, 6.5)
      ),
      paper_bgcolor="rgba(0,0,0,0)",
      plot_bgcolor="rgba(0,0,0,0)", 
      showlegend = FALSE)
  
  return(p2)
}

#### PLOTLY OBJECTS ####
clinic_dropdown = dccDropdown(options = clinic_list,
                              value = 'Madison Center',
                              id = 'clinic-select')

check_in_time_picker = dccDatePickerRange(id = 'date-picker-select',
                                          start_date = as.Date("2014-1-1"),
                                          end_date = as.Date("2014-1-15"),
                                          min_date_allowed = as.Date("2014-1-1"),
                                          max_date_allowed = as.Date("2014-12-31"),
                                          initial_visible_month = as.Date("2014-1-1"))

admit_multi_select = dccDropdown(options = admit_list1,
                                 id="admit-select",
                                 value = admit_list,
                                 multi = TRUE)

hm_click = list()

starter_heat_map = gen_heatmap1("Madison Center", hm_click, "2014-1-1", "2014-1-15", admit_list)

starter_wait_graph = generate_wait_time_graph(df,"Madison Center", admit_list, "2014-1-1", "2014-1-15")

starter_care_graph = generate_care_score_graph(df,"Madison Center", admit_list, "2014-1-1", "2014-1-15")

wait_care_graph = subplot(starter_wait_graph,starter_care_graph)


### APP LAYOUT ###
app$layout(
  htmlDiv(
    id = 'app-container',
    list(
      
      #Banner
      htmlDiv(
        id = 'banner',
        list(
          htmlImg(src="https://dash.plot.ly/assets/images/logo.png")
        ),className = "banner"
      ),
      
      #Left column
      htmlDiv(
        id = 'left-column',
        list(
          htmlDiv(
            id = 'control-card',
            list(
              htmlH5("Clinical Analytics"),
              htmlH3("Welcome to the Clinical Analytics Dashboard"),
              htmlDiv(
                id ='intro',
                list(
                  "Explore clinic patient volume by time of day, waiting time, and care score. Click on the heatmap to visualize patient experience at different time points."
                )
              ),
              htmlP("Select Clinic"),
              clinic_dropdown,
              htmlBr(),
              htmlP("Select Check-in Time"),
              check_in_time_picker,
              htmlBr(),
              htmlBr(),
              htmlP("Select Admit Source"),
              admit_multi_select,
              htmlBr(),
              htmlButton(id='reset-btn', children = "Reset", n_clicks=0)
            )
          )
        ),className = 'four columns'
      ), 
      
      #Right Column
      htmlDiv(
        id = 'right-column',
        list(
          #Heatmap
          htmlDiv(
            id = 'patient_volume_card',
            list(
              htmlB("Patient Volume"),
              htmlHr(),
              htmlDiv(
                list(
                  dccGraph(
                    id = "patient_volume_hm",
                    figure = starter_heat_map,
                    clickData = list(points =list(), range =NULL)
                  )
                )
              )
            )
          ),
          #wait time + care score
          htmlBr(),
          htmlDiv(
            id = 'wait_time_card',
            list(
              htmlB("Patient Wait Time and Satisfactory Scores"),
              htmlHr(),
              #Header
              htmlDiv(
                id = 'header',
                list(
                  htmlDiv(
                    id = "header_department",
                    list(
                      htmlB("Department")
                    ),style = list(display = "table", height = "100%"), className = "two columns row-department"
                  ),
                  htmlDiv(
                    id = "header_wait_time_min",
                    list(
                      htmlB("Wait Time Minutes")
                    ), style = list(textAlign = "center", height = "100%"),className = "five columns"
                  ),
                  htmlDiv(
                    id = "header_care_score",
                    list(
                      htmlB("Care Score")
                    ),style = list(textAlign = "center", height = "100%"), className = "five columns"
                  )
                ), className = "row table-row"
              ),
              htmlDiv(
                list(
                  dccGraph(
                    id = "wait_time_graph",
                    figure = wait_care_graph
                  )
                )
              )
            )
          )
        ), className = 'eight columns'
      )
    )
  )
)

#### CALLBACKS ####

#callback for heatmap
app$callback(
  output = list(id = "patient_volume_hm", property = "figure"),
  params = list(input(id="patient_volume_hm", property = "clickData"),
                input(id = "date-picker-select", property = 'start_date'),
                input(id = 'date-picker-select', property = 'end_date'),
                input(id = "clinic-select", property = 'value'),
                input(id = "admit-select", property = "value"),
                input(id ='reset-btn', property = 'n_clicks')),
  function(click, start, end, clinic, admit_source, reset){
    if(reset ==1){
      reset=0
      return(gen_heatmap1(clinic, list(), "2014-1-1", "2014-1-15", admit_list))
    }
    else{
      map = gen_heatmap1(clinic, click, start, end, admit_source)
      return(map)
    }
  }
)

#callback for wait time graph
app$callback(
  output = list(id = "wait_time_graph", property = 'figure'),
  params = list(input(id = "patient_volume_hm", property = "clickData"),
                input(id = "date-picker-select", property = 'start_date'),
                input(id = 'date-picker-select', property = 'end_date'),
                input(id = "clinic-select", property = 'value'),
                input(id = "admit-select", property = "value"),
                input(id ='reset-btn', property = 'n_clicks')),
  function(click, start, end, clinic, admit_source, reset){
    
    if(reset==1){
      reset = 0
      care1 = generate_care_score_graph(df,clinic, admit_list, "2014-1-1", "2014-1-15")
      wait1 = generate_wait_time_graph(df,clinic, admit_list, "2014-1-1", "2014-1-15")
      wait_care1 = subplot(wait1, care1)
      return(wait_care1)
    }
    
    if(length(click$points) < 1){
      return(wait_care_graph)
    }
    else{
      df_fil = filter(df, df$Check.In.Hour==click$points[[1]]$x & df$Days.of.Wk==click$points[[1]]$y)
      care = generate_care_score_graph(df_fil,clinic, admit_source, start, end)
      wait = generate_wait_time_graph(df_fil,clinic, admit_source, start, end)
      wait_care = subplot(wait,care)
      return(wait_care)
    }
  }
)


app$run_server(showcase = TRUE)

