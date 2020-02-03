appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashPlayer)
library(plotly)
library(rlist)
library(stringr)
library(data.table)
library(glue)
library(compiler)
library(magrittr)
library(R.cache)


DEBUG <- T
FRAMERATE <- 24.0
#list of optimized function through cmpfun

footage_labels <- c("james_bond", "zebra", "car_show_drone", "car_footage", "DroneCanalFestival","DroneCarFestival2","FarmDrone", "ManCCTV", "RestaurantHoldup")
dat.list <- loadCache(pathname="assets/datList.Rcache")
url.list <- loadCache(pathname="assets/urlList.Rcache")
app <- Dash$new()

markdown.text <- "##### What am I looking at?


This app enhances visualization of objects detected using state-of-the-art Mobile Vision Neural Networks.  
Most user generated videos are dynamic and fast-paced, which might be hard to interpret. A confidence  
heatmap stays consistent through the video and intuitively displays the model predictions. The pie chart  
lets you interpret how the object classes are divided, which is useful when analyzing videos with numerous  
and differing objects.
                                
                              
##### More about this dash app


The purpose of this demo is to explore alternative visualization methods for Object Detection. Therefore,  
the visualizations, prelistions and videos are not generated in real time, but done beforehand. To read  
more about it, please visit the [project repo](https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-object-detection)."

app$layout(
  htmlDiv(children=list( #header-section
    #Top bar
    htmlDiv(
      id='top-bar',
      className='row'
    ),
    htmlDiv(className='container', children=list(
      htmlDiv(
        id='left-side-column',
        className='eight columns',
        children=list(
          htmlImg(
            id="logo-mobile",
            src="assets/logo-old.png"
          ),
          htmlDiv(
            id='header-section',
            #title
            children=list(
              ##
              htmlH4(
                'Object Detection Explorer'
              ), 
              #subheader
              htmlDiv(
                children = list(dccMarkdown('
        To get started, select a footage you want to view, and choose the display mode (with or without  bounding boxes).  
        Then, you can start playing the video, and the visualization will be displayed depending on the current time.'))
                
              ),
              #'learn more' button - will contain the markdown_popup 
              htmlButton("Learn More", id="learn-more-button", n_clicks=0)
            )) ,
          #Video-outer-container
          htmlDiv(
            className='video-outer-container',
            children= 
              htmlDiv(
                className='video-container',
                children=dashPlayer(
                  id='video-display',
                  url='https://www.youtube.com/watch?v=gPtn6hD7o8g',
                  controls=TRUE,
                  playing=FALSE,
                  volume=1,
                  width='100%',
                  height='100%'
                )
              )
          ),
          htmlDiv(
            className='control-section',
            children = list(
              
              #LMNT #1
              htmlDiv(
                className='control-element',
                children=list(
                  htmlDiv(children=list("Minimum Confidence Threshold:"), style=list(width= '40%')),
                  htmlDiv(dccSlider(
                      id='slider-minimum-confidence-threshold',
                      min=20,
                      max=80,
                      step=NULL,
                      marks=list('20'=list(label='20%'), 
                                 '30'=list(label='30%'), 
                                 '40'=list(label='40%'), 
                                 '50'=list(label='50%'), 
                                 '60' = list(label = '60%'), 
                                 '70' = list(label='70%'), 
                                 '80'= list(label='80%')
                                 ),
                      value=30,
                      updatemode='drag'
                    ))
                ))
              , 
              #LMNT #2
              htmlDiv(
                className='control-element',
                children=list(
                  htmlDiv(children=list("Footage Selection:")),
                  dccDropdown(
                    id="dropdown-footage-selection",
                    options=list(
                      list(label = 'Drone recording of canal festival',
                           value = 'DroneCanalFestival'),
                      list(label = 'Drone recording of car festival', value = 'car_show_drone'),
                      list(label = 'Drone recording of car festival #2', value = 'DroneCarFestival2'),
                      list(label = 'Drone recording of a farm', value = 'FarmDrone'),
                      list(label = 'Lion fighting Zebras', value = 'zebra'),
                      list(label = 'Man caught by a CCTV', value = 'ManCCTV'),
                      list(label = 'Man driving expensive car', value = 'car_footage'),
                      list(label = 'Restaurant Robbery', value = 'RestaurantHoldup')
                    ),
                    value='car_show_drone',
                    clearable=F
                  )
                )),    
              
              
              #LMNT #3
              htmlDiv(
                className = 'control-element',
                
                children = list(
                  
                  htmlDiv(children = list("Video Display Mode:")), 
                  dccDropdown(id = 'dropdown-video-display-mode',
                              options =list(
                                list(label='Regular Display', value='regular'),
                                list(label='Bounding Boxes', value = 'bounding_box')
                              ), value= 'bounding_box',
                              searchable = F,
                              clearable = F
                  )   
                )),
              #LMNT #4  
              htmlDiv(
                className='control-element',
                children=list(
                  htmlDiv(children="Graph View Mode:"),
                  dccDropdown(
                    id="dropdown-graph-view-mode",
                    options=list(
                      list(label='Visual Mode', value='visual'),
                      list(label = 'Detection Mode', value = 'detection')
                    ),
                    value='visual',
                    searchable=FALSE,
                    clearable=FALSE
                  )
                ))
            ))   
        )),#LCOL end
      htmlDiv(
        id='right-side-column',
        className='four columns',
        children=list(
          htmlDiv(
            className='img-container',
            children=htmlImg(
              id="logo-web",
              src="assets/logo-old.png")
          ),
          htmlDiv(id="div-visual-mode" ),
          htmlDiv(id="div-detection-mode")
        ))
    )),#CONTAINER end
    htmlDiv(
      id='markdown',
      className='modal',
      style=list(display= 'none'),
      children=list(
        
        htmlDiv(
          className='close-container',
          children=htmlButton(
            'Close',
            id='markdown_close',
            n_clicks=0,
            className='closeButton',
          )
        ),
        htmlDiv(
          className='markdown-text',
          children=list(dccMarkdown(children= markdown.text)))
      ))
  )) 
)

#"Learn more" popup
app$callback(

  output = list(

    id = 'markdown', property = 'style'
  ),

  params = list(

    input(id =  "learn-more-button", property = 'n_clicks'),
    input(id = "markdown_close", property = "n_clicks")

  ),

  function(button_click, close_click){

    if(button_click > close_click){return(list(display = 'block', backgroundColor = '#F9F9F9'))}
    else{return(list(display = 'none'))}

  }

)


#Footage Selection
app$callback(
  
  output = list(id ='video-display' , property='url'),
  params = list(input(id='dropdown-footage-selection', property = 'value' ), input(id = 'dropdown-video-display-mode', property = 'value')),
  
  function(footage, display_mode){
    
    #Locate the selected footage and have player video updated
    
    url = url.list[[display_mode]][[footage]]
    return(url)
    
  }
  
)

#Setting up the graphs for the visual mode
app$callback(
  output = list(id = 'div-visual-mode', property = 'children'),
  
  params = list(input(id='dropdown-graph-view-mode', property='value')), 
  
  function(dropdown_value){
    
    if(dropdown_value=="visual"){
      
      
      return(list(
        
        dccInterval(id = 'interval-visual-mode', interval = 700, n_intervals = 0),
        
        htmlDiv(
          children=
            
            list(
              htmlP(children = 'Confidence Level of Object Presence', className = 'plot_title'),
              
              dccGraph(id = 'heatmap-confidence', style = list(height='40vh', width='100%')),
              
              htmlP(children = 'Object Count', className = 'plot_title'),
              
              dccGraph(id = 'pie-object-count', style = list(height = '40vh', width = '100%') )
              
            ))
        
      ))
      
    } else {
      return(list())
    }}
)

#Setting up the heatmap graph
# app$callback(
#   output = list(id = 'div-confidence-mode', property = 'children'),
#   
#   params = list(input(id='dropdown-graph-view-mode', property='value')), 
#   
#   function(dropdown_value){
#     
#     if(dropdown_value=="confidence"){
#       
#       
#       return(list(
#         
#         dccInterval(id = 'interval-confidence-mode', interval = 700, n_intervals = 0),
#         
#         htmlDiv(children = list( ) 
#         )), 
#         
#         htmlDiv(id='heatmap-labels', children = list() )
#       ))    } else {
#       return(list())
#     }}
# )


#Setting up the graphs for the detection mode
app$callback(
  
  output = list(id = 'div-detection-mode', property = 'children'),
  
  params = list(input(id = 'dropdown-graph-view-mode', property='value')),
  
  function(value){
    
    if(value=="detection"){ 
      
      return(
        
        
        list(
          
          
          dccInterval(id = "interval-detection-mode", interval = 700, n_intervals = 0),
          
          htmlDiv(
            list(
              
              htmlP(children =  "Detection Score of Most Probable Objects", className = 'plot-title'),
              
              dccGraph(id = 'bar-score-graph', style = list(height = '55vh'))
              
              
            )
            
            
          )
          
          
        )
        
        
      )
      
    } else{
      
      return(list()) # return empty list/list
    }}
)

#*****************************************************************************************************************


#Updating figures

#DETECTION MODE - bar graph

app$callback(
  
  output = list(id = "bar-score-graph", property = "figure"),
  
  params = list(input(id = "interval-detection-mode", property = "n_intervals"), 
                
                #State 
                state(id = 'video-display', property='currentTime'), 
                
                state(id = 'dropdown-footage-selection', property = 'value'),
                
                state(id =  'slider-minimum-confidence-threshold', property = 'value')
                
  ),
  
  
  
  function(n, current_time, footage, threshold){
    
    if(!is.null(current_time)){
      
      current_frame = round(current_time * FRAMERATE)
      
      if((n>0)&&(current_frame>0)){
        
        INFO_DT = dat.list[[footage]]$INFO_DT
        
        #Select the subset of the dataset
        
        ##Extract the score
        frame_score <-  INFO_DT$score[(INFO_DT$frame == current_frame)]
        
        
        ##Extract the class str
        frame_string <- INFO_DT$class_str[(INFO_DT$frame == current_frame)]
        
        #Select only the frames above the threshold
        
        threshold_dec <- threshold/100 #threshold in decimal
        
        #Filter them accordingly
        frame_score_filtered <- frame_score[frame_score > threshold_dec]
        frame_string_filtered <- as.character(frame_string[frame_score > threshold_dec])
        
        ##Save the length just in case 
        
        m <- 1:min(8, length(frame_score_filtered))
        
        frame_df <- data.table(cbind(frame_score_filtered, frame_string_filtered))
        colnames(frame_df) <- c("Score","Class")
        
        #Select up to top 8 frames w.r.t scores
        
        frame_ordered <- frame_df[order(frame_score_filtered, decreasing=TRUE),]
        top_frames <- frame_ordered[m, ]
        
        
        #Add count to object names (e,g, person --> person 1, person --> person 2)
        
        Objects <- top_frames$Class
        
        set.Obj <- unique(Objects)
        
        Obj_count_lst <- as.list(numeric(length(set.Obj)))
        
        names(Obj_count_lst) <- set.Obj
        Obj_wc <- c()
      
        for (Obj in Objects){
          
          Obj_count_lst[[Obj]] <- Obj_count_lst[[Obj]] + 1
          Obj_wc <- append(Obj_wc, glue("{Obj}  {Obj_count_lst[[Obj]]}"))
          
        }
        
        colors = rep('rgb(250, 79, 86)', length(Obj_wc))
        
        figure <- plot_ly(type='bar',hoverinfo='x+text', name= "Detection Scores",   x = Obj_wc, marker = list(color=colors), y = top_frames$Score)
        
        figure_xaxis <- list(automargin=TRUE, tickangle = -45)
        figure_yaxis <- list(automargin=TRUE, range=c(0,1), title = list(text = 'Score'))
        
        figure %<>% plotly::layout(showlegend = FALSE, autosize=TRUE) 
        
        return(figure)
        
      }
      
      
    }
    p_xaxis <- list(automargin=TRUE)
    p_yaxis <- list(title = 'Score', automargin = TRUE, range= 0:1)
    p <- plot_ly(type =  'bar', x=c("Empty 1","Empty 2" ), y=c(0,0)) %>%
      plotly::layout(showlegend=FALSE, xaxis = p_xaxis, yaxis = p_yaxis)
    
    
    return(p) # Returns an empty bar graph
    
  }
  
  
  
)

#CONFIDENCE MODE - heatmap graph
app$callback(
  output = list(id = "heatmap-confidence", property = "figure"),
  
  params = list(input(id = "interval-visual-mode", property = "n_intervals"),
                state(id = 'video-display', property = 'currentTime'),
                state(id = 'dropdown-footage-selection', property = 'value'),
                state(id = 'slider-minimum-confidence-threshold', property = 'value')),
  
  function(n, current_time, footage, threshold){
    
    heatmap_margin <- list(l=10, r=10, b=20, t=20, pad=4)
    bgc <- 'rgb(249,249,249)'
    
    if(!is.null(current_time)){
      
      current_frame <- round(current_time * FRAMERATE)
      
      if((n>0)&&(current_frame>0)){
        
        #Load variables from the data list
        
        video_info_df <- dat.list[[footage]]$INFO_DT
        
        clses_padded <- dat.list[[footage]]$WPAD_CLS
        
        root_rnd <- dat.list[[footage]]$RT_RND
        
        clses_mtx <- dat.list[[footage]]$MTX_CLS
        
        titles <- dat.list[[footage]]$TITL
        
        ##Extract the score
        frame_score <-  video_info_df$score[(video_info_df$frame == current_frame)]
        
        
        ##Extract the class str
        frame_string <- video_info_df$class_str[(video_info_df$frame == current_frame)]
        
        #Select only the frames above the threshold
        
        threshold_dec = threshold/100 #threshold in decimal
        
        #Filter them accordingly
        frame_score_filtered <- frame_score[frame_score > threshold_dec]
        frame_string_filtered <- as.character(frame_string[frame_score > threshold_dec])
        frame_index_filtered <- frame_string[frame_score > threshold_dec]
        
        
        # Remove duplicate, keep the top result
        class_wodup <- frame_string_filtered[!duplicated(frame_index_filtered)]
        score_wodup <- frame_score_filtered[!duplicated(frame_index_filtered)]
        index_wodup <- which(unlist(lapply(1:length(titles), function(i){if(is.element(titles[i], class_wodup)){return(1)} else{return(0)}}))>0)
        
        #Building "the checkerboard"
        
        #Checkerboard -> m \times m tibble
        
        #Finding m
        
        l <- length(titles)
        m <- ceiling(sqrt(l))
        checker_vec <- numeric(m^2)
        
        checker_vec[index_wodup] = score_wodup
        
        checker_brd <- matrix(round(checker_vec, 2), nrow=m)
        
        titles_pad <- c(titles, rep("", times=(m^2-l)))
        
        hover_text <- titles_pad
        
        hover_text <- matrix(hover_text, ncol = m, byrow=F)
        
        ##Resorted to ggplot2
        
        figure <- plot_ly(type='heatmap',  x=1:m, y=1:m, z=checker_brd, text=hover_text, colors=colorRamp(c('white', 'red')))
        
        return(figure)
        
      }
      
    }
    p = plot_ly(type='heatmap') %>% 
      layout(showlegend = F, autosize=F, margin = heatmap_margin)
    return(p)
    
  }
  
  
)



#VISUAL MODE - Object count pie graph
app$callback(
  output = list(id = 'pie-object-count', property = 'figure'),
  
  params = list(input(id = 'interval-visual-mode', property = 'n_intervals'), 
                state(id = 'video-display',  property = 'currentTime'), 
                state( id = 'dropdown-footage-selection', property = 'value'), 
                state( id = 'slider-minimum-confidence-threshold', property =  'value') ),
  
  function(n, current_time, footage, threshold){
    
    
    
    pie_margin <- list(l=10, r=10, t=15, b=15)
    
    if(!is.null(current_time)){
      
      current_frame <- round(current_time * FRAMERATE)
      
      if((n>0)&&(current_frame>0)){
        
        INFO_DT <- dat.list[[footage]]$INFO_DT
        
        ##Extract the score
        frame_score <-  INFO_DT$score[(INFO_DT$frame == current_frame)]
        
        
        ##Extract the class str
        frame_string <- INFO_DT$class_str[(INFO_DT$frame == current_frame)]
        
        #Select only the frames above the threshold
        
        threshold_dec <- threshold/100 #threshold in decimal
        
        #Filter them accordingly
        frame_score_filtered <- frame_score[frame_score > threshold_dec]
        frame_string_filtered <- as.character(frame_string[frame_score > threshold_dec])
        frame_index_filtered <- frame_string[frame_score > threshold_dec]
        
        
        #Get the count of each object class
        
        cls_tbl <- table(frame_index_filtered)
        
        cls_counts <- cls_tbl[cls_tbl > 0]
        
        clses <- names(cls_counts)
        counts <-  unlist(cls_counts)
        
        text <- unlist(lapply(counts, function(val){ glue("{val} detected")}))
        
        # Set colorscale to piechart
        
        colorscale <- c('#fa4f56', 
                        '#fe6767', 
                        '#ff7c79', 
                        '#ff908b', 
                        '#ffa39d', 
                        '#ffb6b0', 
                        '#ffc8c3', 
                        '#ffdbd7',
                        '#ffedeb', 
                        '#ffffff')
        pie <- plot_ly(type = "pie", labels = clses, values=counts, hoverinfo="text+percent", textinfo = "label+percent", marker = list(colors=colorscale[1:length(clses)])) 
        pie %<>% plotly::layout( showlegend=TRUE, autosize=FALSE, margin= pie_margin)
        
        return(pie)
      }
      
    }
    
    p <- plot_ly(type="pie") %>%plotly::layout(showlegend=TRUE, autosize=FALSE, margin= pie_margin)
  
    return(p)
    
  }
)

app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))


