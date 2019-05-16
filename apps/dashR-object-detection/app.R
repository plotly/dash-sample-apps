require(tidyverse)
require(dashR)
require(dashCoreComponents)
require(dashHtmlComponents)
require(dashPlayer)
require(plotly)
require(rlist)
require(stringr)
require(data.table)
require(glue)
require(compiler)

DEBUG <<- T
FRAMERATE <<- 24.0
f_read <- cmpfun(fread)
data_table = cmpfun(data.table)

gum <- cmpfun(glue)

dict = cmpfun(list)

aplly <- cmpfun(apply)
laplly <- cmpfun(lapply)


load.data <<- cmpfun(function(path){
  #  # Load data about a specific footage (given by the path). It returns a dictionary of useful variables such as
  #  #  the dataframe containing all the detection and bounds localization, the number of classes inside that footage,
  #  #  the matrix of all the classes in string, the given class with padding, and the root of the number of classes,
  #  #  rounded.
  # 
  
  # Load DT containing all the processed object detections inside the video
  
  Info.DT <- f_read(path)
  
  #cls.vec <- the dict of detected object classes for the .csv file
  
  cls.vec <- Info.DT[, class_str]
  n.cls <- length(cls.vec)
  title.vec <- unique(cls.vec)
  
  # Gets the smallest value needed to add to the end of the classes dict to get a square matrix
  #   # n := argmin_{n : int} (n \geq n.c) \land (\sqrt{n} : int)
  #   # to find n we do
  
  r.rnd <- ceiling(sqrt(n.cls))
  
  pad.ent<- as.integer(r.rnd^2 - n.cls)
  pad.seq <- numeric(pad.ent)
  #Pad the class.vec
  cls.wpad <- c(cls.vec, pad.seq)
  cls.mtx <- aplly(matrix(cls.wpad, nrow=r.rnd, ncol=r.rnd, byrow=T),2,rev)
  
  
  
  dat_dict = dict(INFO_DT = Info.DT, N_CLS=n.cls, MTX_CLS = cls.mtx, WPAD_CLS = cls.wpad, RT_RND = r.rnd, TITL=title.vec)
  
  
  if(DEBUG)
  {
    print(gum('{path} loaded.'))
  }
  
  
  return(dat_dict)
  
})

#dict of optimized function through cmpfun
load_data <- cmpfun(load.data)



footage_labels <- c("james_bond", "zebra", "car_show_drone", "car_footage", "DroneCanalFestival","DroneCarFestival2","FarmDrone", "ManCCTV", "RestaurantHoldup")

dat.dict <<- dict( james_bond = load_data("data/james_bond_object_data.csv"), 
                         
                         zebra =  load_data("data/Zebra_object_data.csv"), 
                         
                         car_show_drone = load_data("data/CarShowDrone_object_data.csv"),  
                         
                         car_footage = load_data("data/CarFootage_object_data.csv"),
                         
                         DroneCanalFestival = load_data("data/DroneCanalFestivalDetectionData.csv"),
                         
                         DroneCarFestival2 =  load_data("data/DroneCarFestival2DetectionData.csv"),
                         
                         FarmDrone = load_data("data/FarmDroneDetectionData.csv"),
                         
                         ManCCTV = load_data("data/ManCCTVDetectionData.csv"),
                         
                         RestaurantHoldup =load_data("data/RestaurantHoldupDetectionData.csv"))




url.dict <<-  dict(regular=data_table(james_bond = 'https://www.youtube.com/watch?v=g9S5GndUhko', 
                                         
                                         zebra =  'https://www.youtube.com/watch?v=TVvtD3AVt10', 
                                         
                                         car_show_drone = 'https://www.youtube.com/watch?v=gPtn6hD7o8g',  
                                         
                                         car_footage = 'https://www.youtube.com/watch?v=qX3bDxHuq6I',
                                         
                                         DroneCanalFestival = 'https://youtu.be/0oucTt2OW7M',
                                         
                                         DroneCarFestival2 =  'https://youtu.be/vhJ7MHsJvwY',
                                         
                                         FarmDrone = 'https://youtu.be/aXfKuaP8v_A',
                                         
                                         ManCCTV =  'https://youtu.be/BYZORBIxgbc',
                                         
                                         RestaurantHoldup = 'https://youtu.be/WDin4qqgpac'), 
                      
                      bounding_box=data_table( james_bond = 'https://www.youtube.com/watch?v=g9S5GndUhko', 
                                               
                                               zebra =  'https://www.youtube.com/watch?v=G2pbZgyWQ5E', 
                                               
                                               car_show_drone = 'https://www.youtube.com/watch?v=9F5FdcVmLOY',  
                                               
                                               car_footage = 'https://www.youtube.com/watch?v=EhnNosq1Lrc',
                                               
                                               DroneCanalFestival = 'https://youtu.be/6ZZmsnwk2HQ',
                                               
                                               DroneCarFestival2 =  'https://youtu.be/2Gr4RQ-JHIs',
                                               
                                               FarmDrone = 'https://youtu.be/pvvW5yZlpyc',
                                               
                                               ManCCTV = 'https://youtu.be/1oMrHLrtOZw',
                                               
                                               RestaurantHoldup ='https://youtu.be/HOIKOwixYEY')
) 



app <- Dash$new()

markdown.text = "##### What am I looking at?


This app enhances visualization of objects detected using state-of-the-art Mobile Vision Neural Networks.
Most user generated videos are dynamic and fast-paced, which might be hard to interpret. A confidence
heatmap stays consistent through the video and intuitively displays the model predictions. The pie chart
lets you interpret how the object classes are divided, which is useful when analyzing videos with numerous
and differing objects.
                                
                              
##### More about this dash app


The purpose of this demo is to explore alternative visualization methods for Object Detection. Therefore,
the visualizations, predictions and videos are not generated in real time, but done beforehand. To read
more about it, please visit the [project repo](https://github.com/plotly/dash-object-detection)."

app$layout(
  
  
  
  htmlDiv(className='container', children=dict( #header-section
    #Top bar
    htmlDiv(
      id='top-bar',
      className='row',
      style=dict(backgroundColor = '#fa4f56', height= '5px')
    ),
    
    
    htmlDiv(
      id='left-side-column',
      className='eight columns',
      style=dict(display = 'flex', flexDirection = 'column', flex = 1, height = 'calc(100vh - 5px)', backgroundColor = '#F2F2F2', 'overflow-y' = 'scroll', marginLeft = '0px', justifyContent = 'flex-start', alignItems = 'center'),
      children=dict(
        
        
        htmlDiv(
          id='header-section',
          #title
          children=dict(
            
            
            ##
            htmlH4(
              'Object Detection Explorer'
            ), 
            
            
            #subheader
            htmlDiv(
              children = dict(dccMarkdown(glue('
        
        To get started, select a footage you want to view, and choose the display mode (with or without  bounding boxes).\ 
        
        Then, you can start playing the video, and the visualization will 
                                      be displayed depending on the current time.')))
              
            ),
            #'learn more' button - will contain the markdown_popup 
            htmlButton("Learn More", id="learn-more-button", n_clicks=0)
            
            
          )) ,
        
        
        #
        htmlDiv(
          id='markdown',
          className='model',
          style=dict(display= 'none'),
          children=dict(
            
            htmlDiv(
              className='close-container',
              children=htmlButton(
                'Close',
                id='markdown_close',
                n_clicks=0,
                className='closeButton',
                style=dict(border = 'none', height = '100%')
              )
            ),
            htmlDiv(
              className='markdown-text',
              children=dict(dccMarkdown(
                children= markdown.text
              )
              )
            ) #dict()
          )
          #dict()
          
        ),
        
        
        
        #Video-outer-container
        htmlDiv(
          className='video-outer-container',
          children= 
            htmlDiv(
              
              style=dict(width = '100%', paddingBottom = '56.25%', position = 'relative'),
              
              children=dashPlayer(
                
                id='video-display',
                
                style= dict(position = 'absolute', width= '100%', height = '100%', top = '0', left = '0', bottom = '0', right = '0'),
                
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
          children = dict(
            
            #LMNT #1
            htmlDiv(
              className='control-element',
              children=dict(
                htmlDiv(children=dict("Minimum Confidence Threshold:"), style=dict(width= '40%')),
                htmlDiv(
                  dccSlider(
                    id='slider-minimum-confidence-threshold',
                    min=20,
                    max=80,
                    step=NULL,
                    marks=dict('20'=dict(label='20%'), '30'=dict(label='30%'), '40'=dict(label='40%'), '50' = dict(label='50%'), '60' = dict(label = '60%'), '70' = dict(label='70%'), '80'= dict(label='80%')),
                    value=30,
                    updatemode='drag'
                  ), style=dict(width = '80%')) #60%
              )
            )
            , 
            
            
            #LMNT #2
            htmlDiv(
              className='control-element',
              children=dict(
                htmlP(children=dict("Footage Selection:"), style=dict(width = '40%')),
                dccDropdown(
                  id="dropdown-footage-selection",
                  options=dict(
                    dict(label = 'Drone recording of canal festival',
                         value = 'DroneCanalFestival'),
                    dict(label = 'Drone recording of car festival', value = 'car_show_drone'),
                    dict(label = 'Drone recording of car festival #2', value = 'DroneCarFestival2'),
                    dict(label = 'Drone recording of a farm', value = 'FarmDrone'),
                    dict(label = 'Lion fighting Zebras', value = 'zebra'),
                    dict(label = 'Man caught by a CCTV', value = 'ManCCTV'),
                    dict(label = 'Man driving expensive car', value = 'car_footage'),
                    dict(label = 'Restaurant Robbery', value = 'RestaurantHoldup')
                  ),
                  value='car_show_drone',
                  clearable=F,
                  style=dict(width = '60%')
                )
              )
            ),    
            
            
            #LMNT #3
            htmlDiv(
              className = 'control-element',
              
              children = dict(
                
                htmlDiv(children = dict("Video Display Mode:"), style = dict( width =  '40%')), dccDropdown(id = 'dropdown-video-display-mode',
                                                                                                            options =dict(
                                                                                                              dict(label='Regular Display', value='regular'), 
                                                                                                              dict(label='Bounding Boxes', value = 'bounding_box')
                                                                                                            ), value= 'bounding_box',
                                                                                                            searchable = F,
                                                                                                            clearable = F,
                                                                                                            style = dict(width='60%')
                )   
                
                
              )),
            
            #LMNT #4  
            htmlDiv(
              className='control-element',
              children=dict(
                htmlDiv(children="Graph View Mode:", style=dict(width = '60%')),
                dccDropdown(
                  id="dropdown-graph-view-mode",
                  options=dict(
                    dict(label = 'Visual Mode', value = 'visual'),
                    dict(label = 'Confidence Mode', value = 'confidence'),
                    dict(label = 'Detection Mode', value = 'detection')
                    
                  ),
                  value='visual',
                  searchable=FALSE,
                  clearable=FALSE,
                  style=dict(width = '60%')
                )
              )
            )
            
          ))   
        
      )),#LCOL end
    
    
    htmlDiv(
      id='right-side-column',
      className='four columns',
      style=dict( height = 'calc(100vh - 5px)',
                  overflowY = 'scroll',
                  marginLeft = '1%',
                  display = 'flex',
                  backgroundColor = '#F9F9F9',
                  flexDirection = 'column'),
      children=dict(
        htmlDiv(
          className='img-container',
          children=htmlImg(
            style=dict(height = '100%', margin = '1px'),
            src="https://s3-us-west-1.amazonaws.com/plotly-tutorials/logo/new-branding/dash-logo-by-plotly-stripe.png")
        ),
        
        htmlDiv(id="div-visual-mode" ),
        htmlDiv(id='div-confidence-mode'),
        htmlDiv(id="div-detection-mode")
        
      )
    )
    
  )) #CONTAINER end
  
)

#"Learn more" popup
app$callback(
  
  output = dict(
    
    id = 'markdown', property = 'style'
  ),
  
  params = dict(
    
    input(id =  "learn-more-button", property = 'n_clicks'),
    input(id = "markdown_close", property = "n_clicks")
    
  ),
  
  function(button_click, close_click){
    
    if(button_click > close_click){return(dict(display = 'block', backgroundColor = '#F9F9F9'))}
    else{return(dict(display = 'none'))}
    
  }
  
)


#Footage Selection
app$callback(
  
  output = dict(id ='video-display' , property='url'),
  params = dict(input(id='dropdown-footage-selection', property = 'value' ), input(id = 'dropdown-video-display-mode', property = 'value')),
  
  function(footage, display_mode){
    
    #Locate the selected footage and have player video updated
    
    url = url.dict[[display_mode]][[footage]]
    return(url)
    
  }
  
)

#Setting up the graphs for the visual mode
app$callback(
  output = dict(id = 'div-visual-mode', property = 'children'),
  
  params = dict(input(id='dropdown-graph-view-mode', property='value')), 
  
  function(dropdown_value){
    
    if(dropdown_value=="visual"){
      
      
      return(dict(
        
        dccInterval(id = 'interval-visual-mode', interval = 700, n_intervals = 0),
        
        htmlDiv(
          children=
            
            dict(
              
              htmlP(children = 'Object Count', className = 'plot_title'),
              
              dccGraph(id = 'pie-object-count', style = dict(height = '40vh', width = '100%') )
              
            ))
        
      ))
      
    } else {
      return(dict())
    }}
)

#Setting up the heatmap graph
app$callback(
  output = dict(id = 'div-confidence-mode', property = 'children'),
  
  params = dict(input(id='dropdown-graph-view-mode', property='value')), 
  
  function(dropdown_value){
    
    if(dropdown_value=="confidence"){
      
      
      return(dict(
        
        dccInterval(id = 'interval-confidence-mode', interval = 700, n_intervals = 0),
        
        htmlDiv(children = dict(htmlP(children = 'Confidence Level of Object Presence', className = 'plot_title'),
                                
                                dccGraph(id = 'heatmap-confidence', style = dict(height='40vh', width='100%') ) 
                                
                                
        )), 
        
        htmlDiv(id='heatmap-labels', children = dict() )
      ))
      
      
      
      
    } else {
      return(dict())
    }}
)


#Setting up the graphs for the detection mode
app$callback(
  
  output = dict(id = 'div-detection-mode', property = 'children'),
  
  params = dict(input(id = 'dropdown-graph-view-mode', property='value')),
  
  function(value){
    
    if(value=="detection"){ 
      
      return(
        
        
        dict(
          
          
          dccInterval(id = "interval-detection-mode", interval = 700, n_intervals = 0),
          
          htmlDiv(
            dict(
              
              htmlP(children =  "Detection Score of Most Probable Objects", className = 'plot-title'),
              
              dccGraph(id = 'bar-score-graph', style = dict(height = '55vh'))
              
              
            )
            
            
          )
          
          
        )
        
        
      )
      
    } else{
      
      return(dict()) # return empty dict/dict
    }}
)

#*****************************************************************************************************************


#Updating figures

#DETECTION MODE - bar graph

app$callback(
  
  output = dict(id = "bar-score-graph", property = "figure"),
  
  params = dict(input(id = "interval-detection-mode", property = "n_intervals"), 
                
                #State 
                state(id = 'video-display', property='currentTime'), 
                
                state(id = 'dropdown-footage-selection', property = 'value'),
                
                state(id =  'slider-minimum-confidence-threshold', property = 'value')
                
  ),
  
  
  
  function(n, current_time, footage, threshold){
    
    if(!is.null(current_time)){
      
      current_frame = round(current_time * FRAMERATE)
      
      if((n>0)&&(current_frame>0)){
        
        INFO_DT = dat.dict[[footage]]$INFO_DT
        
        #Select the subset of the dataset
        
        ##Extract the score
        frame_score <-  INFO_DT$score[(INFO_DT$frame == current_frame)]
        
        
        ##Extract the class str
        frame_string <- INFO_DT$class_str[(INFO_DT$frame == current_frame)]
        
        #Select only the frames above the threshold
        
        threshold_dec = threshold/100 #threshold in decimal
        
        #Filter them accordingly
        frame_score_filtered <- frame_score[frame_score > threshold_dec]
        frame_string_filtered <- as.character(frame_string[frame_score > threshold_dec])
        
        ##Save the length just in case 
        
        m <- 1:min(8, length(frame_score_filtered))
        
        frame_df <- data.table(cbind(frame_score_filtered, frame_string_filtered))
        colnames(frame_df) <- c("Score","Class")
        
        #Select up to top 8 frames w.r.t scores
        
        frame_ordered = frame_df[order(frame_score_filtered, decreasing=TRUE),]
        
        top_frames <- frame_ordered[m, ]
        
        
        #Add count to object names (e,g, person --> person 1, person --> person 2)
        
        Objects <- top_frames$Class
        
        set.Obj <- unique(Objects)
        
        Obj_count_lst <- as.list(numeric(length(set.Obj)))
        
        names(Obj_count_lst) <- set.Obj
        Obj_wc = c()
        
        
        
        
        for (Obj in Objects){
          
          Obj_count_lst[[Obj]]= Obj_count_lst[[Obj]] + 1
          Obj_wc <- append(Obj_wc, glue("{Obj}  {Obj_count_lst[[Obj]]}"))
          
        }
        
        colors = rep('rgb(250, 79, 86)', length(Obj_wc))
        
        #Add text information
        
        #y_text <- unlist(laplly(unlist(top_frames$Score),function(val){return(glue("{round(val*100)} % confidence"))} ))
        
        
        
        figure <- plot_ly(type='bar',hoverinfo='x+text', TITL= "Detection Scores",   x = Obj_wc, marker = dict(color=colors), y = top_frames$Score)
        
        figure_xaxis <- dict(automargin=TRUE, tickangle = -45)
        figure_yaxis <- dict(automargin=TRUE, range=c(0,1), title = dict(text = 'Score'))
        
        figure %>% layout(showlegend = FALSE, autosize=TRUE,  paper_bgcolor = 'rgb(249,249,249)', plot_bgcolor = 'rgb(249,249,249)') 
        
        return(figure)
        
      }
      
      
    }
    #p = plot_ly(type =  'bar', showlegend=FALSE, paper_bgcolor = 'rgb(249,249,249)', plot_bgcolor = 'rgb(249,249,249)',   yaxis = dict(title = 'score', automargin=TRUE), range = dict(0, 1))
    
    p_xaxis <- dict(automargin=TRUE)
    p_yaxis <- dict(title = 'Score', automargin = TRUE, range= 0:1)
    p = plot_ly(type =  'bar', x=c("Empty 1","Empty 2" ), y=c(0,0))
    p %>% layout(showlegend=FALSE, paper_bgcolor='rgb(249,249,249)', plot_bgcolor='rgb(249,249,249)', xaxis = p_xaxis, yaxis = p_yaxis)
    
    
    return(p) # Returns an empty bar graph
    
  }
  
  
  
)

#CONFIDENCE MODE - heatmap graph
app$callback(
  output = dict(id = "heatmap-confidence", property = "figure"),
  
  params = dict(input(id = "interval-confidence-mode", property = "n_intervals"),
                state(id = 'video-display', property = 'currentTime'),
                state(id = 'dropdown-footage-selection', property = 'value'),
                state(id = 'slider-minimum-confidence-threshold', property = 'value')),
  
  function(n, current_time, footage, threshold){
    
    heatmap_margin = list(l=10, r=10, b=20, t=20, pad=4)
    bgc = 'rgb(249,249,249)'
    
    if(!is.null(current_time)){
      
      current_frame = round(current_time * FRAMERATE)
      
      if((n>0)&&(current_frame>0)){
        
        #Load variables from the data list
        
        video_info_df = dat.dict[[footage]]$INFO_DT
        
        clses_padded = dat.dict[[footage]]$WPAD_CLS
        
        root_rnd = dat.dict[[footage]]$RT_RND
        
        clses_mtx = dat.dict[[footage]]$MTX_CLS
        
        titles = dat.dict[[footage]]$TITL
        
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
        class_wodup = frame_string_filtered[!duplicated(frame_index_filtered)]
        score_wodup = frame_score_filtered[!duplicated(frame_index_filtered)]
        index_wodup = which(unlist(lapply(1:length(titles), function(i){if(is.element(titles[i], class_wodup)){return(1)} else{return(0)}}))>0)
        
        #Building "the checkerboard"
        
        #Checkerboard -> m \times m tibble
        
        #Finding m
        
        l <- length(titles)
        m <- ceiling(sqrt(l))
        checker_vec <- numeric(m^2)
        
        checker_vec[index_wodup] = score_wodup
        
        checker_brd <- matrix(round(checker_vec, 2), nrow=m)
        
        titles_pad <- c(titles, rep("", times=(m^2-l)))
        
        hover_text = titles_pad
        
        hover_text = matrix(hover_text, ncol = m, byrow=F)
        
        ##Resorted to ggplot2
        
        checker_tbl <- data.table(checker_brd) %>% rownames_to_column('X') %>% gather(Y, value, -X)%>% mutate(
          X= factor(X, levels=1:m),
          Y = factor(gsub("V", "", Y), levels=1:m),
          Z = hover_text
        )
        
        
        ##plot data
        theme(text=element_text(size=5))
        
        plot_score <- ggplot(checker_tbl, aes(X, Y))  + geom_tile(aes(fill=value)) + geom_text(aes(label=round(value, 2), size=2)) + scale_fill_gradient(low="white", high="red")
        
        
        figure <- ggplotly(plot_score, label.size=0.01, tooltip='text') %>% layout(showlegend = FALSE, autosize = FALSE, paper_bgcolor = bgc, plot_bgcolor = bgc, margin = list(l=10, r=10, b=20, t=20, pad=2), font = list(size=7)) 
        
        
        return(figure)
        
      }
      
    }
    p = plot_ly(type='heatmap')
    
    p %>% layout(showlegend = F, paper_bgcolor=bgc, plot_bgcolor=bgc, autosize=F, margin = heatmap_margin)
    return(p)
    
  }
  
  
)


app$callback(
  
  output = dict(id='heatmap-labels', property='children'),
  params = dict(input(id = "interval-confidence-mode", property = "n_intervals"),
                state(id = 'video-display', property = 'currentTime'),
                state(id = 'dropdown-footage-selection', property = 'value')),
  
  function(n, ct, footage){
    titles = dat.dict[[footage]]$TITL
    l <- length(titles)
    m <- ceiling(sqrt(l))
    titles_pad <- c(titles, rep("", times=(m^2-l)))
    
    hover_text = titles_pad
    
    hover_text = matrix(hover_text, ncol = m, byrow=F)
    return(htmlPre(paste('
    ', laplly(1:m, function(i){return(paste(laplly(1:m, function(j){glue('{i} {j} - {hover_text[i, j]}')}), collapse ='
    '))}), collapse='
    ') 
    ))
  }
  
  
)

#VISUAL MODE - Object count pie graph

app$callback(
  output = dict(id = 'pie-object-count', property = 'figure'),
  
  params = dict(input(id = 'interval-visual-mode', property = 'n_intervals'), 
                state(id = 'video-display',  property = 'currentTime'), 
                state( id = 'dropdown-footage-selection', property = 'value'), 
                state( id = 'slider-minimum-confidence-threshold', property =  'value') ),
  
  function(n, current_time, footage, threshold){
    
    bgc =  'rgb(249,249,249)'
    
    pie_margin = dict(l=10, r=10, t=15, b=15)
    
    if(!is.null(current_time)){
      
      current_frame =  round(current_time * FRAMERATE)
      
      if((n>0)&&(current_frame>0)){
        
        INFO_DT = dat.dict[[footage]]$INFO_DT
        
        ##Extract the score
        frame_score <-  INFO_DT$score[(INFO_DT$frame == current_frame)]
        
        
        ##Extract the class str
        frame_string <- INFO_DT$class_str[(INFO_DT$frame == current_frame)]
        
        #Select only the frames above the threshold
        
        threshold_dec = threshold/100 #threshold in decimal
        
        #Filter them accordingly
        frame_score_filtered <- frame_score[frame_score > threshold_dec]
        frame_string_filtered <- as.character(frame_string[frame_score > threshold_dec])
        frame_index_filtered <- frame_string[frame_score > threshold_dec]
        
        
        #Get the count of each object class
        
        cls_tbl <- table(frame_index_filtered)
        
        cls_counts <- cls_tbl[cls_tbl > 0]
        
        clses =  names(cls_counts)
        counts =  unlist(cls_counts)
        
        text = unlist(laplly(counts, function(val){ glue("{val} detected")}))
        
        # Set colorscale to piechart
        
        colorscale = c('#fa4f56', '#fe6767', '#ff7c79', '#ff908b', '#ffa39d', '#ffb6b0', '#ffc8c3', '#ffdbd7',
                       '#ffedeb', '#ffffff')
        pie =plot_ly(type = "pie", labels = clses, values=counts, hoverinfo="text+percent", textinfo = "label+percent", marker = dict(colors=colorscale[1:length(clses)]))
        
        
        pie %>% layout( showlegend=TRUE, paper_bgcolor=bgc, plot_bgcolor = bgc, autosize=FALSE, margin= pie_margin)
        return(pie)
      }
      
    }
    
    p <- plot_ly(type="pie")
    
    p %>% layout(showlegend=TRUE, paper_bgcolor=bgc, plot_bgcolor = bgc, autosize=FALSE, margin= pie_margin)
    return(p)
    
  }
)




app$run_server(port="8893")


