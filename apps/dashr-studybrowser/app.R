
appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

library(jsonlite)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(data.table)
library(stringr)
library(magrittr)
library(foreach)
library(rlist)
library(glue)
library(plotly)


#"Dict-List" of colors
group_colors <- list('negative'='light blue', 'test'='green', 'positive'='red')


returnOptions <- function(test_articles, study_data){ 
  studies <- (study_data[study_data$test_article==test_articles])$study_id %>% unique(.)
  foreach(i=studies) %do% list('label'=glue('{test_articles} (study: {i})'), 'value'=i)
}

asterisks <- function(n){
  if(0.1<=n){
    return(glue(''))
  } else if(0.05<=n && n<0.1){
    return(glue('. \n '))    
  } else if(0.01<=n && n<0.05){
    return(glue('*  \n '))
  } else if(0.001<=n && n<0.01){
    return(glue(' **  \n '))
  } else{
    return(glue('***  \n  '))
  }
  
}

default_study_data <- fread('assets/study.csv')
default_studyid <- "AS100"
default_study_input <- default_study_data[study_id==default_studyid]
default_test_articles <- (default_study_data[default_study_data$group_type=='test'] )$test_article %>% unique(.)
default_options <- foreach(i=default_test_articles, combine='c') %do% returnOptions(i, default_study_data)
default_options <- default_options[[1]]
app <- Dash$new()

default_plot <- plot_ly(
  default_study_input, 
  type='violin', 
  y= ~reading_value, 
  name=~group_name, 
  text=~subject_id, 
  hoveron='points',  
  meanline=list("visible"= T), 
  showlegend=F,
  points='all', 
  pointpos=0, 
  color=~group_colors)

app$layout(
  
  htmlDiv(children=list(
    htmlDiv(id='error-message'),
    #Title
    htmlDiv(#
      className='study-browser-banner row', 
      children=list(
        htmlH2(className="h2-title", children="ANIMAL STUDY BROWSER"),
        htmlImg(className="logo", src='assets/logo.png')
      ))#
    ,
    #App body
    htmlDiv(#
      className='row app-body', 
      children=list(
        htmlDiv(##
          className='four columns card', 
          children=list(
            htmlDiv(###
              className='bg-white', 
              children=list(
                htmlDiv(####
                  className='padding-top-bot',
                  children=list(
                    htmlEm('Test Article:'),
                    dccDropdown(id='study-dropdown', options=list(), value='')
                  )),####
                htmlDiv(####
                  className='padding-top-bot', 
                  children=list(
                    htmlH6("Choose the type of plot"),
                    dccRadioItems(#####
                                  id='chart-type', 
                                  options=list(######
                                               list('label'='Box Plot', 'value'='box'),  
                                               list('label'='Violin Plot', 'value'='violin')
                                  ),######
                                  value='violin',
                                  labelStyle=list(
                                    "display"= "inline-block",
                                    "padding"= "12px 12px 12px 0px"
                                  )
                    )#####
                  )),####
                htmlDiv(####
                  className='padding-top-bot',
                  children=list(
                    htmlH6("CSV File"),
                    dccUpload(#####
                              id='upload-data',
                              children= 
                                htmlDiv(#6
                                  className="upload",
                                  list(
                                    htmlP("Drag and Drop or "), 
                                    htmlA("Select Files")
                                  )),#6
                              filename = ''
                    ),#####
                    #htmlDiv(id='file-checker', children='None'),
                    htmlH5(" 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ", style=list(display='block', width='100%'))
                    
                  ))####
              ))###
          )),##
        #Graph
        htmlDiv(##
          className="eight columns card-left",
          children=list(
            htmlDiv(###
              className="bg-white",
              children=list(#
                htmlH5("Animal data plot"),
                dccGraph(id="plot", figure=default_plot)
              ))###
          ))##
      ))#
  )))

app$callback(
  output=list(id='error-message', property='children')
  ,
  params=list(input(id='upload-data', property='filename'))
  ,
  function(filename){
    if(is.na(filename)|| is.null(filename)){
      return("No file")
    } else if(filename==""){
      return("")
    } else{
      return(glue('Data {filename} has been successfully uploaded'))
    }
  }
)
#Decoding and parsing the file uploaded:
#helper function
parse_contents <- function(contents, filename){
  if(is.null(filename)){
    return(list(NULL, message='NO FILE'))
  }
  contents_parsed <- strsplit(contents, ',')
  decoded <- jsonlite::base64_dec(contents_parsed[[1]][2]) %>% rawToChar(.)
  if(str_detect(filename, '.csv')){
    DT <- fread(decoded)
    
  } else if(str_detect(filename, '.xls')){
    DT <- read_xls(decoded)
  } else{
    print("Invalid file")
    return(list(NULL, message='INVALID FILE'))
  }
  
  return(list(dt=DT, message=glue('Data {filename} has been successfully uploaded')))
}

app$callback(
  
  output=list(id='study-dropdown', property='options')
  ,
  params=list(input(id='upload-data', property='contents'), input(id='upload-data', property='filename'))
  ,
  function(contents, filename){
    
    
    if(is.na(filename)|| is.null(filename)||(filename=='')){
      
      study_data <- default_study_data
      test_articles <- (study_data[study_data$group_type=='test'] )$test_article %>% unique(.)
      #options = foreach(i=test_articles, combine='c') %do% returnOptions(i, study_data)
      options <- lapply(test_articles, returnOptions, study_data = study_data)
      return(options[[1]])
      
    } else {
      contents_parsed <- parse_contents(contents[[1]], filename)
      
      study_data <- contents_parsed$dt
      
      test_articles <- (study_data[study_data$group_type=='test'] )[['test_article']] %>% unique(.)
      #options = foreach(i=test_articles, combine='c') %do% returnOptions(i, study_data)
      options <- lapply(test_articles, returnOptions, study_data = study_data)
      return(options[[1]])
      
    }
  }
  
)


app$callback(
  
  output=list(id='plot', property='figure')
  ,
  params=list(input(id='chart-type',property='value' ), input(id='study-dropdown', property='value'), input(id='upload-data', property='contents'), input(id='upload-data', property='filename'))
  ,
  function(chart_type, study, contents, filename){
    
    if(is.na(filename)|| is.null(filename)||(filename=='')){
      
      study_data <- default_study_data
      
    }else{
      
      study_data <- parse_contents(contents[[1]], filename)$dt
      
    }
    
    if(is.null(study)||is.na(study)||study==''){
      study = study_data$study_id[1]
    }
    
    
    vehicle <- study_data[study_data$group_type == 'negative']
    study_data <- study_data[study_data$study_id == study]
    
    Y_DATA <- study_data %>% split(., study_data$group_id) %>% mapply(getElement, ., 'reading_value', SIMPLIFY = F)
    vehicle_readings <- vehicle[study_id==study]$reading_value
    
    T_TESTS <- ( Y_DATA ) %>% mapply(t.test,  ., lapply(1:length(.), function(n){return(vehicle_readings)}), SIMPLIFY=FALSE)
    
    T_stats <- T_TESTS %>% mapply(getElement, ., lapply(1:length(.),function(n){return('statistic')}), SIMPLIFY=F)
    
    P_values <- T_TESTS %>% mapply(getElement, ., lapply(1:length(.),function(n){return('p.value')}), SIMPLIFY=F)
    
    
    P_ast<-study_data$group_id %>% lapply(., function(n){return(P_values[[n]])}) %>% lapply(., asterisks) %>% unlist(.)
    
    group_title <- mapply(glue, P_ast, study_data$group_name)
    
    
    study_data$group_title = group_title
    study_data$Total_Scores <- study_data$reading_value
    p1 <- plot_ly(study_data, type='box', y= ~Total_Scores, name=~group_title, text=~subject_id, hoveron='points', boxmean=TRUE, showlegend=FALSE, boxpoints='all', pointpos=0, height=500)
    
    p2 <- plot_ly(study_data, type='violin', y= ~Total_Scores, name=~group_title, text=~subject_id, hoveron='points',  meanline=list("visible"= T), showlegend=F,points='all',
                  pointpos=0, height=500)
    chart_data <- list('box'= p1, 'violin'=p2)
    return(chart_data[[chart_type]])
  }
  
)

 app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
#app$run_server()

