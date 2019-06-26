library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(jsonlite)
library(data.table)

appName <- Sys.getenv("DASH_APP_NAME")
pathPrefix <- sprintf("/%s/", appName)

Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
           DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

setwd(sprintf("/app/apps/%s", appName))

source("utils/demo_utils.R")

demodf <- c("cifar","mnist","fashion")
simmodels <- c("softmax","cnn")

namess <- c("step", "train accuracy", "val accuracy", "train cross entropy", "val cross entropy")
data_dict <- list("softmax" = list("cifar" = fread('/data/cifar_softmax_run_log.csv',col.names=namess),
                                   "mnist" = fread('/data/mnist_softmax_run_log.csv',col.names=namess),
                                   "fashion" = fread('/data/fashion_softmax_run_log.csv',col.names=namess)),
                  
                  "cnn" = list("cifar" = fread('/data/cifar_cnn_run_log.csv',col.names=namess),
                               "mnist" = fread('/data/mnist_cnn_run_log.csv',col.names=namess),
                               "fashion" = fread('/data/fashion_cnn_run_log.csv',col.names=namess))
)

app <- Dash$new(name="DashR Live Model Training")

demo_mode <- TRUE

# Function to plot the 2 subgraphs at the bottom
div_graph <- function(name){
  return (htmlDiv(
    className="row",
    children = list(
      htmlDiv(
        className="two columns",
        style = list("padding-bottom" = "5%"),
        children = list(
          htmlDiv(list(
            htmlP(
              className="graph-checkbox-smoothing",
              children = list("Smoothing:"),
              style= list("font-weight"="bold", "margin-bottom"="0px")
            ),
            dccChecklist(
              options = list(
                list(label = "Training", value = "train"),
                list(label = "Validation", value = "val")
              ),
              value = list(),
              id = paste0("checklist-smoothing-options-", name),
              className = "checklist-smoothing"
            )
          ),
          style = list("margin-top" = "10px")
          ),
          
          htmlDiv(list(
            dccSlider(
              min = 0,
              max = 1,
              step = 0.05,
              marks = as.list(setNames(c(paste((0:5)/5)), ((0:5)/5))),
              value = 0.6,
              updatemode = "drag",
              id = paste0("slider-smoothing-", name)
            )
          ), 
          style = list("margin-bottom" = "40px"),
          className="slider-smoothing"
          ),
          
          htmlDiv(list(
            htmlP(
              "Plot Display mode:",
              style=list("font-weight"="bold",
                         "margin-bottom"="0px")
            ),
            dccRadioItems(
              options= list(
                list(label = "Overlapping", value = "overlap"),
                list(label = "Seperate (Vertical)", value = "seperate_vertical"),
                list(label = "Seperate (Horizontal)", value = "seperate_horizontal")
              ),
              value = "overlap",
              id = paste0("radio-display-mode-",name),
              className="plot-display-radio-items"
            ),
            
            htmlDiv(id= paste0("div-current-", name, "-value"))
          ))
        )
      ),
      htmlDiv(id = paste0("div-",name,"-graph"), className="ten columns")
    )
  ))
} #end of div_graph

app$layout(
  htmlDiv(
    style=list("height"="100%"),
    children=list(
      
      #Banner Display
      htmlDiv(list(
        htmlH2(
          'Live Model Training Viewer',
          id = "title",
          className = "eight columns",
          style = list("margin-left" = "3%")
        ),
        htmlButton(
          id="learn-more-button",
          className="two columns",
          children="Learn More",
          n_clicks=0
        ),
        htmlImg(
          src="/assets/dash-logo.png",
          className="two columns"
        )
      ), className="banner row"),
      
      
      #Learn More popup
      htmlDiv(
        id="markdown",
        className="model",
        style=list(display="none"),
        children=list(
          htmlDiv(
            className="close-container",
            children=htmlButton(
              children="Close",
              id="markdown_close",
              n_clicks=0,
              className='closeButton',
              style=list(height="100%")
            )
          ),
          htmlDiv(
            className="markdown-text",
            children=list(
              dccMarkdown(
                readLines("demo.md")
              )
            )
          )
        )
      ),
      
      #Body
      htmlDiv(
        className="container",
        style = list("padding"="35px 25px"),
        children = list(
          
          # Hidden Div that will store the result of simulating a model run
          dccStore(id="storage-simulated-run",storage_type='memory'),
          
          # Increment the simulation step count at a fixed time interval
          dccInterval(
            id="interval-simulated-step", 
            interval=125,
            n_intervals=0
          ),
          htmlDiv(
            className="row",
            style=list(
              "margin-bottom"="8px 0px"
            ),
            children = list(
              htmlDiv(
                className="twelve columns",
                children=list(
                  htmlDiv(
                    className="eight columns",
                    children = list(
                      htmlDiv(
                        dccDropdown(
                          id="dropdown-demo-dataset",
                          options=list(
                            list(label="CIFAR10", value="cifar"),
                            list(label="MNIST", value="mnist"),
                            list(label="Fashion MNIST", value="fashion")
                          ),
                          placeholder = "Select a demo dataset",
                          value = '',
                          searchable = FALSE
                        ),
                        className="six columns dropdown-box",
                      ),
                      
                      htmlDiv(
                        dccDropdown(
                          id = "dropdown-simulation-model",
                          options=list(
                            list(label = "1-Layer Neural Net", value = "softmax"),
                            list(label="Simple Conv Net", value = "cnn")
                          ),
                          placeholder = "Select Model to Simulate",
                          value = '',
                          searchable = FALSE,
                        ),
                        className="six columns dropdown-box"
                      ),
                      
                      htmlDiv(
                        dccDropdown(
                          id="dropdown-interval-control",
                          options=list(
                            list(label="No Updates", value= "no"),
                            list(label= "Slow Updates", value= "slow"),
                            list(label= "Regular Updates", value= "regular"),
                            list(label= "Fast Updates", value= "fast")
                          ),
                          value="regular",
                          className="twelve columns dropdown-box",
                          clearable=FALSE,
                          searchable=FALSE,
                        )
                      )
                    )
                  ), 
                  htmlDiv(
                    id="div-interval-control", 
                    className="four columns",
                    children = list(
                      htmlDiv(
                        id="div-total-step-count",
                        className="twelve columns",
                      ),
                      htmlDiv(
                        id="div-step-display",
                        className="twelve columns",
                      )
                    )
                  )
                )
              )
            )
          ),
          dccInterval(id="interval-log-update", n_intervals=0),
          
          # Hidden Div Storing JSON-serialized dataframe of run log
          dccStore(id="run-log-storage",storage_type='memory')
        )
      ),
      htmlDiv(
        className="container", 
        children=list(div_graph("accuracy"))
      ),
      
      htmlDiv(
        className="container",
        style=list("margin-bottom"="30px"),
        children = list(div_graph("cross-entropy"))
      )
    )
  )) #end layout

demo_callbacks(app, demo_mode)

# checklist_smoothing_options missing from function below. add when dccChecklist is fixed
update_graph <- function(graph_id,
                         graph_title,
                         y_train_index,
                         y_val_index,
                         run_log_json,
                         display_mode,
                         checklist_smoothing_options,
                         slider_smoothing,
                         yaxis_title){
  
  smooth <- function(scalars, weight=0.6){
    last <- scalars[[1]]
    smoothed <- list()
    for(point in scalars){
      smoothed_val <- last * weight + (1 - weight) * point
      smoothed <- append(smoothed,smoothed_val)
      last <- smoothed_val
    }
    return(as.double(smoothed))
  }
  
  if(!is.null(run_log_json)){
    
    run_log_df <- fromJSON(run_log_json)
    step <- run_log_df$step
    y_train <- run_log_df[[y_train_index]]
    y_val <- run_log_df[[y_val_index]]
    
    
    y_train <- smooth(y_train, weight=slider_smoothing)
    y_val <- smooth(y_val, weight=slider_smoothing)
    
    #Apply Smoothing if needed
    if (isTRUE(checklist_smoothing_options == "train")){
      y_train <- smooth(y_train, weight=slider_smoothing)
    }
    if (isTRUE(checklist_smoothing_options == "val")){
      y_val <- smooth(y_val, weight=slider_smoothing)
    }
    
    if(display_mode == "overlap"){
      fig <- plot_ly(run_log_df,type='scatter',
                     x = ~step,
                     y = y_train,
                     mode='lines',
                     name='Training',
                     line=list(color="rgb(54, 218, 170)"))
      
      fig <- add_trace(fig, run_log_df,
                       x = ~step,
                       y = y_val,
                       mode='lines',
                       name='Validation',
                       line=list(color="rgb(246, 236, 145)")) %>% layout(title=graph_title,
                                                                         margin=list(
                                                                           'l'=50,
                                                                           'r'=50,
                                                                           'b'=50,
                                                                           't'=50
                                                                         ),
                                                                         yaxis=list(title=yaxis_title))
    } else if (display_mode == "seperate_vertical"){
      p1 <- plot_ly(run_log_df,type='scatter',
                    x = ~step,
                    y = y_train,
                    mode='lines',
                    name='Training',
                    line=list(color="rgb(54, 218, 170)"))
      
      p2 <- plot_ly(run_log_df,type='scatter',
                    x = ~step,
                    y = y_val,
                    mode='lines',
                    name='Validation',
                    line=list(color="rgb(246, 236, 145)"))
      
      fig <- subplot(p1, p2, nrows=2)%>% layout(title=graph_title,
                                                margin=list(
                                                  'l'=50,
                                                  'r'=50,
                                                  'b'=50,
                                                  't'=50),
                                                yaxis=list(title=yaxis_title),
                                                yaxis2=list(title=yaxis_title)
      )
      
    }
    else if (display_mode == "seperate_horizontal"){
      
      p1 <- plot_ly(run_log_df,type='scatter',
                    x = ~step,
                    y = y_train,
                    mode='lines',
                    name='Training',
                    line=list(color="rgb(54, 218, 170)"))
      
      p2 <- plot_ly(run_log_df,type='scatter',
                    x = ~step,
                    y = y_val,
                    mode='lines',
                    name='Validation',
                    line=list(color="rgb(246, 236, 145)"))
      
      fig <- subplot(p1, p2, nrows=1)%>% layout(title=graph_title,
                                                margin=list(
                                                  'l'=50,
                                                  'r'=50,
                                                  'b'=50,
                                                  't'=50),
                                                yaxis=list(title=yaxis_title),
                                                yaxis2=list(title=yaxis_title)
                                                
      )
    }
  }
  return(dccGraph(id=graph_id,figure=fig))
}

#"Learn more" popup
app$callback(
  output = list(id = 'markdown', property = 'style'),
  list(input(id =  "learn-more-button", property = 'n_clicks'),
       input(id = "markdown_close", property = "n_clicks")),
  opcl_button <- function(button_click, close_click){
    if(button_click > close_click){
      return(list(display = 'block', backgroundColor = '#F9F9F9'))
    }
    else{
      return(list(display = 'none'))
    }
  }
)

#Update log
app$callback(
  output("interval-log-update", "interval"),
  list(input("dropdown-interval-control","value")),
  
  update_interval_log_update <- function(interval_rate){
    if(interval_rate =="fast"){
      return(500)
    }
    else if(interval_rate=="regular"){
      return(1000)
    }
    else if(interval_rate=="slow"){
      return(5*1000)
    }
    else if(interval_rate=="no"){
      return(24*60*60*1000)
    }
  }
)

#Steps Counter
app$callback(
  output("div-step-display","children"),
  list(input("run-log-storage","data")),
  
  update_div_step_display <- function(run_log_json){
    if(!is.null(run_log_json)){
      run_log_df <- fromJSON(run_log_json)
      return(htmlH6(paste0("Steps:",run_log_df$step[length(run_log_df$step)]),
                    style= list("margin-top" = "3px",
                                "float" = "right")
      ))
    }
  }
)

#Update the Prediction Accuracy graph
#dccChecklist is working again
app$callback(
  output(id="div-accuracy-graph", "children"),
  list(input("run-log-storage", "data"),
       input("radio-display-mode-accuracy", "value"),
       input("checklist-smoothing-options-accuracy", "value"),
       input("slider-smoothing-accuracy", "value")),
  
  update_accuracy_graph <- function(run_log_json, display_mode, checklist_smoothing_options, slider_smoothing){
    
    graph <- update_graph(
      "accuracy-graph",
      "Prediction Accuracy",
      "train accuracy",
      "val accuracy",
      run_log_json,
      display_mode,
      checklist_smoothing_options,
      slider_smoothing,
      "Accuracy"
    )
    
    if(isTRUE(display_mode == "overlap")) {
      graph$props$figure$x$layoutAttrs[[1]]$yaxis[["range"]]=list(0,1)
    }
    else if (isTRUE(display_mode == "seperate_horizontal") || isTRUE(display_mode == "seperate_vertical")){
      graph$props$figure$x$layoutAttrs[[1]]$yaxis[["range"]]=list(0,1)
      graph$props$figure$x$layoutAttrs[[1]]$yaxis2[["range"]]=list(0,1)
    }
    
    return(graph)
    
  }
)

#Callback for Cross entropy graph
app$callback(
  output(id="div-cross-entropy-graph", "children"),
  list(input("run-log-storage", "data"),
       input("radio-display-mode-cross-entropy", "value"),
       input("checklist-smoothing-options-cross-entropy", "value"),
       input("slider-smoothing-cross-entropy", "value")),
  
  update_cross_entropy_graph <- function(run_log_json, display_mode, checklist_smoothing_options, slider_smoothing){
    
    graph <- update_graph(
      "cross-entropy-graph",
      "Cross Entropy Loss",
      "train cross entropy",
      "val cross entropy",
      run_log_json,
      display_mode,
      checklist_smoothing_options,
      slider_smoothing,
      "Loss"
    )
    return(graph)
  }
)

#callback to display training and validation accuracy values
app$callback(
  output("div-current-accuracy-value", "children"),
  list(input("run-log-storage","data")),
  
  update_div_current_accuracy_value <- function(run_log_json){
    if(!is.null(run_log_json)){
      run_log_df <- fromJSON(run_log_json)
      return(list(
        htmlP("Current Accuracy:",style=list("font-weight"="bold",
                                             "margin-top"="15px",
                                             "margin-bottom"="0px")),
        htmlDiv(paste0("Training:",sprintf('%.4f',run_log_df[["train accuracy"]][length(run_log_df[["train accuracy"]])]))),
        htmlDiv(paste0("Validation:",sprintf('%.4f',run_log_df[["val accuracy"]][length(run_log_df[["val accuracy"]])])))
      )
      )
    }
  }
)

#callback to display training and validation cross entropy loss
app$callback(
  output("div-current-cross-entropy-value", "children"),
  list(input("run-log-storage","data")),
  
  update_div_current_cross_entropy_value <- function(run_log_json){
    if(!is.null(run_log_json)){
      run_log_df <- fromJSON(run_log_json)
      return(list(
        htmlP("Current Loss:",style=list("font-weight"="bold",
                                         "margin-top"="15px",
                                         "margin-bottom"="0px")),
        htmlDiv(paste0("Training:",sprintf('%.4f',run_log_df[["train cross entropy"]][length(run_log_df[["train cross entropy"]])]))),
        htmlDiv(paste0("Validation:",sprintf('%.4f',run_log_df[["val cross entropy"]][length(run_log_df[["val cross entropy"]])])))
      )
      )
    }
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}