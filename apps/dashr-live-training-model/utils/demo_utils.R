library(dashCoreComponents)
library(dashHtmlComponents)
library(data.table)
library(jsonlite)

demo_callbacks <- function(app, demo_mode){
  if(demo_mode){
    
    demodf <- c("cifar","mnist","fashion")
    simmodels <- c("softmax","cnn")
    
    setwd("C:/Users/MatthiasL/Desktop/DATA/Plotly/RApps/dash-live-model-training-r")
    
    namess <- c("step", "train accuracy", "val accuracy", "train cross entropy", "val cross entropy")
    data_dict <- list(
      "softmax" = list("cifar" = fread('./data/cifar_softmax_run_log.csv',col.names=namess),
                       "mnist" = fread('./data/mnist_softmax_run_log.csv',col.names=namess),
                       "fashion" = fread('./data/fashion_softmax_run_log.csv',col.names=namess)),
      
      "cnn" = list("cifar" = fread('./data/cifar_cnn_run_log.csv',col.names=namess),
                   "mnist" = fread('./data/mnist_cnn_run_log.csv',col.names=namess),
                   "fashion" = fread('./data/fashion_cnn_run_log.csv',col.names=namess))
    )
    
    app$callback(
      output(id="storage-simulated-run","data"),
      list(input("interval-simulated-step","n_intervals"),
           state("dropdown-demo-dataset","value"),
           state("dropdown-simulation-model","value")),
      
      simulate_run <- function(n_intervals, demo_dataset, simulation_model) {
        if((is.element(demo_dataset ,demodf)) &&
           (is.element(simulation_model, simmodels))
           && n_intervals > 0){
          steps <- n_intervals * 5
          run_logs <- data_dict[[simulation_model]][[demo_dataset]]
          run_below_steps <- subset(run_logs, step <= steps)
          json <- toJSON(run_below_steps)
          return(json)
        }
      }
    )
    
    app$callback(
      output("interval-simulated-step","n_intervals"),
      list(input("dropdown-demo-dataset","value"),
           input("dropdown-simulation-model","value")),
      reset_interval_simulated_step <- function(demo_ds,demo_sm){
        return(0)
      }
    )
    
    
    app$callback(
      output("run-log-storage","data"),
      list(input("interval-log-update","n_intervals"),
           state("storage-simulated-run","data")),
      
      get_run_log <- function(n_intervs, simulated_run){
        if(!is.null(simulate_run)){
          return(simulated_run)
        }
      }
    )
    
    #Total steps counter
    app$callback(
      output("div-total-step-count","children"),
      list(input("dropdown-demo-dataset","value")),
      
      total_step_count <- function(dataset_name){
        if(!is.null(dataset_name)){
          dataset <- data_dict$softmax[[dataset_name]]
          return(htmlH6(paste0("Total steps: ",dataset$step[length(dataset$step)]),
                        style= list("margin-top" = "3px",
                                    "float" = "right")
          ))
        }
      }
    )
  }
}
