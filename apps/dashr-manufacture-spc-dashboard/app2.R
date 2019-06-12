library(rlist)
library(glue)
library(plotly)
library(dashR)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dashTable)
library(dashDaq)
library(data.table)

app <- Dash$new("DashR Test3")

data_dir <- "data/"
data_file <- "spc_data.csv"
data_path <- paste(data_dir, data_file, sep="")

df <- read.table(data_path, header = TRUE, sep = ",")

params = names(df)
max_length = nrow(df)

ud_usl_input <- daqNumericInput(
  id = "ud_usl_input",
  className = "setting-input",
  size = 200,
  max = 9999999
)
ud_lsl_input <- daqNumericInput(
  id = "ud_lsl_input",
  className = "setting-input",
  size = 200,
  max = 9999999
)
ud_ucl_input <- daqNumericInput(
  id = "ud_ucl_input",
  className = "setting-input",
  size = 200,
  max = 9999999
)
ud_lcl_input <- daqNumericInput(
  id = "ud_lcl_input",
  className = "setting-input",
  size = 200,
  max = 9999999
)

populate_ooc <- function(data, ucl, lcl) {
  ooc_count <- 0
  output <- list()
  for (i in 1:max_length) {
    if (data[[i]] >= ucl | data[[i]] <= lcl) {
      ooc_count <- ooc_count + 1
    }
    output[[length(output)+1]] <- ooc_count / (i+1)
  }
  return(output)
}

init_df <- function() {
  output <- list()
  for (i in 2:length(params)) {
    col <- params[i]
    data <- df[[col]]
    stats <- summary(data)
    ucl <- stats[["Mean"]]+3*sd(data)
    lcl <- stats[["Mean"]]-3*sd(data)
    usl <- stats[["Mean"]]+sd(data)
    lsl <- stats[["Mean"]]-sd(data)
    output[[col]] <- list(
      "count" = length(data),
      "data" = data,
      "mean" = stats[["Mean"]],
      "std" = sd(data),
      "ucl" = round(ucl, 3),
      "lcl" = round(lcl, 3),
      "usl" = round(usl, 3),
      "lsl" = round(lsl, 3),
      "min" = stats[["Min."]],
      "max" = stats[["Max."]],
      "ooc" = populate_ooc(data, ucl, lcl)
    )
  }
  return(output)
}

state_dict <- init_df()

build_value_setter_line <- function(line_num, label, value, col3) {
  return(
    htmlDiv(
      id = line_num,
      className = "row",
      children = list(
        htmlLabel(
          className = "four columns",
          children = label
        ),
        htmlLabel(
          className = "four columns",
          children = value
        ),
        htmlDiv(
          className = "four columns",
          children = col3 
        )
      )
    )
  )
}

app$layout(
  htmlDiv(
    id = "value-setter-panel",
    children = list(
      build_value_setter_line(
        "value-setter-panel-header",
        "Specs",
        "Historical Value",
        "Set new value"
      ),
      build_value_setter_line(
        "value-setter-panel-usl",
        "Upper Specification limit",
        state_dict[["Diameter"]][["usl"]],
        ud_usl_input #unstring
      ),
      build_value_setter_line(
        "value-setter-panel-lsl",
        "Lower Specification limit",
        state_dict[["Diameter"]][["lsl"]],
        ud_lsl_input #unstring
      ),
      build_value_setter_line(
        "value-setter-panel-ucl",
        "Upper Control limit",
        state_dict[["Diameter"]][["ucl"]],
        ud_ucl_input #unstring
      ),
      build_value_setter_line(
        "value-setter-panel-lcl",
        "Lower Control limit",
        state_dict[["Diameter"]][["lcl"]],
        ud_lcl_input #unstring
      )
    )
  )
)

app$run_server()