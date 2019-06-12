library(rlist)
library(glue)
library(plotly)
library(dashR)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dashTable)
library(dashDaq)
library(data.table)

########################################################################################################################
# INITIATE APP

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  Sys.setenv(
    DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
    DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix
  )
  setwd(sprintf("/app/apps/%s", appName))
}

app <- Dash$new(name = "DashR Manufacture SPC Dashboard")

########################################################################################################################
# IMPORT DATA

data_dir <- "data/"
data_file <- "spc_data.csv"
data_path <- paste(data_dir, data_file, sep="")

df <- read.table(data_path, header = TRUE, sep = ",")

params = names(df)
max_length = nrow(df)

########################################################################################################################
# DEFINE GLOBAL VARIABLES

suffix_row <- "_row"
suffix_button_id <- "_button"
suffix_sparkline_graph <- "_sparkline_graph"
suffix_count <- "_count"
suffix_ooc_n <- "_OOC_number"
suffix_ooc_g <- "_OOC_graph"
suffix_indicator <- "_indicator"

theme <- list(
  "dark" = TRUE,
  "detail" = "#2d3038",  # Background-card
  "primary" = "#007439",  # Green
  "secondary" = "#FFD15F"  # Accent
)

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

########################################################################################################################
# DEFINE FUNCTIONS FOR APP LAYOUT AND CALLBACKS

build_banner <- function() {
  return(
    htmlDiv(
      id = "banner",
      className = "banner",
      children = list(
        htmlDiv(
          id = "banner-text",
          children = list(
            htmlH5("Manufacturing SPC Dashboard"),
            htmlH6("Process Control and Exception Reporting")
          )
        ),
        htmlDiv(
          id = "banner-logo",
          children = list (
            htmlButton (
              id = "learn-more-button",
              children = "LEARN MORE",
              n_clicks = 0
            ),
            htmlImg(
              id = "logo",
              src = "https://s3-us-west-1.amazonaws.com/plotly-tutorials/logo/new-branding/dash-logo-by-plotly-stripe-inverted.png"
            )
          )
        )
      )
    )
  )
}

build_tabs <- function() {
  return(
    htmlDiv(
      id = "tabs",
      className = "row container scalable",
      children = list(
        dccTabs(
          id = "app-tabs",
          value = "tab1",
          className = "custom-tabs",
          children = list(
            dccTab(
              id = "Specs-tab",
              label = "Specification Settings",
              value = "tab1",
              className = "custom-tab",
              selected_className = "custom-tab--selected",
              disabled_style = list (
                "backgroundColor" = "#1d202d",
                "color" = "white",
                "border" = "#23262E solid 3px",
                "display" = "flex",
                "flex-direction" = "column",
                "alignItems" = "center",
                "justifyContent" = "center",
                "cursor" = "default"
              ),
              disabled = FALSE
            ),
            dccTab(
              id = "Control-chart-tab",
              label = "Control Charts Dashboard",
              value = "tab2",
              className = "custom-tab",
              selected_className = "custom-tab--selected",
              disabled_style = list(
                "backgroundColor" = "#1d202d",
                "color" = "white",
                "border" = "#23262E solid 3px",
                "display" = "flex",
                "flex-direction" = "column",
                "alignItems" = "center",
                "justifyContent" = "center",
                "cursor" = "default"
              ),
              disabled = FALSE
            )
          )
        )
      )
    )
  )
}

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

init_value_setter_store <- function() {
  state_dict <- init_df()
  return(state_dict)
}

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

########################################################################################################################
# GENERATE APP LAYOUT

app$layout(
  htmlDiv(
    id = "big-app-container",
    children = list(
      build_banner(),
      htmlDiv(
        id = "app-container",
        children = list(
          build_tabs(),
          htmlDiv(
            id = "app-content", 
            className = "container scalable"
          ),
          htmlButton(
            children = "Proceed to Measurement", 
            id = "tab-trigger-btn", 
            n_clicks = 0
          )
        )
      ),
      dccStore(
        id = "value-setter-store",
        data = init_value_setter_store()
      )
      #generate_modal()
    )
  )
)

########################################################################################################################
# DEFINE APP CALLBACKS

# app$layout(
#   htmlDiv(
#     id = "value-setter-panel",
#     children = list(
#       build_value_setter_line(
#         "value-setter-panel-header",
#         "Specs",
#         "Historical Value",
#         "Set new value"
#       ),
#       build_value_setter_line(
#         "value-setter-panel-usl",
#         "Upper Specification limit",
#         state_dict[["Diameter"]][["usl"]],
#         ud_usl_input #unstring
#       ),
#       build_value_setter_line(
#         "value-setter-panel-lsl",
#         "Lower Specification limit",
#         state_dict[["Diameter"]][["lsl"]],
#         ud_lsl_input #unstring
#       ),
#       build_value_setter_line(
#         "value-setter-panel-ucl",
#         "Upper Control limit",
#         state_dict[["Diameter"]][["ucl"]],
#         ud_ucl_input #unstring
#       ),
#       build_value_setter_line(
#         "value-setter-panel-lcl",
#         "Lower Control limit",
#         state_dict[["Diameter"]][["lcl"]],
#         ud_lcl_input #unstring
#       )
#     )
#   )
# )

app$run_server()