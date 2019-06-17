library(rlist)
library(stringr)
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

app <- Dash$new(name = "DashR Manufacture SPC Dashboard", suppress_callback_exceptions = TRUE)

########################################################################################################################
# IMPORT DATA

data_dir <- "data/"
data_file <- "spc_data.csv"
data_path <- paste(data_dir, data_file, sep="")

df <- read.table(data_path, header = TRUE, sep = ",")

########################################################################################################################
# HANDLE & STORE DATA

params = names(df) # list of parameters (params[-1] removes "Batch")
max_length = nrow(df) # number of data entries

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
  for (i in 1:length(params)) {
    col <- params[i]
    data <- df[[col]]
    stats <- summary(data)
    ucl <- stats[["Mean"]]+3*sd(data)
    lcl <- stats[["Mean"]]-3*sd(data)
    usl <- stats[["Mean"]]+sd(data)
    lsl <- stats[["Mean"]]-sd(data)
    output[[col]] <- list(
      "count" = as.numeric(length(data)),
      "data" = data,
      "mean" = as.numeric(stats[["Mean"]]),
      "std" = as.numeric(sd(data)),
      "ucl" = round(ucl, 3),
      "lcl" = round(lcl, 3),
      "usl" = round(usl, 3),
      "lsl" = round(lsl, 3),
      "min" = as.numeric(stats[["Min."]]),
      "max" = as.numeric(stats[["Max."]]),
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
  "detail" = "#2d3038",  # background-card
  "primary" = "#007439",  # green
  "secondary" = "#FFD15F"  # accent
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

generate_modal <- function() {
  return(
    htmlDiv(
      id = "markdown",
      className = "modal",
      children = list(
        htmlDiv(
          id = "markdown-container",
          className = "markdown-container",
          children = list(
            htmlDiv(
              className = "close-container",
              children = htmlButton(
                children = "Close",
                id = "markdown_close",
                n_clicks = 0,
                className = "closeButton"
              )
            ),
            htmlDiv(
              className = "markdown-text",
              children = dccMarkdown(
                # CHANGE TO:
                # children = list(
                #   dccMarkdown("What is this mock app about?"),
                #   htmlPre("'dash-manufacture-spc-dashboard' is a dashboard for monitoring 
                #   real-time process quality along manufacture production line."),
                #   dccMarkdown("What does this app show?"),
                #   htmlPre("Click on buttons in 'Parameter' column to visualize
                #   details of measurement trendlines on the bottom panel.
                #   The Sparkline on top panel and Control chart on bottom panel show
                #   Shewhart process monitor using mock data.
                #   The trend is updated every other second to simulate real-time measurements.
                #   Data falling outside of six-sigma control limit are signals indicating 'Out of Control(OOC)',
                #   and will trigger alerts instantly for a detailed checkup.")
                # )
                children = "**What is this mock app about?**
                'dash-manufacture-spc-dashboard' is a dashboard for monitoring read-time process quality along manufacture production line.
                **What does this app show?**
                Click on buttons in 'Parameter' column to visualize details of measurement trendlines on the bottom panel.
                The Sparkline on top panel and Control chart on bottom panel show Shewhart process monitor using mock data.
                The trend is updated every other second to simulate real-time measurements. Data falling outside of six-sigma control 
                limit are signals indicating 'Out of Control(OOC)', and will trigger alerts instantly for a detailed checkup."
              )
            )
          )
        )
      )
    )
  )
}

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
          children = list(
            htmlButton(
              className = "trigger-button", # DELETE?
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

# generate layout for tab 1

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

build_tab_1 <- function() {
  return(
    list(
      htmlDiv(
        id = "set-specs-intro-container",
        children = htmlP("Use historical control limits to establish a benchmark, or set new values.")
      ),
      htmlDiv(
        id = "settings-menu",
        children = list(
          htmlDiv(
            id = "metric-select-menu",
            children = list(
              htmlLabel(
                id = "metric-select-title", 
                children = "Select Metrics"
              ),
              htmlBr(),
              dccDropdown(
                id = "metric-select-dropdown",
                options = lapply(params[-1], function(x) {return(list(label = x, value = x))}),
                value = params[[2]]
              )
            )
          ),
          htmlDiv(
            id = "value-setter-menu",
            children = list(
              htmlDiv(
                id = "value-setter-panel"
              ),
              htmlBr(),
              htmlButton(
                id = "value-setter-set-btn",
                children = "Update"
              ),
              htmlButton(
                id = "value-setter-view-btn",
                children = "View current setup",
                n_clicks = 0
              ),
              htmlDiv(
                id = "value-setter-view-output", 
                className = "output-datatable"
              )
            )
          )
        )
      )
    )
  )
}

# generate layout for tab 2

build_quick_stats_panel <- function() {
  return(
    htmlDiv(
      id = "quick-stats",
      className = "row",
      children = list(
        htmlDiv(
          id = "card-1",
          children = list(
            htmlH5("Operator ID"),
            daqLEDDisplay(
              value = "1704",
              color = "#91dfd2",
              backgroundColor = "#242633",
              size = "55"
            )
          )
        ),
        htmlDiv(
          id = "card-2",
          children = list(
            htmlH5("Time to completion"),
            daqGauge(
              id = "progress-gauge",
              value = 0,
              max = max_length * 2,
              min = 0,
              showCurrentValue = TRUE
            )
          )
        ),
        htmlDiv(
          id = "utility-card",
          children = list(
            daqStopButton(
              id = "stop-button", 
              size = 160,
              buttonText = "start"
            )
          )
        )
      )
    )
  )
}

generate_section_banner <- function(title) {
  return(
    htmlDiv(
      className = "section-banner", 
      children = title
    )
  )
}

generate_metric_row <- function(id, style, col1, col2, col3, col4, col5, col6) {
  if (length(style) == 0) {
    style = list(
      "height" = "10rem", 
      "width" = "100%"
    )
  }
  return(
    htmlDiv(
      id = id,
      className = "row metric-row",
      style = style,
      children = list(
        htmlDiv(
          id = col1[["id"]],
          className = "one column",
          style = list("margin-right" = "2.5rem"),
          children = col1[["children"]]
        ),
        htmlDiv(
          id = col2[["id"]],
          className = "one column",
          style = list("textAlign" = "center"),
          children = col2[["children"]]
        ),
        htmlDiv(
          id = col3[["id"]],
          className = "four columns",
          style = list("height" = "100%"),
          children = col3[["children"]]
        ),
        htmlDiv(
          id = col4[["id"]],
          className = "one column",
          style = list(),
          children = col4[["children"]]
        ),
        htmlDiv(
          id = col5[["id"]],
          className = "three columns",
          style = list("height" = "100%"),
          children = col5[["children"]]
        ),
        htmlDiv(
          id = col6[["id"]],
          className = "one column",
          style = list("display" = "flex", "justifyContent" = "center"),
          children = col6[["children"]]
        )
      )
    )
  )
}

generate_metric_list_header <- function() {
  return(
    generate_metric_row(
      "metric_header",
      list(
        "height" = "3rem", 
        "margin" = "1rem 0", 
        "textAlign" = "center"
      ),
      list(
        "id" = "m_header_1", 
        "children" = htmlDiv("Parameter")
      ),
      list(
        "id" = "m_header_2", 
        "children" = htmlDiv("Count")
      ),
      list(
        "id" = "m_header_3", 
        "children" = htmlDiv("Sparkline")
      ),
      list(
        "id" = "m_header_4", 
        "children" = htmlDiv("OOC%")
      ),
      list(
        "id" = "m_header_5", 
        "children" = htmlDiv("%OOC")
      ),
      list(
        "id" = "m_header_6", 
        "children" = "Pass/Fail"
      )
    )
  )
}



########################################################################################################################
# GENERATE APP LAYOUT



########################################################################################################################
# DEFINE APP CALLBACKS



########################################################################################################################
# RUN SERVER

# if (appName != "") {
#   app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
# } else {
#   app$run_server()
# }

