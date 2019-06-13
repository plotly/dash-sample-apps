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
              className = "trigger-button",
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
              selected_className = "custom-tab--selected"
            ),
            dccTab(
              id = "Control-chart-tab",
              label = "Control Charts Dashboard",
              value = "tab2",
              className = "custom-tab",
              selected_className = "custom-tab--selected"
            )
          )
        )
      )
    )
  )
}

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
                children = "**What is this mock app about?**
                'dash-manufacture-spc-dashboard' is a dashboard for monitoring read-time process quality along manufacture production line.
                **What does this app show?**
                Click on buttons in 'Parameter' column to visualize details of measurement trendlines on the bottom panel.
                The Sparkline on top panel and Control chart on bottom panel show Shewhart process monitor using mock data.
                The trend is updated every other second to simulate real-time measurements. Data falling outside of six-sigma control limit are signals indicating 'Out of Control(OOC)', and will
                trigger alerts instantly for a detailed checkup."
              )
            )
          )
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
      # dccStore(
      #   id = "value-setter-store",
      #   data = init_value_setter_store()
      # ),
      generate_modal()
    )
  )
)

########################################################################################################################
# DEFINE APP CALLBACKS

# callbacks for modal popup when clicking "Learn More" 
# app$callback(
#   output = list(id = "markdown", property = "style"),
#   params = list(
#     input(id = "learn-more-button", property = "n_clicks"),
#     input(id = "markdown_close", property = "n_clicks")
#   ),
#   function(button_click, close_click) {
#     if (button_click > close_click) {
#       return(list("display" = "block"))
#     }
#     return(list("display" = "none"))
#   }
# )

# callbacks for switching tabs when button "Proceed To Measurement" clicked
# app$callback(
#   output = list(id = "app-tabs", property = "value"),
#   params = list(input(id = "tab-trigger-btn", property = "n_clicks")),
#   function(tab_switch) {
#     if (tab_switch == 0) {
#       return("tab1")
#     } else if (tab_switch) {
#       return("tab2")
#     }
#   }
# )

########################################################################################################################
# RUN SERVER

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
