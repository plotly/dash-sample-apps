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

# TAB 2: useful functions

generate_section_banner <- function(title) {
  return(
    htmlDiv(
      className = "section-banner", 
      children = title
    )
  )
}

# TAB 2: build_quick_stats_panel

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

# TAB 2: build_top_panel (not done)

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
          style = list("textAlign" = "center"),
          className = "one column",
          children = col2[["children"]]
        ),
        htmlDiv(
          id = col3[["id"]],
          style = list("height" = "100%"),
          className = "four columns",
          children = col3[["children"]]
        ),
        htmlDiv(
          id = col4[["id"]],
          style = list(),
          className = "one column",
          children = col4[["children"]]
        ),
        htmlDiv(
          id = col5[["id"]],
          style = list("height" = "100%"),
          className = "three columns",
          children = col5[["children"]]
        ),
        htmlDiv(
          id = col6[["id"]],
          style = list("display" = "flex", "justifyContent" = "center"),
          className = "one column",
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
generate_metric_row_helper <- function(index) {
  item <- gsub("[.]", "-", params[-1][[index]])
  div_id <- glue("{item}{suffix_row}")
  button_id <- glue("{item}{suffix_button_id}")
  count_id <- glue("{item}{suffix_count}")
  sparkline_graph_id <- glue("{item}{suffix_sparkline_graph}")
  ooc_percentage_id <- glue("{item}{suffix_ooc_n}")
  ooc_graph_id <- glue("{item}{suffix_ooc_g}")
  indicator_id <- glue("{item}{suffix_indicator}")
  return(
    generate_metric_row(
      div_id,
      list(),
      list(
        "id" = item,
        "className" = "metric-row-button-text",
        "children" = htmlButton(
          id = button_id,
          className = "metric-row-button",
          children = item,
          title = "Click to visualize live SPC chart",
          n_clicks = 0
        )
      ),
      list(
        "id" = count_id,
        "children" = "0"
      ),
      list(
        "id" = glue("{item}_sparkline"),
        "children" = dccGraph(
          id = sparkline_graph_id,
          style = list(
            "width" = "100%", 
            "height" = "95%"
          ),
          config = list(
            "staticPlot" = FALSE,
            "editable" = FALSE,
            "displayModeBar" = FALSE
          ),
          figure = list(
            "data" = list(
              list(
                "type" = "scatter",
                "x" = list(),
                "y" = list(),
                "mode" = "lines+markers",
                "name" = item,
                "line" = list("color" = "#f4d44d")
              )
            ),
            "layout" = list(
              "uirevision" = TRUE,
              "margin" = list(l=0, r=0, t=4, b=4, pad=0),
              "paper_bgcolor" = "#1d202d",
              "plot_bgcolor" = "#1d202d"
            )
          )
        )
      ),
      list(
        "id" = ooc_percentage_id,
        "children" = "0.00%"
      ),
      list(
        "id" = glue("{ooc_graph_id}_container"),
        "children" = daqGraduatedBar(
          id = ooc_graph_id,
          color = list(
            "ranges" = list(
              "#91dfd2" = list(0, 3),
              "#f4d44d " = list(3, 7),
              "#f45060" = list(7, 15)
            )
          ),
          showCurrentValue = FALSE,
          size = 146.16,
          max = 15,
          value = 0
        )
      ),
      list(
        "id" = glue("{item}_pf"),
        "children" = daqIndicator(
          id = indicator_id, 
          value = TRUE, 
          color = "#91dfd2"
        )
      )
    )
  )
}
generate_piechart <- function() {
  return(
    dccGraph(
      id = "piechart",
      figure = list(
        "data" = list(
          list(
            "labels" = gsub("[.]", "-", params[-1]),
            "values" = list(1, 1, 1, 1, 1, 1, 1),
            #"values" = list(),
            "type" = "pie",
            "marker" = list(
              "colors" = list("#91dfd2","#91dfd2","#91dfd2","#91dfd2","#91dfd2","#91dfd2","#91dfd2"),
              "line" = list("color" = "white", "width" = 2)
            ),
            "hoverinfo" = "label",
            "textinfo" = "label"
          )
        ),
        "layout" = list(
          "uirevision" = TRUE,
          "showlegend" = FALSE,
          "paper_bgcolor" = "#1d202d",
          "plot_bgcolor" = "#1d202d",
          "font" = list("color" = "white"),
          "autosize" = TRUE
        )
      )
    )
  )
}
build_top_panel <- function() {
  return(
    htmlDiv(
      id = "top-section-container",
      className = "row",
      children = list(
        htmlDiv(
          id = "metric-summary-session",
          className = "eight columns",
          children = list(
            generate_section_banner("Process Control Metrics Summary"),
            htmlDiv(
              id = "metric-div",
              children = list(
                generate_metric_list_header(),
                htmlDiv(
                  id = "metric-rows",
                  children = list(
                    generate_metric_row_helper(1),
                    generate_metric_row_helper(2),
                    generate_metric_row_helper(3),
                    generate_metric_row_helper(4),
                    generate_metric_row_helper(5),
                    generate_metric_row_helper(6),
                    generate_metric_row_helper(7)
                  )
                )
              )
            )
          )
        ),
        htmlDiv(
          id = "ooc-piechart-outer",
          className = "four columns",
          children = list(
            generate_section_banner("% OOC per Parameter"),
            generate_piechart()
          )
        )
      )
    )
  )
}

# TAB 2: build_chart_panel (not done?)

build_chart_panel <- function() {
  return(
    htmlDiv(
      id = "control-chart-container",
      className = "twelve columns",
      children = list(
        generate_section_banner("Live SPC Chart"),
        dccInterval(
          id = "interval-component",
          interval = 2 * 1000,  # in milliseconds
          n_intervals = 0,
          disabled = TRUE
        ),
        dccStore(
          id = "control-chart-state"
        ),
        dccGraph(
          id = "control-chart-live",
          figure = list(
            "data" = list(
              list(
                "type" = "scatter",
                "x" = list(),
                "y" = list(),
                "mode" = "lines+markers",
                "name" = params[[2]]
              )
            ),
            "layout" = list(
              "paper_bgcolor" = "#1d202d",
              "plot_bgcolor" = "#1d202d",
              "legend" = list("font" = list("color" = "darkgray")),
              "font" = list("color" = "darkgray"),
              "autosize" = TRUE
            )
          )
        )
      )
    )   
  )
}

build_tab_2 <- function() {
  return(
    htmlDiv(
      id = "status-container",
      children = list(
        build_quick_stats_panel(),
        htmlDiv(
          id = "graphs-container",
          children = list(
            build_top_panel(),
            build_chart_panel()
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
          #build_tabs(),
          htmlDiv(
            id = "app-content", 
            className = "container scalable",
            children = build_tab_2()
          )
          # htmlButton(
          #   children = "Proceed to Measurement", 
          #   id = "tab-trigger-btn", 
          #   n_clicks = 0
          # )
        )
      ),
      dccStore(
        id = "value-setter-store",
        data = init_value_setter_store()
      ),
      generate_modal()
    )
  )
)

########################################################################################################################
# DEFINE APP CALLBACKS

# callbacks for modal popup
app$callback(
  output = list(id = "markdown", property = "style"),
  params = list(
    input(id = "learn-more-button", property = "n_clicks"),
    input(id = "markdown_close", property = "n_clicks")
  ),
  function(button_click, close_click) {
    if (button_click > close_click) {
      return(list("display" = "block"))
    }
    return(list("display" = "none"))
  }
)

########################################################################################################################
# RUN SERVER

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}

