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

print(params) #delete
print(params[2:length(params)]) #delete
print(max_length) #delete

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
  ret <- list()
  for (i in 1:max_length) {
    if (data[[i]] >= ucl | data[[i]] <= lcl) {
      ooc_count <- ooc_count + 1
    }
    ret[[length(ret)+1]] <- ooc_count / (i+1)
  }
  return(ret)
}

# not done
init_df <- function() {
  return('')
}

state_dict <- init_df()

init_value_setter_store <- function() {
  return (init_df())
}

build_tab_1 <- function() {
  tmp_list <- lapply(params[2:length(params)], function(x) {return (list(label = x, value = x))})
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
                options = tmp_list,
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
                children = "Update", 
                id = "value-setter-set-btn"
                ),
              htmlButton(
                children = "View current setup", 
                id = "value-setter-view-btn", 
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

build_tab_2 <- function() {
  return(
    htmlDiv(
      id = "status-container",
      children = list(
        build_quick_stats_panel(),
        html.Div(
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
        ),
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
                **What does this app shows**
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
# DEFINE FUNCTIONS FOR APP LAYOUT AND CALLBACKS (ORGANIZE!)

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

generate_piechart <- function() {
  return(
    dccGraph(
      id = "piechart",
      figure = list(
        "data" = list(
          list(
            "labels" = params[2:length(params)],
            "values" = list(1, 1, 1, 1, 1, 1, 1),
            "type" = "pie",
            "marker" = list("line" = list("color" = "white", "width" = 1)),
            "hoverinfo" = "label",
            "textinfo" = "label"
          )
        ),
        "layout" = list(
          "showlegend" = TRUE,
          "paper_bgcolor" = "#1d202d",
          "plot_bgcolor" = "#1d202d",
          "font" = list("color" = "white"),
          "autosize" = TRUE
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
  item <- params[[index]]
  div_id <- glue("{item}{suffix_row}")
  button_id <- glue("{item}{suffix_button_id}")
  sparkline_graph_id <- glue("{item}{suffix_sparkline_graph}")
  count_id <- glue("{item}{suffix_count}")
  ooc_percentage_id <- glue("{item}{suffix_ooc_n}")
  ooc_graph_id <- glue("{item}{suffix_ooc_g}")
  indicator_id <- glue("{item}{suffix_indicator}")
  return(
    generate_metric_row(
      div_id,
      NA,
      list(
        "id" = item,
        "className" = "metric-row-button-text",
        "children" = html.Button(
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
        "id" = item + "_sparkline",
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
          figure = plot_ly(
            df,
            "x" = list(),
            "y" = list(),
            "mode" = "lines+markers",
            "name" = item,
            "line" = list("color" = "#f4d44d")
          ) %>% layout(
            "uirevision" = TRUE,
            "margin" = list(l=0, r=0, t=4, b=4, pad=0),
            "paper_bgcolor" = "#1d202d",
            "plot_bgcolor" = "#1d202d"
          )
          # figure = go.Figure(
          #   list(
          #     "data" = list(
          #       list(
          #         "x" = list(),
          #         "y" = list(),
          #         "mode" = "lines+markers",
          #         "name" = item,
          #         "line" = list("color" = "#f4d44d"),
          #       )
          #     ),
          #     "layout" = list(
          #       "uirevision" = TRUE,
          #       "margin" = list(l=0, r=0, t=4, b=4, pad=0),
          #       "paper_bgcolor" = "#1d202d",
          #       "plot_bgcolor" = "#1d202d",
          #     )
          #   )
          # )
        )
      ),
      list(
        "id" = ooc_percentage_id,
        "children" = "0.00%"
      ),
      list(
        "id" = ooc_graph_id + "_container",
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
          max = 15,
          value = 0
        )
      ),
      list(
        "id" = item + "_pf",
        "children" = daqIndicator(
          id = indicator_id, 
          value = TRUE, 
          color = "#91dfd2"
        )
      )
    )
  )
}

generate_metric_row <- function(id, style, col1, col2, col3, col4, col5, col6) {
  if (is.na(style)) {
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

# not done
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
          figure = plot_ly(
            df,
            "x" = list(),
            "y" = list(),
            "mode" = "lines+markers",
            "name" = params[[2]]
          ) %>% layout(
            "paper_bgcolor" = "#1d202d",
            "plot_bgcolor" = "#1d202d",
            "autosize" = TRUE
          )
          # figure = go.Figure(
          #   list(
          #     "data" = list(
          #       list(
          #         "x" = list(),
          #         "y" = list(),
          #         "mode" = "lines+markers",
          #         "name" = params[[2]]
          #       )
          #     ),
          #     "layout" = list(
          #       "paper_bgcolor" = "#1d202d",
          #       "plot_bgcolor" = "#1d202d",
          #       "autosize" = TRUE
          #     )
          #   )
          # )
        )
      )
    )   
  )
}

# not done
generate_graph <- function(interval, specs_dict, col) {
  return(
    ""
  )
}

# not done
update_sparkline <- function(interval, param) {
  return(
    ""
  )
}

# not done
update_count <- function(interval, col, data) {
  return(
    ""
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
      ),
      generate_modal()
    )
  )
)

########################################################################################################################
# DEFINE APP CALLBACKS

# callback for switching tabs
app$callback(
  output = list(id = "app-tabs", property = "value"),
  params = list(input(id = "tab-trigger-btn", property = "n_clicks")),
  function(tab_switch) {
    if (tab_switch == 0) {
      return("tab1")
    }
    if (tab_switch) {
      return("tab2")
    }
  }
)
app$callback(
  output = list(id = "app-content", property = "children"),
  params = list(input(id = "tab-trigger-btn", property = "n_clicks")),
  function(tab_switch) {
    if (tab_switch == 0) {
      return(build_tab_1())
    }
    if (tab_switch) {
      return(build_tab_2())
    }
  }
)
app$callback(
  output = list(id = "Specs-tab", property = "disabled"),
  params = list(input(id = "tab-trigger-btn", property = "n_clicks")),
  function(tab_switch) {
    if (tab_switch == 0) {
      return(FALSE)
    }
    if (tab_switch) {
      return(TRUE)
    }
  }
)
app$callback(
  output = list(id = "Control-chart-tab", property = "disabled"),
  params = list(input(id = "tab-trigger-btn", property = "n_clicks")),
  function(tab_switch) {
    if (tab_switch == 0) {
      return(TRUE)
    }
    if (tab_switch) {
      return(FALSE)
    }
  }
)
app$callback(
  output = list(id = "tab-trigger-btn", property = "style"),
  params = list(
    input(id = "tab-trigger-btn", property = "n_clicks")
  ),
  function(tab_switch) {
    if (tab_switch == 0) {
      return(list("display" = "inline-block", "float" = "right"))
    }
    if (tab_switch) {
      return(list("display" = "none"))
    }
  }
)

# callback for modal popup
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
