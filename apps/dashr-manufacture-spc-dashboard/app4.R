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
# DEFINE GLOBAL VARIABLES

params = names(df)
max_length = nrow(df)

#print(params) # DELETE!
#print(params[-1]) # DELETE!
#print(max_length) # DELETE!

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

# GENERATE TAB 1 LAYOUT

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

state_dict_all <- init_df()
state_dict <- init_df()[-1]

init_value_setter_store <- function() {
  state_dict <- init_df()[-1]
  return(state_dict)
}

ud_usl_input <- daqNumericInput(
  id = "ud_usl_input",
  className = "setting-input",
  #value = state_dict[["Diameter"]][["usl"]], # DELETE!
  size = 200,
  max = 9999999
)
ud_lsl_input <- daqNumericInput(
  id = "ud_lsl_input",
  className = "setting-input",
  #value = state_dict[["Diameter"]][["lsl"]], # DELETE!
  size = 200,
  max = 9999999
)
ud_ucl_input <- daqNumericInput(
  id = "ud_ucl_input",
  className = "setting-input",
  #value = state_dict[["Diameter"]][["ucl"]], # DELETE!
  size = 200,
  max = 9999999
)
ud_lcl_input <- daqNumericInput(
  id = "ud_lcl_input",
  className = "setting-input",
  #value = state_dict[["Diameter"]][["lcl"]], # DELETE!
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
        )
      )
    )
  )
}

build_tab_1 <- function() {
  dropdown_list <- lapply(params[-1], function(x) {return(list(label = x, value = x))})
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
                options = dropdown_list,
                value = params[[2]]
              )
            )
          ),
          htmlDiv(
            id = "value-setter-menu",
            children = list(
              htmlDiv(
                id = "value-setter-panel"
                #children = list("HERE1")
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
                #children = list("HERE2")
              )
            )
          )
        )
      )
    )
  )
}

# GENERATE TAB 2 LAYOUT

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
            "values" = list(1, 1, 1, 1, 1, 1, 1), # CHANGE TO: "values" = list(),
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
          #"uirevision" = TRUE,
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
          disabled = FALSE # CHANGE BACK: disabled = TRUE
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

generate_graph <- function(interval, specs_dict, col) {
  mean <- state_dict[[col]][["mean"]]
  ucl <- specs_dict[[col]][["ucl"]]
  lcl <- specs_dict[[col]][["lcl"]]
  usl <- specs_dict[[col]][["usl"]]
  lsl <- specs_dict[[col]][["lsl"]]
  x_array <- as.list(state_dict_all[["Batch"]][["data"]])
  y_array <- as.list(state_dict[[col]][["data"]])
  total_count <- 0
  if (interval > max_length) {
    total_count <- max_length - 1
  } else if (interval > 0) {
    total_count <- interval
  }
  x_data <- x_array[1:total_count]
  y_data <- y_array[1:total_count]
  ooc_trace <- list(
    "x": list(),
    "y": list(),
    "name": "Out of Control",
    "mode": "markers",
    "marker": list(color = "rgba(210, 77, 87, 0.7)", symbol = "square", size = 11)
  )
  for (i in 1:length(y_data)) {
    if (y_data[[i]] >= ucl | y_data[[i]] <= lcl) {
      ooc_trace[["x"]][[length(ooc_trace[["x"]])+1]] <- i
      ooc_trace[["y"]][[length(ooc_trace[["y"]])+1]] <- y_data[[i]]
    }
  }
  histo_trace <- list(
    "x" = x_data,
    "y" = y_data,
    "type" = "histogram",
    "orientation" = "h",
    "name" = "Distribution",
    "xaxis" = "x2",
    "yaxis" = "y2",
    "marker" = list("color" = "#f4d44d"),
  )
  fig <- list(
    "data" = list(
      list(
        "x" = x_data,
        "y" = y_data,
        "mode" = "lines+markers",
        "name" = col,
        "line" = list("color" = "#f4d44d"),
      ),
      ooc_trace,
      histo_trace,
    )
  )
  len_figure <- len(fig[["data"]][[0]][["x"]])
  fig[["layout"]] <- list(
    hovermode = "closest",
    uirevision = col,
    paper_bgcolor = "#1d202d",
    plot_bgcolor = "#1d202d",
    legend = list("font" = list("color" = "darkgray")),
    font = list("color" = "darkgray"),
    showlegend = TRUE,
    xaxis = list(
      "zeroline" = FALSE,
      "title" = "Batch_Num",
      "showline" = FALSE,
      "domain" = list(0, 0.8),
      "titlefont" = list("color" = "darkgray")
    ),
    yaxis = list(
      "title" = col, 
      "autorange" = TRUE, 
      "titlefont" = list("color" = "darkgray")
    ),
    annotations = list(
      list(
        "x" = 0.75,
        "y" = lcl,
        "xref" = "paper",
        "yref" = "y",
        #"text" = "LCL:" + str(round(lcl, 3)),
        "showarrow" = FALSE,
        "font" = list("color" = "white"),
      ),
      list(
        "x" = 0.75,
        "y" = ucl,
        "xref" = "paper",
        "yref" = "y",
        #"text" = "UCL:" + str(round(ucl, 3)),
        "showarrow" = FALSE,
        "font" = list("color" = "white"),
      ),
      list(
        "x" = 0.75,
        "y" = usl,
        "xref" = "paper",
        "yref" = "y",
        #"text" = "USL:" + str(round(usl, 3)),
        "showarrow" = FALSE,
        "font" = list("color" = "white"),
      ),
      list(
        "x" = 0.75,
        "y" = lsl,
        "xref" = "paper",
        "yref" = "y",
        #"text" = "LSL:" + str(round(lsl, 3)),
        "showarrow" = FALSE,
        "font" = list("color" = "white"),
      ),
      list(
        "x" = 0.75,
        "y" = mean,
        "xref" = "paper",
        "yref" = "y",
        #"text" = "Targeted mean:" + str(round(mean, 3)),
        "showarrow" = FALSE,
        "font" = list("color" = "white"),
      ),
    ),
    shapes = list(
      list(
        "type" = "line",
        "xref" = "x",
        "yref" = "y",
        "x0" = 1,
        "y0" = usl,
        "x1" = len_figure + 1,
        "y1" = usl,
        "line" = list("color" = "#91dfd2", "width" = 1, "dash" = "dashdot"),
      ),
      list(
        "type" = "line",
        "xref" = "x",
        "yref" = "y",
        "x0" = 1,
        "y0" = lsl,
        "x1" = len_figure + 1,
        "y1" = lsl,
        "line" = list("color" = "#91dfd2", "width" = 1, "dash" = "dashdot"),
      ),
      list(
        "type" = "line",
        "xref" = "x",
        "yref" = "y",
        "x0" = 1,
        "y0" = ucl,
        "x1" = len_figure + 1,
        "y1" = ucl,
        "line" = list("color" = "rgb(255,127,80)", "width" = 1, "dash" = "dashdot"),
      ),
      list(
        "type" = "line",
        "xref" = "x",
        "yref" = "y",
        "x0" = 1,
        "y0" = mean,
        "x1" = len_figure + 1,
        "y1" = mean,
        "line" = list("color" = "rgb(255,127,80)", "width" = 2),
      ),
      list(
        "type" = "line",
        "xref" = "x",
        "yref" = "y",
        "x0" = 1,
        "y0" = lcl,
        "x1" = len_figure + 1,
        "y1" = lcl,
        "line" = list("color" = "rgb(255,127,80)", "width" = 1, "dash" = "dashdot"),
      ),
    ),
    xaxis2 = list(
      "title" = "count",
      "domain" = list(0.8, 1),  # 70 to 100 % of width
      "titlefont" = list("color" = "darkgray"),
    ),
    yaxis2 = list(
      "anchor" = "free",
      "overlaying" = "y",
      "side" = "right",
      "showticklabels" = FALSE,
      "titlefont" = list("color" = "darkgray")
    )
  )
  return(fig)
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

# callbacks for modal popup when clicking "Learn More"
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

# DEFINE TAB 1 CALLBACKS

# callbacks to update values based on stored data and dropdown selection
app$callback(
  output = list(id = "ud_usl_input", property = "value"),
  params = list(
    input(id = "metric-select-dropdown", property = "value"),
    state(id = "value-setter-store", property = "data")
  ),
  function(dd_select, state_value) {
    return(state_value[[dd_select]][["usl"]])
  }
)
app$callback(
  output = list(id = "ud_lsl_input", property = "value"),
  params = list(
    input(id = "metric-select-dropdown", property = "value"),
    state(id = "value-setter-store", property = "data")
  ),
  function(dd_select, state_value) {
    return(state_value[[dd_select]][["lsl"]])
  }
)
app$callback(
  output = list(id = "ud_ucl_input", property = "value"),
  params = list(
    input(id = "metric-select-dropdown", property = "value"),
    state(id = "value-setter-store", property = "data")
  ),
  function(dd_select, state_value) {
    return(state_value[[dd_select]][["ucl"]])
  }
)
app$callback(
  output = list(id = "ud_lcl_input", property = "value"),
  params = list(
    input(id = "metric-select-dropdown", property = "value"),
    state(id = "value-setter-store", property = "data")
  ),
  function(dd_select, state_value) {
    return(state_value[[dd_select]][["lcl"]])
  }
)
app$callback(
  output = list(id = "value-setter-panel", property = "children"),
  params = list(input(id = "metric-select-dropdown", property = "value")),
  function(dd_select) {
    return(
      list(
        build_value_setter_line(
          "value-setter-panel-header",
          "Specs",
          "Historical Value",
          "Set new value"
        ),
        build_value_setter_line(
          "value-setter-panel-usl",
          "Upper Specification limit",
          state_dict[[dd_select]][["usl"]],
          ud_usl_input
        ),
        build_value_setter_line(
          "value-setter-panel-lsl",
          "Lower Specification limit",
          state_dict[[dd_select]][["lsl"]],
          ud_lsl_input
        ),
        build_value_setter_line(
          "value-setter-panel-ucl",
          "Upper Control limit",
          state_dict[[dd_select]][["ucl"]],
          ud_ucl_input
        ),
        build_value_setter_line(
          "value-setter-panel-lcl",
          "Lower Control limit",
          state_dict[[dd_select]][["lcl"]],
          ud_lcl_input
        )
      )
    )
  }
)

# callbacks to update stored data when clicking "Update" button
app$callback(
  output = list(id = "value-setter-store", property = "data"),
  params = list(
    input(id = "value-setter-set-btn", property = "n_clicks"),
    state(id = "metric-select-dropdown", property = "value"),
    state(id = "value-setter-store", property = "data"),
    state(id = "ud_usl_input", property = "value"),
    state(id = "ud_lsl_input", property = "value"),
    state(id = "ud_ucl_input", property = "value"),
    state(id = "ud_lcl_input", property = "value")
  ),
  function(set_btn, param, stats, usl, lsl, ucl, lcl) {
    if (is.integer(set_btn)) {
      stats[[param]][["usl"]] <- usl
      stats[[param]][["lsl"]] <- lsl
      stats[[param]][["ucl"]] <- ucl
      stats[[param]][["lcl"]] <- lcl
      stats[[param]][["ooc"]] <- populate_ooc(df[[param]], ucl, lcl)
    }
    return(stats)
  }
)

# callbacks to show current stored data when clicking "View Current Setup" button
app$callback(
  output = list(id = "value-setter-view-output", property = "children"),
  params = list(
    input(id = "value-setter-view-btn", property = "n_clicks"),
    input(id = "metric-select-dropdown", property = "value"),
    input(id = "value-setter-store", property = "data")
  ),
  function(n_clicks, dd_select, store_data) {
    if (n_clicks > 0) {
      new_df <- data.frame(
        "Specs" = c(
          "Upper Specification Limit",
          "Lower Specification Limit",
          "Upper Control Limit",
          "Lower Control Limit"
        ),
        "Current Setup" = c(
          store_data[[dd_select]][["usl"]],
          store_data[[dd_select]][["lsl"]],
          store_data[[dd_select]][["ucl"]],
          store_data[[dd_select]][["lcl"]]
        )
      )
      return(
        dashTable::dashDataTable(
          style_header = list(
            "backgroundColor" = "#2d3038", 
            "fontWeight" = "bold"
          ),
          style_as_list_view = TRUE,
          style_cell_conditional = lapply(list("Specs"), function(x) {return(list("if" = list("column_id" = x), "textAlign" = "left"))}),
          style_cell = list(
            "backgroundColor" = "#2d3038",
            "color" = "darkgray",
            "border" = "darkgray"
          ),
          data = dashTable::df_to_list(new_df),
          columns = lapply(list("Specs", "Current Setup"), function(x) {return(list("id" = x, "name" = x))})
        )
      )
    }
  }
)

# callbacks for switching tabs when clicking "Proceed To Measurement" button
app$callback(
  output = list(id = "app-tabs", property = "value"),
  params = list(input(id = "tab-trigger-btn", property = "n_clicks")),
  function(tab_switch) {
    if (tab_switch == 0) {
      return("tab1")
    } else if (tab_switch) {
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
    } else if (tab_switch) {
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
    } else if (tab_switch) {
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
    } else if (tab_switch) {
      return(FALSE)
    }
  }
)
app$callback(
  output = list(id = "tab-trigger-btn", property = "style"),
  params = list(input(id = "tab-trigger-btn", property = "n_clicks")),
  function(tab_switch) {
    if (tab_switch == 0) {
      return(list("display" = "inline-block", "float" = "right"))
    } else {
      return(list("display" = "none"))
    }
  }
)

# DEFINE TAB 2 CALLBACKS

# callbacks for stopping interval update when clicking "Stop/Start" button
app$callback(
  output = list(id = "interval-component", property = "disabled"),
  params = list(
    input(id = "stop-button", property = "n_clicks"),
    state(id = "interval-component", property = "disabled")
  ),
  function(n_clicks, current) {
    return(!current)
  }
)
app$callback(
  output = list(id = "stop-button", property = "buttonText"),
  params = list(
    input(id = "stop-button", property = "n_clicks"),
    state(id = "interval-component", property = "disabled")
  ),
  function(n_clicks, current) {
    ifelse(current, text <- "stop", text <- "start")
    return(text)
  }
)

# callbacks to update progress gauge i.e. "Time to completion"
app$callback(
  output = list(id = "progress-gauge", property = "value"),
  params = list(input(id = "interval-component", property = "n_intervals")),
  function(interval) {
    if (interval < max_length) {
      total_count <- interval
    } else {
      total_count <- max_length
    }
    return(as.integer(total_count))
  }
)

########################################################################################################################
# RUN SERVER

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
