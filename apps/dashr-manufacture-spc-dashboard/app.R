library(rlist)
library(stringr)
library(plotly)
library(dash)
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
# DEFINE GLOBAL VARIABLES

theme <- list(              # old      new
  "dark1" = "#1E2130",      # #2D3038  #1E2130
  "dark2" = "#171A29",      # #1D202D  #171A29
  "gray" = "darkgray",      # darkgray #A7ADB9
  "white" = "#FFFFFF",      # #FFFFFF  #F9FEFE
  "blue" = "#91DFD2",       # #91DFD2  #90E0D4
  "yellow" = "#F4D44D",     # #F4D44D  #F5D54C
  "red" = "#F45060"         # #F45060  #F45060
)

suffix <- list(
  row = "_row",
  button_id = "_button",
  count = "_count",
  sparkline_id = "_sparkline_id",
  sparkline_graph = "_sparkline_graph",
  ooc_n = "_OOC_number",
  ooc_g = "_OOC_graph",
  indicator = "_indicator"
)

########################################################################################################################
# HANDLE & STORE DATA

df <- read.table(paste("data", "spc_data.csv", sep="/"), header = TRUE, sep = ",")

names(df) <- gsub("[.]", "-", names(df))
params = names(df) # list of parameters (params[-1] removes "Batch")
max_length = nrow(df) # number of data entries

populate_ooc <- function(data, ucl, lcl) {
  ooc_count <- 0
  output <- list()
  for (i in 1:max_length) {
    if (data[[i]] >= ucl | data[[i]] <= lcl) {
      ooc_count <- ooc_count + 1
    }
    output[[i]] <- ooc_count / (i+1)
  }
  return(output)
}

init_df <- function() {
  output <- list()
  for (param in params) {
    data <- df[[param]]
    stats <- summary(data)
    ucl <- stats[["Mean"]]+3*sd(data)
    lcl <- stats[["Mean"]]-3*sd(data)
    usl <- stats[["Mean"]]+sd(data)
    lsl <- stats[["Mean"]]-sd(data)
    output[[param]] <- list(
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
init_value_setter_store <- function() {return(init_df())} # TODO: check if this needs to called this way

########################################################################################################################
# DEFINE FUNCTIONS FOR APP LAYOUT AND CALLBACKS

generate_modal <- function() {
  return(
    htmlDiv(
      id = "markdown",
      className = "modal",
      children = list(
        htmlDiv(
          className = "markdown-container",
          children = list(
            htmlDiv(
              className = "close-container",
              children = htmlButton(
                id = "close-button",
                className = "button",
                children = "Close",
                n_clicks = 0
              )
            ),
            htmlDiv(
              className = "markdown-text",
              children = dccMarkdown("
              **What is this mock app about?**

              'dash-manufacture-spc-dashboard' is a dashboard for monitoring read-time process quality along manufacture production line.

              **What does this app show?**

              Click on buttons in 'Parameter' column to visualize details of measurement trendlines on the bottom panel.

              The Sparkline on top panel and Control chart on bottom panel show Shewhart process monitor using mock data.
              The trend is updated every other second to simulate real-time measurements. Data falling outside of six-sigma control
              limit are signals indicating 'Out of Control(OOC)', and will trigger alerts instantly for a detailed checkup.
              ")
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
              id = "learn-more-button",
              children = "LEARN MORE",
              n_clicks = 0
            ),
            htmlImg(
              id = "logo",
              src = paste("assets", "dash-logo.png", sep="/")
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
      children = list(
        dccTabs(
          id = "app-tabs",
          value = "tab2",
          children = list(
            dccTab(
              id = "specs-tab",
              label = "Specification Settings",
              value = "tab1",
              className = "tab",
              selected_className = "tab--selected"
            ),
            dccTab(
              id = "control-chart-tab",
              label = "Control Charts Dashboard",
              value = "tab2",
              className = "tab",
              selected_className = "tab--selected"
            )
          )
        )
      )
    )
  )
}

# TAB 1 ----------------------------------------------------------------------------------------------------------------

build_tab_1 <- function() {
  return(
    htmlDiv(
      id = "specs-container",
      children = list(
        htmlP("Use historical control limits to establish a benchmark, or set new values."),
        htmlDiv(
          id = "settings-menu",
          children = list(
          htmlDiv(
              id = "metric-select-menu",
              children = list(
                htmlDiv(
                  className = "row",
                  children = htmlLabel(
                    children = htmlH6("Select Metrics")
                  )
                ),
                htmlDiv(
                  className = "row",
                  children = htmlDiv(
                    children = dccDropdown(
                      id = "metric-select-dropdown",
                      options = lapply(params[-1], function(x) {return(list(label = x, value = x))}),
                      value = params[[2]]
                    )
                  )
                )
              )
            ),
            htmlDiv(
              id = "value-setter-menu",
              children = list(
                htmlDiv(id = "value-setter-panel"),
                htmlBr(),
                htmlDiv(
                  id = "value-setter-buttons",
                  children = list(
                    htmlButton(
                      id = "value-setter-set-btn",
                      children = "Update"
                    ),
                    htmlButton(
                      id = "value-setter-view-btn",
                      children = "View current setup",
                      n_clicks = 0
                    )
                  )
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
  )
}

build_value_setter_line <- function(label, hist_value, new_value) {
  return(
    htmlDiv(
      className = "row",
      children = list(
        htmlLabel(
          className = "four columns",
          children = label
        ),
        htmlLabel(
          className = "four columns",
          children = hist_value
        ),
        htmlDiv(
          className = "four columns",
          children = new_value
        )
      )
    )
  )
}

# TAB 2 ----------------------------------------------------------------------------------------------------------------

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

build_quick_stats_panel <- function() {
  return(
    htmlDiv(
      id = "quick-stats",
      className = "row",
      children = list(
        htmlDiv(
          id = "card-1",
          children = list(
            htmlH6("Operator ID"),
            daqLEDDisplay(
              id = "operator-led",
              value = "1704",
              color = theme$blue,
              backgroundColor = theme$dark1,
              size = "50"
            )
          )
        ),
        htmlDiv(
          id = "card-2",
          children = list(
            htmlH6("Time to completion"),
            daqGauge(
              id = "progress-gauge",
              value = 0,
              max = max_length,
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
      children = htmlH6(title)
    )
  )
}

generate_metric_row <- function(id, style, col1, col2, col3, col4, col5, col6) {
  if (length(style) == 0) {
    style <- list(
      height = "10rem",
      width = "100%"
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
      id = "metric_header",
      style = list(
        height = "3rem",
        margin = "1rem 0",
        textAlign = "center"
      ),
      list(
        id = "m_header_1",
        children = htmlDiv("Parameter")
      ),
      list(
        id = "m_header_2",
        children = htmlDiv("Count")
      ),
      list(
        id = "m_header_3",
        children = htmlDiv("Sparkline")
      ),
      list(
        id = "m_header_4",
        children = htmlDiv("OOC%")
      ),
      list(
        id = "m_header_5",
        children = htmlDiv("%OOC")
      ),
      list(
        id = "m_header_6",
        children = "Pass/Fail"
      )
    )
  )
}

generate_metric_row_helper <- function(index) {
  param <- params[-1][[index]]
  ooc_colors <- list(list(0, 3), list(3, 7), list(7, 15))
  names(ooc_colors) <- list(theme$blue, theme$yellow, theme$red)
  return(
    generate_metric_row(
      id = sprintf("%s%s", param, suffix$row),
      style = list(),
      list(
        id = param,
        className = "metric-row-button-text",
        children = htmlButton(
          id = sprintf("%s%s", param, suffix$button_id),
          className = "metric-row-button",
          children = param,
          title = "Click to visualize live SPC chart",
          n_clicks = 0
        )
      ),
      list(
        id = sprintf("%s%s", param, suffix$count),
        children = "0"
      ),
      list(
        id = sprintf("%s%s", param, suffix$sparkline_id),
        children = dccGraph(
          id = sprintf("%s%s", param, suffix$sparkline_graph),
          style = list(
            width = "100%",
            height = "95%"
          ),
          config = list(
            staticPlot = FALSE,
            editable = FALSE,
            displayModeBar = FALSE
          ),
          figure = list(
            data = list(
              list(
                type = "scatter",
                x = list(),
                y = list(),
                mode = "lines+markers",
                name = param,
                line = list("color" = theme$yellow)
              )
            ),
            layout = list(
              uirevision = TRUE,
              margin = list(l=0, r=0, t=4, b=4, pad=0),
              paper_bgcolor = theme$dark2,
              plot_bgcolor = theme$dark2
            )
          )
        )
      ),
      list(
        id = sprintf("%s%s", param, suffix$ooc_n),
        children = "0.00%"
      ),
      list(
        id = sprintf("%s%s", param, "_OOC_graph_container"),
        children = daqGraduatedBar(
          id = sprintf("%s%s", param, suffix$ooc_g),
          color = list(ranges = ooc_colors),
          showCurrentValue = FALSE,
          #size = 146.16,
          max = 15,
          value = 0
        )
      ),
      list(
        id = sprintf("%s%s", param, "_pf"),
        children = daqIndicator(
          id = sprintf("%s%s", param, suffix$indicator),
          value = TRUE,
          color = theme$blue
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
            "labels" = params[-1],
            "values" = lapply(1:7, function(x) {return(1)}),
            "type" = "pie",
            "marker" = list(
              "colors" = lapply(1:7, function(x) {return(theme$blue)}),
              "line" = list("color" = "white", "width" = 2)
            ),
            "hoverinfo" = "label",
            "textinfo" = "label"
          )
        ),
        "layout" = list(
          "font" = list("color" = "white"),
          "paper_bgcolor" = theme$dark2,
          "plot_bgcolor" = theme$dark2,
          "margin" = list(t = 40, b = 40, r = 60, l = 60), # TODO: edit margins
          "uirevision" = TRUE,
          "showlegend" = FALSE,
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

build_chart_panel <- function() {
  return(
    htmlDiv(
      id = "control-chart-container",
      className = "twelve columns",
      children = list(
        generate_section_banner("Live SPC Chart"),
        dccGraph(
          id = "control-chart-live",
          figure = generate_graph(0, state_dict, params[-1][[1]])
        )
      )
    )
  )
}

generate_graph <- function(interval, specs_dict, col) {
  stats <- list(
    mean = list(data = state_dict[[col]][["mean"]], label = "Targeted Mean", color = "rgb(255,127,80)", style = "solid", width = 2),
    ucl = list(data = specs_dict[[col]][["ucl"]], label = "ucl", color = "rgb(255,127,80)", style = "dot", width = 1),
    lcl = list(data = specs_dict[[col]][["lcl"]], label = "lcl", color = "rgb(255,127,80)", style = "dot", width = 1),
    usl = list(data = specs_dict[[col]][["usl"]], label = "usl", color = theme$blue, style = "dot", width = 1),
    lsl = list(data = specs_dict[[col]][["lsl"]], label = "lsl", color = theme$blue, style = "dot", width = 1)
  )
  total_count <- ifelse(interval < max_length, interval, max_length)
  x_data <- list()
  y_data <- list()
  x_data_ooc <- list()
  y_data_ooc <- list()
  if (total_count != 0) {
    x_data <- state_dict[["Batch"]][["data"]][1:total_count]
    y_data <- state_dict[[col]][["data"]][1:total_count]
    x_data_ooc <- lapply(1:total_count, function(i) {return(ifelse(y_data[[i]] > stats$lcl$data && y_data[[i]] < stats$ucl$data, NA, i))})
    y_data_ooc <- lapply(1:total_count, function(i) {return(ifelse(y_data[[i]] > stats$lcl$data && y_data[[i]] < stats$ucl$data, NA, y_data[[i]]))})
  }
  # define figure
  fig <- list(
    data = list(
      # parameter trace
      list(
        type = "scatter",
        x = x_data,
        y = y_data,
        mode = "lines+markers",
        name = col,
        line = list("color" = theme$yellow)
      ),
      # ooc trace
      list(
        type = "scatter",
        x = x_data_ooc,
        y = y_data_ooc,
        mode = "markers",
        name = "Out of Control",
        marker = list(color = theme$red, symbol = "square", size = 11)
      ),
      # histogram
      list(
        type = "histogram",
        x = x_data,
        y = y_data,
        orientation = "h",
        xaxis = "x2",
        yaxis = "y2",
        name = "Distribution",
        marker = list("color" = theme$yellow)
      )
    ),
    layout = list(
      hovermode = "closest",
      uirevision = col,
      paper_bgcolor = theme$dark2,
      plot_bgcolor = theme$dark2,
      margin = list(t = 40), # TODO: fix margins
      legend = list("font" = list("color" = theme$gray)),
      font = list("color" = theme$gray),
      showlegend = TRUE,
      shapes = lapply(
        unname(stats),
        function(line) {
          return(
            list(
              type = "line",
              xref = "x",
              yref = "y",
              x0 = 1,
              y0 = line$data,
              x1 = total_count+1,
              y1 = line$data,
              line = list(color = line$color, width = line$width, dash = line$style)
            )
          )
        }
      ),
      annotations = lapply(
        unname(stats),
        function(line) {
          return(
            list(
              x = 0.75,
              y = line$data,
              xref = "paper",
              yref = "y",
              text = sprintf("%s: %.3f", ifelse(line$label == "Targeted Mean", line$label, toupper(line$label)), line$data),
              showarrow = FALSE,
              font = list("color" = "white")
            )
          )
        }
      ),
      xaxis = list(
        zeroline = FALSE,
        title = "Batch_Num",
        showline = FALSE,
        domain = list(0, 0.8),
        titlefont = list("color" = theme$gray)
      ),
      yaxis = list(
        title = col,
        autorange = TRUE,
        titlefont = list("color" = theme$gray)
      ),
      xaxis2 = list(
        title = "Count",
        domain = list(0.8, 1),  # 70 to 100 % of width
        titlefont = list("color" = theme$gray)
      ),
      yaxis2 = list(
        anchor = "free",
        overlaying = "y",
        side = "right",
        showticklabels = FALSE,
        titlefont = list("color" = theme$gray)
      )
    )
  )
  return(fig)
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
          htmlDiv(id = "app-content")
        )
      ),
      dccInterval(
        id = "interval-component",
        interval = 2*1000, # in milliseconds
        n_intervals = 10, # TODO: change this to 50 as final, to 10 for trials
        disabled = TRUE
      ),
      dccStore(
        id = "interval-component-store",
        data = 10 # TODO: change this to 50 as final, to 10 for trials
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
    input(id = "close-button", property = "n_clicks")
  ),
  function(button_click, close_click) {
    if (button_click > close_click) {
      return(list("display" = "block"))
    }
    return(list("display" = "none"))
  }
)

# callbacks for switching tabs
app$callback(
  output = list(id = "app-content", property = "children"),
  params = list(input(id = "app-tabs", property = "value")),
  function(tab) {
    if (tab == "tab1") {
      return(build_tab_1())
    } else if (tab == "tab2") {
      return(build_tab_2())
    }
  }
)
app$callback(
  output = list(id = "interval-component", property = "n_intervals"),
  params = list(
    input(id = "app-tabs", property = "value"),
    state(id = "interval-component-store", property = "data")
  ),
  function(tab, stopped_interval) {
    return(stopped_interval)
  }
)
app$callback(
  output = list(id = "interval-component-store", property = "data"),
  params = list(
    input(id = "app-tabs", property = "value"),
    state(id = "interval-component", property = "disabled"),
    state(id = "interval-component", property = "n_intervals"),
    state(id = "interval-component-store", property = "data")
  ),
  function(tab, disabled, cur_interval, cur_stage) {
    if (disabled) {
      return(cur_interval)
    }
    if (tab == "tab1") {
      return(cur_interval)
    }
    return(cur_stage)
    # TODO: change back to the following:
    # if (disabled) {
    #   return(cur_interval)
    # } else if (tab == "tab1") {
    #   return(cur_interval)
    # } else if (tab == "tab2") {
    #   return(cur_stage)
    # }
  }
)

# TAB 1 ----------------------------------------------------------------------------------------------------------------

# callbacks to update values based on stored data and dropdown selection
app$callback(
  output = list(id = "value-setter-panel", property = "children"),
  params = list(
    input(id = "metric-select-dropdown", property = "value"),
    state(id = "value-setter-store", property = "data")
  ),
  function(dd_select, state_value) {
    return(
      list(
        build_value_setter_line(
          htmlH6("Specs"),
          htmlH6("Historical Value"),
          htmlH6("Set new value")
        ),
        build_value_setter_line(
          "Upper Specification Limit",
          state_dict[[dd_select]][["usl"]],
          daqNumericInput(
            id = "ud_usl_input",
            className = "setting-input",
            size = 200,
            max = 9999999,
            value = state_value[[dd_select]][["usl"]]
          )
        ),
        build_value_setter_line(
          "Lower Specification Limit",
          state_dict[[dd_select]][["lsl"]],
          daqNumericInput(
            id = "ud_lsl_input",
            className = "setting-input",
            size = 200,
            max = 9999999,
            value = state_value[[dd_select]][["lsl"]]
          )
        ),
        build_value_setter_line(
          "Upper Control Limit",
          state_dict[[dd_select]][["ucl"]],
          daqNumericInput(
            id = "ud_ucl_input",
            className = "setting-input",
            size = 200,
            max = 9999999,
            value = state_value[[dd_select]][["ucl"]]
          )
        ),
        build_value_setter_line(
          "Lower Control Limit",
          state_dict[[dd_select]][["lcl"]],
          daqNumericInput(
            id = "ud_lcl_input",
            className = "setting-input",
            size = 200,
            max = 9999999,
            value = state_value[[dd_select]][["lcl"]]
          )
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
        dashDataTable(
          style_header = list(
            "backgroundColor" = theme$dark1,
            "fontWeight" = "bold"
          ),
          style_as_list_view = TRUE,
          style_cell_conditional = lapply(list("Specs"), function(x) {return(list("if" = list("column_id" = x), "textAlign" = "left"))}),
          style_cell = list(
            "backgroundColor" = theme$dark1,
            "color" = theme$gray,
            "border" = theme$gray
          ),
          data = df_to_list(new_df),
          columns = lapply(colnames(new_df), function(x) {return(list("id" = x, "name" = x))})
        )
      )
    }
  }
)

# TAB 2 ----------------------------------------------------------------------------------------------------------------

# callbacks for pausing interval when clicking "Stop/Start" button
app$callback(
  output = list(id = "interval-component", property = "disabled"),
  params = list(
    input(id = "stop-button", property = "n_clicks"),
    state(id = "interval-component", property = "disabled")
  ),
  function(n_clicks, current) {
    return(ifelse(n_clicks == 0, TRUE, !current))
  }
)
app$callback(
  output = list(id = "stop-button", property = "buttonText"),
  params = list(
    input(id = "stop-button", property = "n_clicks"),
    state(id = "interval-component", property = "disabled")
  ),
  function(n_clicks, current) {
    return(ifelse(n_clicks == 0, "start", ifelse(current, "stop", "continue")))
  }
)

# callbacks to update progress gauge i.e. "Time to completion"
app$callback(
  output = list(id = "progress-gauge", property = "value"),
  params = list(input(id = "interval-component", property = "n_intervals")),
  function(interval) {
    total_count <- ifelse(interval < max_length, interval, max_length)
    return(as.integer(total_count))
  }
)

# callbacks to update each metric row according to interval
update_metric_summary <- function(param) {
  app$callback(
    output = list(id = sprintf("%s%s", param, suffix$count), property = "children"),
    params = list(input(id = "interval-component", property = "n_intervals")),
    function(interval, stored_data) {
      total_count <- ifelse(interval < max_length, interval, max_length)
      return(as.character(total_count))
    }
  )
  app$callback(
    output = list(id = sprintf("%s%s", param, suffix$sparkline_graph), property = "extendData"),
    params = list(input(id = "interval-component", property = "n_intervals")),
    function(interval) {
      total_count <- ifelse(interval < max_length, interval, max_length)
      x_new <- ifelse(total_count == 0, NULL, state_dict[["Batch"]][["data"]][[total_count]])
      y_new <- ifelse(total_count == 0, NULL, state_dict[[param]][["data"]][[total_count]])
      return(list(x = list(list(x_new)), y = list(list(y_new))))
    }
  )
  app$callback(
    output = list(id = sprintf("%s%s", param, suffix$ooc_n), property = "children"),
    params = list(
      input(id = "interval-component", property = "n_intervals"),
      state(id = "value-setter-store", property = "data")
    ),
    function(interval, stored_data) {
      total_count <- ifelse(interval < max_length, interval, max_length)
      ooc_n <- ifelse(total_count == 0, 0, stored_data[[param]][["ooc"]][[total_count]]*100) # ERROR?
      return(sprintf("%.2f%%", ooc_n))
    }
  )
  app$callback(
    output = list(id = sprintf("%s%s", param, suffix$ooc_g), property = "value"),
    params = list(
      input(id = "interval-component", property = "n_intervals"),
      state(id = "value-setter-store", property = "data")
    ),
    function(interval, stored_data) {
      total_count <- ifelse(interval < max_length, interval, max_length)
      ooc_n <- ifelse(total_count == 0, 0, stored_data[[param]][["ooc"]][[total_count]]*100)
      ooc_g <- ifelse(ooc_n == 0, 0.00001, ifelse(ooc_n <= 15, ooc_n, 15))
      return(ooc_g)
    }
  )
  app$callback(
    output = list(id = sprintf("%s%s", param, suffix$indicator), property = "color"),
    params = list(
      input(id = "interval-component", property = "n_intervals"),
      state(id = "value-setter-store", property = "data")
    ),
    function(interval, stored_data) {
      total_count <- ifelse(interval < max_length, interval, max_length)
      ooc_n <- ifelse(total_count == 0, 0, stored_data[[param]][["ooc"]][[total_count]]*100)
      color <- ifelse(ooc_n <= 5, theme$blue, "#FF0000")
      return(color)
    }
  )
}
for (param in params[-1]) {update_metric_summary(param)}

# callbacks to update pie chart
app$callback(
  output = list(id = "piechart", property = "figure"),
  params = list(
    input(id = "interval-component", property = "n_intervals"),
    state(id = "value-setter-store", property = "data"),
    state(id = "piechart", property = "figure")
  ),
  function(interval, stored_data, fig) {
    if (interval != 0) {
      total_count <- ifelse(interval < max_length, interval, max_length)
      values <- list()
      colors <- list()
      for (param in params[-1]) {
        ooc_param <- stored_data[[param]][["ooc"]][[total_count]]*100+1
        values <- append(values, ooc_param)
        colors <- append(colors, ifelse(ooc_param > 6, theme$red, theme$blue))
      }
      fig$data[[1]]$values <- values
      fig$data[[1]]$marker$colors <- colors
    }
    return(fig)
  }
)

# callbacks to update graph according to interval and selected parameter
app$callback(
  output = list(id = "control-chart-live", property = "figure"),
  params = list(
    input(id = "interval-component", property = "n_intervals"),
    input(id = sprintf("%s%s", params[-1][[1]], suffix$button_id), property = "n_clicks"),
    input(id = sprintf("%s%s", params[-1][[2]], suffix$button_id), property = "n_clicks"),
    input(id = sprintf("%s%s", params[-1][[3]], suffix$button_id), property = "n_clicks"),
    input(id = sprintf("%s%s", params[-1][[4]], suffix$button_id), property = "n_clicks"),
    input(id = sprintf("%s%s", params[-1][[5]], suffix$button_id), property = "n_clicks"),
    input(id = sprintf("%s%s", params[-1][[6]], suffix$button_id), property = "n_clicks"),
    input(id = sprintf("%s%s", params[-1][[7]], suffix$button_id), property = "n_clicks"),
    state(id = "value-setter-store", property = "data"),
    state(id = "control-chart-live", property = "figure")
  ),
  function(interval, n1, n2, n3, n4, n5, n6, n7, data, cur_fig) {
    # find which one was triggered
    ctx <- app$callback_context()
    if(!ctx$triggered$value) {
      id <- params[-1][[1]]
    } else {
      # get most recently triggered input id and property
      input <- unlist(strsplit(ctx$triggered$prop_id, "[.]"))
      if (input[[2]] == "n_clicks") {
        id <- gsub('.{7}$', '', input[[1]])
      } else if (input[[2]] == "n_intervals") {
        id <- cur_fig$data[[1]][["name"]]
      }
    }
    return(generate_graph(interval, data, id))
  }
)

########################################################################################################################
# RUN SERVER

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(port = 8080)
}
