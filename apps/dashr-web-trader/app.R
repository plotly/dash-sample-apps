library(plotly)
library(xts)
library(httr)
library(jsonlite)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(data.table)
library(fasttime)

appName <- Sys.getenv("DASH_APP_NAME")

if (!appName == ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

source("utils/helper-functions.R")

# load the R object with currency pairs data
load(file = "pairs/currencyPairs.RData")

currencyPairs <- list(EURUSD, USDCHF, USDJPY, GBPUSD)

# ensure timezone is correct in the data tables
myTZ <- Sys.timezone(location = TRUE)

for (df in currencyPairs){
  df[, Date := fastPOSIXct(Date, tz = myTZ)]
}

app <- Dash$new(name = "DashR Web Trader")

app$layout(
  htmlDiv(
    list(
      # interval component for live clock
      dccInterval(id = "interval", interval = 1 * 1000, n_intervals = 0),
      # interval component for ask bid updates
      dccInterval(id = "i_bis", interval = 1 * 2000, n_intervals = 0),
      # interval component for graph updates
      dccInterval(id = "i_tris", interval = 1 * 5000, n_intervals = 0),
      # interval for news
      dccInterval(id="i_news", interval = 24 * 60 * 60 * 1000, n_intervals=0),
      # left div
      htmlDiv(
        list(
          getLogo(),
          htmlDiv(
            list(
              htmlDiv(
                list(getHeader()),
                style = list(
                  backgroundColor = "#262a30",
                  paddingBottom = "15"
                ),
                id = "ask_bid_header",
                className = "row"
              ),
              htmlDiv(
                children = getFirstPairs(),
                style = list(
                  maxHeight = "45%",
                  backgroundColor = "#262a30",
                  color = "white",
                  paddingBottom = "15"
                ),
                id = "pairs"
              )
            ),
            id = "prices_div"
          ),
          # news div
          htmlDiv(
            list(
              htmlDiv(
                updateNews(),
                id = "news"
              )
            ),
            style = list(
              height = "33%",
              backgroundColor = "#262a30",
              color = "white",
              #fontSize = 12,
              padding = "1rem 1rem 0 1rem",
              marginTop = "0.5rem"
            )
          )
        ),
        className = "three columns",
        style = list(
          backgroundColor = "#262a30",
          padding = "1rem",
          height = "100%",
          overflow = "auto"
        )
      ),
      # center Div:
      htmlDiv(
        list(
          htmlDiv(
            getTopBar(),
            id = "top_bar",
            style = list(
              height = "15%",
              color = "white",
              backgroundColor = "#202228"
            )
          ),
          htmlDiv(
            lapply(currencyPairs, chartDiv),
            style = list(
              height = "65%",
              margin = "0px 5px"
            ),
            id = "charts",
            className = "row"
          ),
          bottomPanel()
        ),
        className = "nine columns",
        id = "rightpanel",
        style = list(
          backgroundColor = "#202228",
          height = "100vh",
          color = "white"
        )
      ),
      # hidden divs:
      htmlDiv(id = "charts_clicked", style = list(display = "none")),
      htmlDiv(ordersDiv()),
      htmlDiv(lapply(currencyPairs, modal)),
      htmlDiv(id = "orders", style = list(display = "none"))
    ),
    style = list(
      overflow = "auto",
      padding = "0",
      height = "110vh",
      backgroundColor = "#202228"
    )
  )
)

# Generate ask_bid_rows for each:
app$callback(
  output("EURUSDrow", "children"),
  list(
    input("i_bis", "n_intervals"),
    state("EURUSDindex", "children"),
    state("EURUSDbid", "children"),
    state("EURUSDask", "children")
  ),
  function(n_intervals, index, bid, ask){
    replaceRow(EURUSD, index, bid, ask)
  }
)
app$callback(
  output("GBPUSDrow", "children"),
  list(
    input("i_bis", "n_intervals"),
    state("GBPUSDindex", "children"),
    state("GBPUSDbid", "children"),
    state("GBPUSDask", "children")
  ),
  function(n_intervals, index, bid, ask){
    replaceRow(GBPUSD, index, bid, ask)
  }
)
app$callback(
  output("USDCHFrow", "children"),
  list(
    input("i_bis", "n_intervals"),
    state("USDCHFindex", "children"),
    state("USDCHFbid", "children"),
    state("USDCHFask", "children")
  ),
  function(n_intervals, index, bid, ask){
    replaceRow(USDCHF, index, bid, ask)
  }
)
app$callback(
  output("USDJPYrow", "children"),
  list(
    input("i_bis", "n_intervals"),
    state("USDJPYindex", "children"),
    state("USDJPYbid", "children"),
    state("USDJPYask", "children")),
  function(n_intervals, index, bid, ask){
    replaceRow(USDJPY, index, bid, ask)
  }
)

# generate_chart_button_callback
app$callback(
  output("charts_clicked", "children"),
  list(
    input("EURUSDButton_chart", "n_clicks"),
    input("USDCHFButton_chart", "n_clicks"),
    input("USDJPYButton_chart", "n_clicks"),
    input("GBPUSDButton_chart", "n_clicks"),
    state("charts_clicked", "children")
  ),
  function(eurusd, usdchf, usdjpy, gbpusd, clicked){
    pairs <- ""
    ifelse(
      (eurusd > 0 & nchar(pairs) > 1),
      pairs <- paste(pairs, "EURUSD", sep = ","),
      ifelse(
        eurusd > 0,
        pairs <- "EURUSD",
        pairs <- pairs)
    )
    ifelse(
      (usdchf > 0 & nchar(pairs) > 1),
      pairs <- paste(pairs, "USDCHF", sep = ","),
      ifelse(
        usdchf > 0,
        pairs <- "USDCHF",
        pairs <- pairs)
    )
    ifelse(
      (usdjpy > 0 & nchar(pairs) > 1),
      pairs <- paste(pairs, "USDJPY", sep = ","),
      ifelse(
        usdjpy > 0,
        pairs <- "USDJPY",
        pairs <- pairs)
    )
    ifelse(
      (gbpusd > 0 & nchar(pairs) > 1),
      pairs <- paste(pairs, "GBPUSD", sep = ","),
      ifelse(
        gbpusd > 0,
        pairs <- "GBPUSD",
        pairs <- pairs)
    )
    return(pairs)
  }
)

# generate_show_hide_graph_div_callback
app$callback(
  output("EURUSDgraph_div", "style"),
  list(input("charts_clicked", "children")),
  function(charts){
    charts <- unlist(strsplit(charts, ","))[1:4]
    if ("EURUSD" %in% charts){
      s <- list(
        position = "relative",
        float = "left",
        padding = "0px 10px",
        overflow = "hidden",
        marginBottom = "30px"
      )
      s["marginLeft"] <- "0px"
      ifelse(
        length(
          charts[!is.na(charts)]) == 4,
          s["height"] <- "50%",
          s["height"] <- "100%"
      )
    } else {
      s <- list(display = "none")
    }
  return(s)
  }
)
app$callback(
  output("USDCHFgraph_div", "style"),
  list(input("charts_clicked", "children")),
  function(charts){
    charts <- unlist(strsplit(charts, ","))[1:4]
    if ("USDCHF" %in% charts){
      s <- list(
        position = "relative",
        float = "left",
        padding = "0px 10px",
        overflow = "hidden",
        marginBottom = "30px"
      )
      s["marginLeft"] <- "0px"
      ifelse(
        length(charts[!is.na(charts)]) == 4,
        s["height"] <- "50%",
        s["height"] <- "100%"
      )
    }
    else {
      s <- list(display = "none")
    }
    return(s)
  }
)
app$callback(
  output("USDJPYgraph_div", "style"),
  list(input("charts_clicked", "children")),
  function(charts){
    charts <- unlist(strsplit(charts, ","))[1:4]
    if ("USDJPY" %in% charts){
      s <- list(
        position = "relative",
        float = "left",
        padding = "0px 10px",
        overflow = "hidden",
        marginBottom = "30px"
      )
      if (length(charts[!is.na(charts)]) == 4){
        s["marginLeft"] <- "0px"
      }
      ifelse(
        length(charts[!is.na(charts)]) == 4,
        s["height"] <- "50%",
        s["height"] <- "100%"
      )
    }
    else {
      s <- list(display = "none")
    }
    return(s)
  }
)
app$callback(
  output("GBPUSDgraph_div", "style"),
  list(input("charts_clicked", "children")),
  function(charts){
    charts <- unlist(strsplit(charts, ","))[1:4]
    if ("GBPUSD" %in% charts){
      s <- list(
        position = "relative",
        float = "left",
        padding = "0px 10px",
        overflow = "hidden",
        marginBottom = "30px"
      )
      s["marginLeft"] <- "0px"
      ifelse(
        length(charts[!is.na(charts)]) == 4,
        s["height"] <- "50%",
        s["height"] <- "100%"
      )
    }
    else {
      s <- list(display = "none")
    }
    return(s)
  }
)

# generate_size_graph_div
app$callback(
  output("EURUSDgraph_div", "className"),
  list(input("charts_clicked", "children")),
  function(charts){
    charts <- unlist(strsplit(charts, ","))[1:4]
    len_list <- length(charts[!is.na(charts)])
    if (len_list == 0 | !"EURUSD" %in% charts){
      return("")
    }
    width <- ifelse(
      len_list %% 2 == 0,
      "six columns",
      ifelse(
        len_list == 1,
        "twelve columns",
        "four columns"
      )
    )
    return(width)
  }
)
app$callback(
  output("USDCHFgraph_div", "className"),
  list(input("charts_clicked", "children")),
  function(charts){
    charts <- unlist(strsplit(charts, ","))[1:4]
    len_list <- length(charts[!is.na(charts)])
    if (len_list == 0 | !"USDCHF" %in% charts){
      return("")
    }
    width <- ifelse(len_list %% 2 == 0, "six columns",
      ifelse(
        len_list == 1,
        "twelve columns",
        "four columns"
      )
    )
    return(width)
  }
)
app$callback(
  output("USDJPYgraph_div", "className"),
  list(input("charts_clicked", "children")),
  function(charts){
    charts <- unlist(strsplit(charts, ","))[1:4]
    len_list <- length(charts[!is.na(charts)])
    if (len_list == 0 | !"USDJPY" %in% charts){
      return("")
    }
    width <- ifelse(
      len_list %% 2 == 0,
      "six columns",
      ifelse(
        len_list == 1,
        "twelve columns",
        "four columns"
      )
    )
    return(width)
  }
)
app$callback(
  output("GBPUSDgraph_div", "className"),
  list(input("charts_clicked", "children")),
  function(charts){
    charts <- unlist(strsplit(charts, ","))[1:4]
    len_list <- length(charts[!is.na(charts)])
    if (len_list == 0 | !"GBPUSD" %in% charts){
      return("")
    }
    width <- ifelse(
      len_list %% 2 == 0,
      "six columns",
      ifelse(
        len_list == 1,
        "twelve columns",
        "four columns"
        )
    )
    return(width)
  }
)

### generate_figure_callback
app$callback(
  output("EURUSDchart", "figure"),
  list(
    input("i_tris", "n_intervals"),
    input("EURUSDdropdown_period", "value"),
    input("EURUSDchart_type", "value"),
    input("EURUSDstudies", "value"),
    input("charts_clicked", "children"),
    state("EURUSDask", "children"),
    state("EURUSDbid", "children"),
    state("EURUSDchart", "figure"),
    state("EURUSDgraph_div", "style")
  ),
  function(n_i, p, t, s, pairs, a, b, old_fig, display){
    pairs <- unlist(strsplit(pairs, ","))
    if ("EURUSD" %in% pairs){
      f <- getFig(EURUSD, a, b, t, s, p)
      data <- f$x$data
      layout <- f$x$layout
      layout["uirevision"] <- display$position
      layout[["xaxis"]]["tickformat"] <- "%H:%M"
      layout[["margin"]] <- list(b = "20%", r = "0", l = "50", t = "0")
      layout[["xaxis"]]["showgrid"] <- FALSE
      layout["paper_bgcolor"] <- "#262a30"
      layout["plot_bgcolor"] <- "#262a30"
      layout["autosize"] <- TRUE
      return(list(data = data, layout = layout))
    }
    return(list(data = list(), layout = old_fig$layout))
  }
)
app$callback(
  output("USDCHFchart", "figure"),
  list(
    input("i_tris", "n_intervals"),
    input("USDCHFdropdown_period", "value"),
    input("USDCHFchart_type", "value"),
    input("USDCHFstudies", "value"),
    input("charts_clicked", "children"),
    state("USDCHFask", "children"),
    state("USDCHFbid", "children"),
    state("USDCHFchart", "figure"),
    state("USDCHFgraph_div", "style")
  ),
  function(n_i, p, t, s, pairs, a, b, old_fig, display){
    pairs <- unlist(strsplit(pairs, ","))
    if ("USDCHF" %in% pairs){
      f <- getFig(USDCHF, a, b, t, s, p)
      data <- f$x$data
      layout <- f$x$layout
      layout["uirevision"] <- display$position
      layout[["xaxis"]]["tickformat"] <- "%H:%M"
      layout[["margin"]] <- list(b = "20%", r = "0", l = "50", t = "0")
      layout[["xaxis"]]["showgrid"] <- FALSE
      layout["paper_bgcolor"] <- "#262a30"
      layout["plot_bgcolor"] <- "#262a30"
      layout["autosize"] <- TRUE
      return(list(data = data, layout = layout))
    }
    return(list(data = list(), layout = old_fig$layout))
  }
)
app$callback(
  output("USDJPYchart", "figure"),
  list(
    input("i_tris", "n_intervals"),
    input("USDJPYdropdown_period", "value"),
    input("USDJPYchart_type", "value"),
    input("USDJPYstudies", "value"),
    input("charts_clicked", "children"),
    state("USDJPYask", "children"),
    state("USDJPYbid", "children"),
    state("USDJPYchart", "figure"),
    state("USDJPYgraph_div", "style")
  ),
  function(n_i, p, t, s, pairs, a, b, old_fig, display){
    pairs <- unlist(strsplit(pairs, ","))
    if ("USDJPY" %in% pairs){
      f <- getFig(USDJPY, a, b, t, s, p)
      data <- f$x$data
      layout <- f$x$layout
      layout["uirevision"] <- display$position
      layout[["xaxis"]]["tickformat"] <- "%H:%M"
      layout[["margin"]] <- list(b = "20%", r = "0", l = "50", t = "0")
      layout[["xaxis"]]["showgrid"] <- FALSE
      layout["paper_bgcolor"] <- "#262a30"
      layout["plot_bgcolor"] <- "#262a30"
      layout["autosize"] <- TRUE
      return(list(data = data, layout = layout))
    }
    return(list(data = list(), layout = old_fig$layout))
  }
)
app$callback(
  output("GBPUSDchart", "figure"),
  list(
    input("i_tris", "n_intervals"),
    input("GBPUSDdropdown_period", "value"),
    input("GBPUSDchart_type", "value"),
    input("GBPUSDstudies", "value"),
    input("charts_clicked", "children"),
    state("GBPUSDask", "children"),
    state("GBPUSDbid", "children"),
    state("GBPUSDchart", "figure"),
    state("GBPUSDgraph_div", "style")
  ),
  function(n_i, p, t, s, pairs, a, b, old_fig, display){
    pairs <- unlist(strsplit(pairs, ","))
    if ("GBPUSD" %in% pairs){
      f <- getFig(GBPUSD, a, b, t, s, p)
      data <- f$x$data
      layout <- f$x$layout
      layout["uirevision"] <- display$position
      layout[["xaxis"]]["tickformat"] <- "%H:%M"
      layout[["margin"]] <- list(b = "20%", r = "0", l = "50", t = "0")
      layout[["xaxis"]]["showgrid"] <- FALSE
      layout["paper_bgcolor"] <- "#262a30"
      layout["plot_bgcolor"] <- "#262a30"
      layout["autosize"] <- TRUE
      return(list(data = data, layout = layout))
    }
    return(list(data = list(), layout = old_fig$layout))
  }
)

# generate_close_graph_callback
app$callback(
  output("EURUSDButton_chart", "n_clicks"),
  list(
    input("EURUSDclose", "n_clicks"),
    state("EURUSDButton_chart", "n_clicks")
  ),
  function(n, n2){
    if (n == 0){
      if (n2 == 1){
        return(1)
      }
      return(0)
    }
    return(0)
  }
)
app$callback(
  output("USDCHFButton_chart", "n_clicks"),
  list(
    input("USDCHFclose", "n_clicks"),
    state("USDCHFButton_chart", "n_clicks")
  ),
  function(n, n2){
    if (n == 0){
      if (n2 == 1){
        return(1)
      }
      return(0)
    }
    return(0)
  }
)
app$callback(
  output("USDJPYButton_chart", "n_clicks"),
  list(
    input("USDJPYclose", "n_clicks"),
    state("USDJPYButton_chart", "n_clicks")
  ),
  function(n, n2){
    if (n == 0){
      if (n2 == 1){
        return(1)
      }
      return(0)
    }
    return(0)
  }
)
app$callback(
  output("GBPUSDButton_chart", "n_clicks"),
  list(
    input("GBPUSDclose", "n_clicks"),
    state("GBPUSDButton_chart", "n_clicks")
  ),
  function(n, n2){
    if (n == 0){
      if (n2 == 1){
        return(1)
      }
      return(0)
    }
    return(0)
  }
)

# generate_active_menu_tab_callback
# updates hidden div that stores the last clicked menu tab
app$callback(
  output("EURUSDmenu_tab", "children"),
  list(
    input("EURUSDstyle_header", "n_clicks_timestamp"),
    input("EURUSDstudies_header", "n_clicks_timestamp")
  ),
  function(n_style, n_studies){
    ifelse(n_style >= n_studies, "Style", "Studies")
  }
)
app$callback(
  output("USDCHFmenu_tab", "children"),
  list(
    input("USDCHFstyle_header", "n_clicks_timestamp"),
    input("USDCHFstudies_header", "n_clicks_timestamp")
  ),
  function(n_style, n_studies){
    ifelse(n_style >= n_studies, "Style", "Studies")
  }
)
app$callback(
  output("USDJPYmenu_tab", "children"),
  list(
    input("USDJPYstyle_header", "n_clicks_timestamp"),
    input("USDJPYstudies_header", "n_clicks_timestamp")
  ),
  function(n_style, n_studies){
    ifelse(n_style >= n_studies, "Style", "Studies")
  }
)
app$callback(
  output("GBPUSDmenu_tab", "children"),
  list(
    input("GBPUSDstyle_header", "n_clicks_timestamp"),
    input("GBPUSDstudies_header", "n_clicks_timestamp")
  ),
  function(n_style, n_studies){
    ifelse(n_style >= n_studies, "Style", "Studies")
  }
)

# Generate Update Style Header Callback
app$callback(
  output("EURUSDstyle_header", "style"),
  list(
    input("EURUSDmenu_tab", "children"),
    state("EURUSDstyle_header", "style")
  ),
  function(current_tab, old_style){
    if (current_tab == "Style"){
      old_style["borderBottom"] <- "2px solid #45df7e"
    } else {
      old_style["borderBottom"] <- "2px solid rgba(68,149,209,.9)"
    }
    return(old_style)
  }
)
app$callback(
  output("USDCHFstyle_header", "style"),
  list(
    input("USDCHFmenu_tab", "children"),
    state("USDCHFstyle_header", "style")
  ),
  function(current_tab, old_style){
    if (current_tab == "Style"){
      old_style["borderBottom"] <- "2px solid #45df7e"
    } else {
      old_style["borderBottom"] <- "2px solid rgba(68,149,209,.9)"
    }
    return(old_style)
  }
)
app$callback(
  output("USDJPYstyle_header", "style"),
  list(
    input("USDJPYmenu_tab", "children"),
    state("USDJPYstyle_header", "style")
  ),
  function(current_tab, old_style){
    if (current_tab == "Style"){
      old_style["borderBottom"] <- "2px solid #45df7e"
    } else {
      old_style["borderBottom"] <- "2px solid rgba(68,149,209,.9)"
    }
    return(old_style)
  }
)
app$callback(
  output("GBPUSDstyle_header", "style"),
  list(
    input("GBPUSDmenu_tab", "children"),
    state("GBPUSDstyle_header", "style")
  ),
  function(current_tab, old_style){
    if (current_tab == "Style"){
      old_style["borderBottom"] <- "2px solid #45df7e"
    } else {
      old_style["borderBottom"] <- "2px solid rgba(68,149,209,.9)"
    }
    return(old_style)
  }
)

# generate_update_studies_header_callback
app$callback(
  output("EURUSDstudies_header", "style"),
  list(
    input("EURUSDmenu_tab", "children"),
    state("EURUSDstudies_header", "style")
  ),
  function(current_tab, old_studies){
    if (current_tab == "Studies"){
      old_studies["borderBottom"] <- "2px solid #45df7e"
    } else {
      old_studies["borderBottom"] <- "2px solid rgba(68,149,209,.9)"
    }
    return(old_studies)
  }
)
app$callback(
  output("USDCHFstudies_header", "style"),
  list(
    input("USDCHFmenu_tab", "children"),
    state("USDCHFstudies_header", "style")
  ),
  function(current_tab, old_studies){
    if (current_tab == "Studies"){
      old_studies["borderBottom"] <- "2px solid #45df7e"
    } else {
      old_studies["borderBottom"] <- "2px solid rgba(68,149,209,.9)"
    }
    return(old_studies)
  }
)
app$callback(
  output("USDJPYstudies_header", "style"),
  list(
    input("USDJPYmenu_tab", "children"),
    state("USDJPYstudies_header", "style")
  ),
  function(current_tab, old_studies){
    if (current_tab == "Studies"){
      old_studies["borderBottom"] <- "2px solid #45df7e"
    } else {
      old_studies["borderBottom"] <- "2px solid rgba(68,149,209,.9)"
    }
    return(old_studies)
  }
)
app$callback(
  output("GBPUSDstudies_header", "style"),
  list(
    input("GBPUSDmenu_tab", "children"),
    state("GBPUSDstudies_header", "style")
  ),
  function(current_tab, old_studies){
    if (current_tab == "Studies"){
      old_studies["borderBottom"] <- "2px solid #45df7e"
    } else {
      old_studies["borderBottom"] <- "2px solid rgba(68,149,209,.9)"
    }
    return(old_studies)
  }
)

#generate_studies_content_tab_callback
app$callback(
  output("EURUSDstudies_tab", "style"),
  list(input("EURUSDmenu_tab", "children")),
  function(current_tab){
    if (current_tab == "Studies"){
      return(list(display = "block", textAlign = "left", marginTop = "30"))
    } else {
      return(list(display = "none"))
    }
  }
)
app$callback(
  output("USDCHFstudies_tab", "style"),
  list(input("USDCHFmenu_tab", "children")),
  function(current_tab){
    if (current_tab == "Studies"){
      return(list(display = "block", textAlign = "left", marginTop = "30"))
    } else {
      return(list(display = "none"))
    }
  }
)
app$callback(
  output("USDJPYstudies_tab", "style"),
  list(input("USDJPYmenu_tab", "children")),
  function(current_tab){
    if (current_tab == "Studies"){
      return(list(display = "block", textAlign = "left", marginTop = "30"))
    } else {
      return(list(display = "none"))
    }
  }
)
app$callback(
  output("GBPUSDstudies_tab", "style"),
  list(input("GBPUSDmenu_tab", "children")),
  function(current_tab){
    if (current_tab == "Studies"){
      return(list(display = "block", textAlign = "left", marginTop = "30"))
    } else {
      return(list(display = "none"))
    }
  }
)

# generate_style_content_tab_callback
app$callback(
  output("EURUSDstyle_tab", "style"),
  list(input("EURUSDmenu_tab", "children")),
  function(current_tab){
    if (current_tab == "Style"){
      return(list(display = "block", textAlign = "left", marginTop = "30"))
    } else {
      return(list(display = "none"))
    }
  }
)
app$callback(
  output("USDCHFstyle_tab", "style"),
  list(input("USDCHFmenu_tab", "children")),
  function(current_tab){
    if (current_tab == "Style"){
      return(list(display = "block", textAlign = "left", marginTop = "30"))
    } else {
      return(list(display = "none"))
    }
  }
)
app$callback(
  output("USDJPYstyle_tab", "style"),
  list(input("USDJPYmenu_tab", "children")),
  function(current_tab){
    if (current_tab == "Style"){
      return(list(display = "block", textAlign = "left", marginTop = "30"))
    } else {
      return(list(display = "none"))
    }
  }
)
app$callback(
  output("GBPUSDstyle_tab", "style"),
  list(input("GBPUSDmenu_tab", "children")),
  function(current_tab){
    if (current_tab == "Style"){
      return(list(display = "block", textAlign = "left", marginTop = "30"))
    } else {
      return(list(display = "none"))
    }
  }
)

# generate_open_close_menu_callback
app$callback(
  output("EURUSDmenu", "className"),
  list(
    input("EURUSDmenu_button", "n_clicks"),
    input("EURUSDclose_menu", "n_clicks"),
    state("EURUSDmenu", "className")
  ),
  function(n, n2, currentClass){
    if (n == 0){
      return("not_visible")
    } else if (currentClass == "visible"){
      return("not_visible")
    } else {
      return("visible")
    }
  }
)
app$callback(
  output("USDCHFmenu", "className"),
  list(
    input("USDCHFmenu_button", "n_clicks"),
    input("USDCHFclose_menu", "n_clicks"),
    state("USDCHFmenu", "className")
  ),
  function(n, n2, currentClass){
    if (n == 0){
      return("not_visible")
    } else if (currentClass == "visible"){
      return("not_visible")
    } else {
      return("visible")
    }
  }
)
app$callback(
  output("USDJPYmenu", "className"),
  list(
    input("USDJPYmenu_button", "n_clicks"),
    input("USDJPYclose_menu", "n_clicks"),
    state("USDJPYmenu", "className")
  ),
  function(n, n2, currentClass){
    if (n == 0){
      return("not_visible")
    } else if (currentClass == "visible"){
      return("not_visible")
    } else {
      return("visible")
    }
  }
)
app$callback(
  output("GBPUSDmenu", "className"),
  list(
    input("GBPUSDmenu_button", "n_clicks"),
    input("GBPUSDclose_menu", "n_clicks"),
    state("GBPUSDmenu", "className")
  ),
  function(n, n2, currentClass){
    if (n == 0){
      return("not_visible")
    } else if (currentClass == "visible"){
      return("not_visible")
    } else {
      return("visible")
    }
  }
)

# generate_modal_open_callback
app$callback(
  output("EURUSDmodal", "style"),
  list(input("EURUSDBuy", "n_clicks")),
  function(n){
    if (n > 0){
      return(list(display = "block"))
    } else {
      return(list(display = "none"))
    }
  }
)
app$callback(
  output("USDCHFmodal", "style"),
  list(input("USDCHFBuy", "n_clicks")),
  function(n){
    if (n > 0){
      return(list(display = "block"))
    } else {
      return(list(display = "none"))
    }
  }
)
app$callback(
  output("USDJPYmodal", "style"),
  list(input("USDJPYBuy", "n_clicks")),
  function(n){
    if (n > 0){
      return(list(display = "block"))
    } else {
      return(list(display = "none"))
    }
  }
)
app$callback(
  output("GBPUSDmodal", "style"),
  list(input("GBPUSDBuy", "n_clicks")),
  function(n){
    if (n > 0){
      return(list(display = "block"))
    } else {
      return(list(display = "none"))
    }
  }
)

# generate_clean_sl_callback
app$callback(
  output("EURUSDSL", "value"),
  list(input("EURUSDBuy", "n_clicks")),
  function(n){
    return(0)
  }
)
app$callback(
  output("USDCHFSL", "value"),
  list(input("USDCHFBuy", "n_clicks")),
  function(n){
    return(0)
  }
)
app$callback(
  output("USDJPYSL", "value"),
  list(input("USDJPYBuy", "n_clicks")),
  function(n){
    return(0)
  }
)
app$callback(
  output("GBPUSDSL", "value"),
  list(input("GBPUSDBuy", "n_clicks")),
  function(n){
    return(0)
  }
)

# generate_clean_tp_callback
app$callback(
  output("EURUSDTP", "value"),
  list(input("EURUSDBuy", "n_clicks")),
  function(n){
    return(0)
  }
)
app$callback(
  output("USDCHFTP", "value"),
  list(input("USDCHFBuy", "n_clicks")),
  function(n){
    return(0)
  }
)
app$callback(
  output("USDJPYTP", "value"),
  list(input("USDJPYBuy", "n_clicks")),
  function(n){
    return(0)
  }
)
app$callback(
  output("GBPUSDTP", "value"),
  list(input("GBPUSDBuy", "n_clicks")),
  function(n){
    return(0)
  }
)

# generate_modal_close_callback
app$callback(
  output("EURUSDBuy", "n_clicks"),
  list(input("EURUSDcloseModal", "n_clicks"),
    input("EURUSDbutton_order", "n_clicks")),
  function(n, n2){
    return(0)
  }
)
app$callback(
  output("USDCHFBuy", "n_clicks"),
  list(input("USDCHFcloseModal", "n_clicks"),
    input("USDCHFbutton_order", "n_clicks")),
  function(n, n2){
    return(0)
  }
)
app$callback(
  output("USDJPYBuy", "n_clicks"),
  list(input("USDJPYcloseModal", "n_clicks"),
    input("USDJPYbutton_order", "n_clicks")),
  function(n, n2){
    return(0)
  }
)
app$callback(
  output("GBPUSDBuy", "n_clicks"),
  list(input("GBPUSDcloseModal", "n_clicks"),
    input("GBPUSDbutton_order", "n_clicks")),
  function(n, n2){
    return(0)
  }
)

# generate_modal_figure_callback
app$callback(
  output("EURUSDmodal_graph", "figure"),
  list(
    input("EURUSDindex", "children"),
    input("EURUSDBuy", "n_clicks"),
    state("EURUSDmodal_graph", "figure")
  ),
  function(index, n, old_fig){
    if (n == 1){
      return(getModalFig(EURUSD, index))
    }
    return(list(data = list(), layout = old_fig$layout))
  }
)
app$callback(
  output("USDCHFmodal_graph", "figure"),
  list(
    input("USDCHFindex", "children"),
    input("USDCHFBuy", "n_clicks"),
    state("USDCHFmodal_graph", "figure")
  ),
  function(index, n, old_fig){
    if (n == 1){
      return(getModalFig(USDCHF, index))
    }
    return(list(data = list(), layout = old_fig$layout))
  }
)
app$callback(
  output("USDJPYmodal_graph", "figure"),
  list(
    input("USDJPYindex", "children"),
    input("USDJPYBuy", "n_clicks"),
    state("USDJPYmodal_graph", "figure")
  ),
  function(index, n, old_fig){
    if (n == 1){
      return(getModalFig(USDJPY, index))
    }
    return(list(data = list(), layout = old_fig$layout))
  }
)
app$callback(
  output("GBPUSDmodal_graph", "figure"),
  list(
    input("GBPUSDindex", "children"),
    input("GBPUSDBuy", "n_clicks"),
    state("GBPUSDmodal_graph", "figure")
  ),
  function(index, n, old_fig){
    if (n == 1){
      return(getModalFig(GBPUSD, index))
    }
    return(list(data = list(), layout = old_fig$layout))
  }
)

# generate_order_button_callback
app$callback(
  output("EURUSDorders", "children"),
  list(
    input("EURUSDbutton_order", "n_clicks"),
    state("EURUSDvolume", "value"),
    state("EURUSDtrade_type", "value"),
    state("EURUSDSL", "value"),
    state("EURUSDTP", "value"),
    state("EURUSDorders", "children"),
    state("EURUSDask", "children"),
    state("EURUSDbid", "children")
  ),
  function(n, vol, type_order, sl, tp, pair_order, ask, bid){
    if (n > 0){
      j <- fromJSON(pair_order, flatten = TRUE)
      price <- ifelse(type_order == "sell", bid, ask)
      if (tp != 0){
        tp <- price + (tp * 0.001)
      }
      if (sl != 0){
        sl <- price - (tp * 0.001)
      }
      order <- list(
        id = paste0("EURUSD", as.character(n)),
        time = strftime(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
        type = type_order,
        volume = vol,
        symbol = "EURUSD",
        tp = tp,
        sl = sl,
        price = price,
        profit = 0.00,
        status = "open"
      )
      newlist <- list(j, data.frame(order))
      # remove empty element
      newlist <- newlist[lapply(newlist, length) > 0]
      json <- toJSON(rbindlist(newlist), auto_unbox = TRUE)
      return(json)
    }
    return(toJSON(list()))
  }
)
app$callback(
  output("USDCHForders", "children"),
  list(
    input("USDCHFbutton_order", "n_clicks"),
    state("USDCHFvolume", "value"),
    state("USDCHFtrade_type", "value"),
    state("USDCHFSL", "value"),
    state("USDCHFTP", "value"),
    state("USDCHForders", "children"),
    state("USDCHFask", "children"),
    state("USDCHFbid", "children")
  ),
  function(n, vol, type_order, sl, tp, pair_order, ask, bid){
    if (n > 0){
      j <- fromJSON(pair_order, flatten = TRUE)
      price <- ifelse(type_order == "sell", bid, ask)
      if (tp != 0){
        tp <- price + (tp * 0.001)
      }
      if (sl != 0){
        sl <- price - (tp * 0.001)
      }
      order <- list(
        id = paste0("USDCHF", as.character(n)),
        time = strftime(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
        type = type_order,
        volume = vol,
        symbol = "USDCHF",
        tp = tp,
        sl = sl,
        price = price,
        profit = 0.00,
        status = "open"
      )
      newlist <- list(j, data.frame(order))
      # remove empty element
      newlist <- newlist[lapply(newlist, length) > 0]
      json <- toJSON(rbindlist(newlist), auto_unbox = TRUE)
      return(json)
    }
    return(toJSON(list()))
  }
)
app$callback(
  output("USDJPYorders", "children"),
  list(
    input("USDJPYbutton_order", "n_clicks"),
    state("USDJPYvolume", "value"),
    state("USDJPYtrade_type", "value"),
    state("USDJPYSL", "value"),
    state("USDJPYTP", "value"),
    state("USDJPYorders", "children"),
    state("USDJPYask", "children"),
    state("USDJPYbid", "children")
  ),
  function(n, vol, type_order, sl, tp, pair_order, ask, bid){
    if (n > 0){
      j <- fromJSON(pair_order, flatten = TRUE)
      price <- ifelse(type_order == "sell", bid, ask)
      if (tp != 0){
        tp <- price + (tp * 0.001)
      }
      if (sl != 0){
        sl <- price - (tp * 0.001)
      }
      order <- list(
        id = paste0("USDJPY", as.character(n)),
        time = strftime(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
        type = type_order,
        volume = vol,
        symbol = "USDJPY",
        tp = tp,
        sl = sl,
        price = price,
        profit = 0.00,
        status = "open"
      )
      newlist <- list(j, data.frame(order))
      # remove empty element
      newlist <- newlist[lapply(newlist, length) > 0]
      json <- toJSON(rbindlist(newlist), auto_unbox = TRUE)
      return(json)
    }
    return(toJSON(list()))
  }
)
app$callback(
  output("GBPUSDorders", "children"),
  list(
    input("GBPUSDbutton_order", "n_clicks"),
    state("GBPUSDvolume", "value"),
    state("GBPUSDtrade_type", "value"),
    state("GBPUSDSL", "value"),
    state("GBPUSDTP", "value"),
    state("GBPUSDorders", "children"),
    state("GBPUSDask", "children"),
    state("GBPUSDbid", "children")
  ),
  function(n, vol, type_order, sl, tp, pair_order, ask, bid){
    if (n > 0){
      j <- fromJSON(pair_order, flatten = TRUE)
      price <- ifelse(type_order == "sell", bid, ask)
      if (tp != 0){
        tp <- price + (tp * 0.001)
      }
      if (sl != 0){
        sl <- price - (tp * 0.001)
      }
      order <- list(
        id = paste0("GBPUSD", as.character(n)),
        time = strftime(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
        type = type_order,
        volume = vol,
        symbol = "GBPUSD",
        tp = tp,
        sl = sl,
        price = price,
        profit = 0.00,
        status = "open"
      )
      newlist <- list(j, data.frame(order))
      # remove empty element
      newlist <- newlist[lapply(newlist, length) > 0]
      json <- toJSON(rbindlist(newlist), auto_unbox = TRUE)
      return(json)
    }
    return(toJSON(list()))
  }
)

# generate_update_orders_div_callback
app$callback(
  output("orders", "children"),
  list(
    input("EURUSDorders", "children"),
    input("USDCHForders", "children"),
    input("USDJPYorders", "children"),
    input("GBPUSDorders", "children"),
    input("EURUSDbid", "children"),
    input("USDCHFbid", "children"),
    input("USDJPYbid", "children"),
    input("GBPUSDbid", "children"),
    input("EURUSDask", "children"),
    input("USDCHFask", "children"),
    input("USDJPYask", "children"),
    input("GBPUSDask", "children"),
    input("closable_orders", "value"),
    state("orders", "children")
  ),
  function(...){
    arguments <- list(...)
    open_orders <- data.table()
    closed_orders <- data.table()
    current_orders <- unlist(arguments[length(arguments)])
    close_id <- unlist(arguments[length(arguments) - 1])
    arguments <- arguments[1:(length(arguments) - 2)]
    len_args <- length(arguments)
    bid_start <- floor(len_args / 3 + 1)
    current_bids <- arguments[bid_start:(bid_start + 3)]
    ask_start <- (floor(len_args / 3 ) * 2) + 1
    current_asks <- arguments[ask_start:(ask_start + 3)]
    arguments <- arguments[1:4]
    closed_id <- list()
    if (!is.null(current_orders)){
      co_df <- setDT(fromJSON(current_orders))
      closed_id <- c(closed_id, co_df[status == "closed", id])
      closed_orders <- rbind(closed_orders, co_df[status == "closed"])
    }
    for (list_order in arguments){
      if (list_order != "[]"){
        lo_df <- setDT(fromJSON(list_order))
        open_orders <- rbind(open_orders, lo_df[!id %in% closed_id])
      }
    }
    if (nrow(open_orders) == 0){
      return(current_orders)
    }
    orders <- updateOrders(open_orders, current_bids, current_asks, close_id)
    orders <- rbindlist(list(orders, closed_orders), fill = TRUE)
    return(toJSON(orders, auto_unbox = TRUE))
  }
)

# update_orders_table
app$callback(
  output("orders_table", "children"),
  list(
    input("orders", "children"),
    input("bottom_tab", "pathname")
  ),
  function(orders, url){
    if (url == "/"){
      url <- "open"
    } else{
      url <- "closed"
    }
    if (is.null(unlist(orders)) | orders == "[]"){
      return(ordersRows(NULL, url))
    }
    return(ordersRows(fromJSON(orders), url))
  }
)

# update link based on open / close dropdown
app$callback(
  output("bottom_tab", "pathname"),
  list(input("open_close_dropdown", "value")),
  function(val){
    val
  }
)

# update open / close dropdown values to reflect # of orders
app$callback(
  output("open_close_dropdown", "options"),
  list(input("orders", "children")),
  function(orders){
    n_open <- 0
    n_closed <- 0
    if (!is.null(unlist(orders))){
      orders <- fromJSON(orders)
      n_open <- sum(orders[, "status"] == "open")
      n_closed <- sum(orders[, "status"] == "closed")
    }
    list(
      list(
        label = sprintf("Open positions (%s)", n_open),
        value = "/"
      ),
      list(
        label = sprintf("Closed positions (%s)", n_closed),
        value = "/closed"
      )
    )
  }
)

# show_hide_close_orders
app$callback(
  output("close_orders_div", "style"),
  list(input("bottom_tab", "pathname")),
  function(url){
    if (url == "/"){
      style <- list(
        float = "right",
        right = "5",
        top = "3",
        width = "15%",
        textAlign = "center"
      )
    } else {
      style <- list(display = "none")
    }
    return(style)
  }
)

# update_close_dropdown
app$callback(
  output("closable_orders", "options"),
  list(input("orders", "children")),
  function(orders){
    options <- list()
    if (!is.null(unlist(orders))){
      orders <- fromJSON(orders)
      for (i in seq(1, nrow(orders))){
        if (orders[i, "status"] == "open"){
          options[[length(options) + 1]] <- list(
            label = orders[i, "id"],
            value = orders[i, "id"]
          )
        }
      }
    }
    return(options)
  }
)

# update_top_bar
app$callback(
  output("top_bar", "children"),
  list(input("orders", "children")),
  function(orders){
    if (is.null(unlist(orders))){
      return(getTopBar())
    }
    orders <- fromJSON(orders)
    open_pl <- 0
    balance <- 50000
    free_margin <- 50000
    margin <- 0
    for (i in 1:nrow(orders)){
      if (orders[i, "status"] == "open"){
        open_pl  <- open_pl + orders[i, "profit"]
        conversion_price <- ifelse(
          substr(orders[i, "symbol"], 1, 3) == "USD",
          1,
          orders[i, "price"]
        )
        margin <- margin +
          (orders[i, "volume"] * 100000) /
          (200 * conversion_price)
      } else {
        balance <- balance + orders[i, "profit"]
      }
    }
    equity <- balance - open_pl
    free_margin <- equity - margin
    margin_level <- ifelse(
      margin == 0,
      "%",
      paste0(as.character(round(equity / margin * 100, 2)), "%")
    )
    equity <- round(equity, 2)
    balance <- round(balance, 2)
    open_pl <- round(open_pl, 2)
    free_margin <- round(free_margin, 2)
    margin <- round(margin, 2)
    return(
      getTopBar(
        balance,
        equity,
        margin,
        free_margin,
        margin_level,
        open_pl
      )
    )
  }
)

# Callback for live clock
app$callback(
  output("live_clock", "children"),
  list(input("interval", "n_intervals")),
  function(n){
    strftime(Sys.time(), format = "%H:%M:%S")
  }
)

app$callback(
  # update news every 24 hours
  output("news", "children"),
  list(input("i_news", "n_intervals")),
  function(n){
    updateNews()
  }
)

if (appName == ""){
  app$run_server(debug = TRUE)
} else {app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8080))}
