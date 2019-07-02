# This file contains functions used to generate most of
# the divs and plots used in app.R
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(data.table)

# Generate plotly logo
getLogo <- function(){
  htmlImg(
    src = "https://user-images.githubusercontent.com/37411533/56820500-10b13c80-681a-11e9-9519-4236cee3551f.png",
    className = "sept columns"
  )
}

# Generate the empty news table (content filled by updateNews function)
generateNewsTable <- function(dataframe, max_rows = 10){
  n <- min(max_rows, nrow(dataframe))
  rows <- lapply(seq(1, n), function(i){
    htmlTr(
      htmlTd(
        htmlA(
          dataframe[i, "title"],
          href = as.character(dataframe[i, "url"]),
          style = list(
            textDecoration = "none"
          )
        ),
        className = "newsRow",
        style = list(
          fontSize = 15,
          padding = "0px 10px 0px 10px"
        ),
      )
    )
  })
  htmlDiv(
    list(
      htmlP(
        list(
          "Headlines",
          htmlBr(),
          htmlSpan(
            id = "news_update",
            children = paste(
              "Last updated :",
              strftime(Sys.time(), format = "%Y-%m-%d %H:%M:%S")
            ),
            style = list(
              fontSize = 12,
              marginTop = 4
            )
          )
        ),
        style = list(marginTop = "4", color = "#b4b4b4")
      ),
      htmlBr(),
      htmlDiv(
        list(
          htmlTable(
            rows
          )
        ),
        style = list(
          height = "400px",
          overflowY = "scroll",
          border = "1px solid #56585d"
        )
      )
    ),
    style = list(height = "100%")
  )
}

# Get the news from newsapi.org
updateNews <- function(){
  apiURL <- "https://newsapi.org/v2/top-headlines?sources=bbc-news&apiKey=da8e2e705b914f9f86ed2e9692e66012" 
  r <- GET(apiURL)
  jsonData <- content(r, "parsed")$articles
  # include tryCatch for when the api breaks
  out <- tryCatch({
    d <- data.frame(do.call("rbind", jsonData))
    d <- d[, c("title", "url")]
  },
    error = function(x){
      d <- data.frame(title = "news api not available",
                      url = apiURL)
    }
  )
  generateNewsTable(out)
}

# Header for left div (table header with live clock)
getHeader <- function(t = Sys.time()){
  htmlDiv(
    list(
      htmlP(
        strftime(Sys.time(), format = "%H:%M:%S"),
        id = "live_clock",
        className = "four columns",
        style = list(color = "#b4b4b4", textAlign = "center")
      ),
      htmlP(
        "Bid",
        className = "four columns",
        style = list(color = "white", textAlign = "center")
      ),
      htmlP(
        "Ask",
        className = "four columns",
        style = list(color = "white", textAlign = "center")
      )
    )
  )
}

# color of Bid & Ask rates (depending on previous bids and asks)
getColor <- function(a, b){
  ifelse(
    a == b,
    "white",
    ifelse(
      a > b,
      "#45df7e",
      "#da5657"
    )
  )
}

# Row or previous 10 rows of currency pair dataset
getAskBid <- function(currency_pair, index, modal = FALSE){
  if (modal == FALSE){
    currency_pair[index, ]
  } else {
    currency_pair[c( (index - 10):index), ]
  }
}

# returns dataset row with nearest datetime to current time
firstAskBid <- function(pair){
  t <- as.POSIXct(
    paste("2016-01-05", strftime(Sys.time(), format = "%H:%M:%OS3")),
    tz = myTZ
  )
  dfRow <- pair[, .SD[which.min(abs(difftime(pair$Date, t)))]]
  intIndex <- which.min(abs(difftime(pair$Date, t)))
  list(dfRow, intIndex)
}

# creates a Bid & Ask row for a currency pair + buttons
getRow <- function(data){
  currentRow <- data[[1]]
  index <- data[[2]]
  htmlDetails(
    list(
      htmlSummary(
        list(
          htmlDiv(
            list(
              htmlP(
                children = currentRow$Symbol,
                id = currentRow$Symbol,
                className = "four columns",
                style = list(textAlign = "center")
              ),
              htmlP(
                children = currentRow$Bid,
                id = paste0(currentRow$Symbol, "bid"),
                className = "four columns",
                style = list(textAlign = "center", color = "white")
              ),
              htmlP(
                children = currentRow$Ask,
                id = paste0(currentRow$Symbol, "ask"),
                className = "four columns",
                style = list(textAlign = "center", color = "white")
              ),
              htmlDiv(children = index,
                id = paste0(currentRow$Symbol, "index"),
                style = list(display = "none")
              )
            ),
            id = paste0(currentRow$Symbol, "row"),
            className = "row eleven columns",
            style = list(height = "40", float = "right")
          )
        ),
        className = "row askBidRow",
        style = list(paddingLeft = "10")
      ),
      htmlDiv(
        list(
          htmlButton(
            children = "Buy/Sell",
            id = paste0(currentRow$Symbol, "Buy"),
            n_clicks = 0,
            style = list(
              margin = "0px 7px 7px 10px",
              textAlign = "center"
            ),
          ),
          htmlButton(
            children = "Chart",
            id = paste0(currentRow$Symbol, "Button_chart"),
            n_clicks = ifelse(
              currentRow$Symbol %in% c("EURUSD", "USDCHF"), 1, 0
            ),
            style = list(
              margin = "0 3rem 1rem 3rem",
              textAlign = "center"
            )
          )
        ),
        style = list(textAlign = "right")
      )
    ),
    id = paste0(currentRow$Symbol, "row_div"),
    n_clicks = 0,
    style = list(textAlign = "center", paddingTop = "4")
  )
}

# Return an updated ask and bid row for given currency pair
replaceRow <- function(currency_pair, index, bid, ask){
  index <- index + 1
  newRow <- ifelse(
    index != nrow(currency_pair),
    list(getAskBid(currency_pair, index)),
    firstAskBid(currency_pair)
  )
  list(
    htmlP(
      children = newRow[[1]]$Symbol,
      id = newRow[[1]]$Symbol,
      className = "four columns",
      style = list(textAlign = "center")
    ),
    htmlP(
      children = newRow[[1]]$Bid,
      id = paste0(newRow[[1]]$Symbol, "bid"),
      className = "four columns",
      style = list(
        textAlign = "center",
        color = getColor(newRow[[1]]$Bid, bid)
      )
    ),
    htmlP(
      children = newRow[[1]]$Ask,
      id = paste0(newRow[[1]]$Symbol, "ask"),
      className = "four columns",
      style = list(
        textAlign = "center",
        color = getColor(newRow[[1]]$Ask, ask)
      )
    ),
    htmlDiv(
      children = index,
      id = paste0(newRow[[1]]$Symbol, "index"),
      style = list(display = "none")
    )
  )
}

# Get the initial ask and bids
getFirstPairs <- function(){
  lapply(currencyPairs,
         function(x){
           getRow(firstAskBid(x))
         }
  )
}

# return one cell of top bar
getTopBarCell <- function(up, down, color = "white"){
  htmlDiv(
    list(
      htmlP(
        children = up,
        style = list(
          marginBottom = "1px",
          color = "#b4b4b4",
          fontSize = 14
        )
      ),
      htmlP(
        children = down,
        id = up,
        style = list(marginBottom = "3", fontSize = 18)
      )
    ),
    className = "one columns",
    style = list(
      textAlign = "center",
      width = "16%",
      color = color,
      paddingBottom = "15px"
    )
  )
}

# return top bar with updated values
getTopBar <- function(balance=50000.00,
                      equity=50000,
                      margin=0,
                      fm=50000,
                      m_level="%",
                      open_pl=0){
  colorOpenPL <- getColor(open_pl, 0)
  htmlDiv(
    list(
      getTopBarCell("Balance", balance),
      getTopBarCell("Equity", equity),
      getTopBarCell("Margin", margin),
      getTopBarCell("Free Margin", fm),
      getTopBarCell("Margin Level", m_level),
      getTopBarCell("Open P/L", open_pl, color = colorOpenPL)
    ),
    style = list(padding = "55px 0px 15px 0px")
  )
}

# Get OHLC data for given pair
getOHLCData <- function(currency_pair,
                        period = c("5Min", "15Min", "30Min")){
  t <- as.POSIXct(
    strftime(Sys.time(), format = "2016-01-05 %H:%M:%S")
  )
  d <- currency_pair[Date < t, c("Date", "Bid")]
  if (period == "5Min"){
    as.data.table(to.minutes5(d))
  } else if (period == "15Min"){
    as.data.table(to.minutes15(d))
  } else {
    as.data.table(to.minutes30(d))
  }
}

## Traces For Studies
# moving average
movingAverageTrace <- function(df, fig){
  df <- df[, c("closeRA") := lapply(.SD, rollmean, k = 5, fill = NA),
             .SDcols = "d.Close"]
  fig %>%
    add_trace(inherit = FALSE,
              data = df, x = ~index, y = ~closeRA,
              type = "scatter", mode = "lines",
              showlegend = FALSE,
              name = "MA",
              fill = "none"
    ) %>%
    layout(
      yaxis = list(title = "")
    )
}

# exponential moving average (copied from python app)
eMovingAverageTrace <- function(df, fig){
  df <- df[, c("closeERA") := lapply(.SD, rollmean, k = 20, fill = NA),
             .SDcols = "d.Close"]
  fig %>%
    add_trace(inherit = FALSE, data = df, x = ~index, y = ~closeERA,
              type = "scatter", mode = "lines",
              showlegend = FALSE,
              name = "EMA",
              fill = "none"
    ) %>%
    layout(
      yaxis = list(title = "")
    )
}

# bollinger Bands
bollingerBands <- function(df, fig, window_size = 10, num_of_std = 5){
  df <- df[, c("RollingMean") := lapply(
                  .SD, rollmean, k =
                    window_size, fill = NA
             ),
             .SDcols = "d.Close"]
  df <- df[, c("RollingSTD") := lapply(
                  .SD, function(x){
                      rollapply(x, width = window_size, FUN = sd, fill = NA)
                  }
             ),
             .SDcols = "d.Close"]
  df$upperBand <- df[, "RollingMean"] + (df[, "RollingSTD"] * num_of_std)
  df$lowerBand <- df[, "RollingMean"] - (df[, "RollingSTD"] * num_of_std)
  fig %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~upperBand,
      type = "scatter", mode = "lines",
      name = "BB_upper",
      showlegend = FALSE,
      fill = "none"
    ) %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~RollingMean,
      type = "scatter", mode = "lines",
      name = "BB_mean",
      showlegend = FALSE,
      fill = "none"
    ) %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~lowerBand,
      type = "scatter", mode = "lines",
      name = "BB_lower",
      showlegend = FALSE,
      fill = "none"
    ) %>%
    layout(
      yaxis = list(title = "")
    )
}

# Accumulation Distribution
accumulationTrace <- function(df){
  df$volume <- ( (df$d.Close - df$d.Low) -
                (df$d.High - df$d.Close)) /
                (df$d.High - df$d.Low)
  plot_ly(df, x = ~index, y = ~volume,
    mode = "lines", type = "scatter",
    showlegend = FALSE,
    name = "Accumulation"
  )
}

# Commodity Channel Index
cciTrace <- function(df){
  df[, c("TP") := (df[, "d.High"] + df[, "d.Close"] + df[, "d.Low" ]) / 3]
  df[, c("TP_RM") := lapply(
         .SD, rollmean, k = 10, fill = NA, align = "right"
         ),
      .SDcols = "TP"]
  df[, c("TP_RSD") := lapply(
         .SD, function(x){
           rollapply(x, width = 10, FUN = sd, fill = NA, align = "right")
         }),
       .SDcols = "TP"]
  df$CCI <- (df$TP - df$TP_RM) / (0.015 * df$TP_RSD)
  plot_ly(
    df, x = ~index, y = ~CCI,
    mode = "lines", type = "scatter",
    showlegend = FALSE,
    name = "CCI"
  )
}

# Price Rate of Change
rocTrace <- function(df, n_days = 5){
  df[, c("N") := d.Close - shift(d.Close, 5)]
  df[, c("D") := shift(d.Close, 5)]
  df[, c("ROC") := N / D]
  plot_ly(
    df, x = ~index, y = ~ROC,
    mode = "lines", type = "scatter",
    showlegend = FALSE,
    name = "ROC"
  )
}

# Stochastic Oscillator:
stocTrace <- function(df){
  df[, c("SOk") := (d.Close - d.Low) / (d.High - d.Low)]
  plot_ly(
    df, x = ~index, y = ~SOk,
    mode = "lines", type = "scatter",
    showlegend = FALSE,
    name = "Stochastic Oscillator"
  )
}

# momentum Trace
momTrace <- function(df, n = 5){
  df[, c("M") := d.Close - shift(d.Close, 5)]
  plot_ly(
    df, x = ~index, y = ~M,
    mode = "lines", type = "scatter",
    showlegend = FALSE,
    name = "MOM"
  )
}

# Pivot Points
ppTrace <- function(df, fig){
  df[, c("PP") := ( (df[, "d.High"] + df[, "d.Close"] + df[, "d.Low" ]) / 3)]
  df[, c("R1") := 2 * PP - d.Low]
  df[, c("S1") := 2 * PP - d.High]
  df[, c("R2") := PP + d.High - d.Low]
  df[, c("S2") := PP - d.High + d.Low]
  df[, c("R3") := d.High + 2 * (PP - d.Low)]
  df[, c("S3") := d.Low - 2 * (d.High - PP)]
  fig %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~PP,
      type = "scatter", mode = "lines",
      showlegend = FALSE,
      name = "PP",
      fill = "none"
    ) %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~R1,
      type = "scatter", mode = "lines",
      showlegend = FALSE,
      name = "R1",
      fill = "none"
    ) %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~R2,
      type = "scatter", mode = "lines",
      showlegend = FALSE,
      name = "R2",
      fill = "none"
    ) %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~R3,
      type = "scatter", mode = "lines",
      showlegend = FALSE,
      name = "R3",
      fill = "none"
    ) %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~S1,
      type = "scatter", mode = "lines",
      showlegend = FALSE,
      name = "S1",
      fill = "none"
    ) %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~S2,
      type = "scatter", mode = "lines",
      showlegend = FALSE,
      name = "S2",
      fill = "none"
    ) %>%
    add_trace(
      inherit = FALSE,
      data = df, x = ~index, y = ~S3,
      type = "scatter", mode = "lines",
      showlegend = FALSE,
      name = "S3",
      fill = "none"
    ) %>%
    layout(
      yaxis = list(title = "")
    )
}

# Main chart traces
lineTrace <- function(df){
  plot_ly(data = df, x = ~index, y = ~d.Close,
    type = "scatter",
    mode = "lines",
    showlegend = FALSE,
    name = "line"
  )
}

areaTrace <- function(df){
  plot_ly(data = df, x = ~index, y = ~d.Close,
    type = "scatter",
    mode = "lines",
    fill = "toself",
    showlegend = FALSE,
    name = "area"
  )
}

barTrace <- function(df){
  plot_ly(data = df, type = "ohlc", x = ~index,
    open = ~d.Open, high = ~d.High, low = ~d.Low, close = ~d.Close,
    increasing = list(line = list(color = "#888888")),
    decreasing = list(line = list(color = "#888888")),
    showlegend = FALSE,
    name = "bar"
  ) %>%
    layout(xaxis = list(rangeslider = list(visible = FALSE)))
}

coloredBarTrace <- function(df){
  plot_ly(data = df, type = "ohlc", x = ~index,
    open = ~d.Open, high = ~d.High, low = ~d.Low, close = ~d.Close,
    showlegend = FALSE,
    name = "colored bar",
    tickwidth = 0.5
  ) %>%
    layout(xaxis = list(rangeslider = list(visible = FALSE)))
}

candlestickTrace <- function(df){
  plot_ly(data = df, type = "candlestick", x = ~index,
    open = ~d.Open, high = ~d.High, low = ~d.Low, close = ~d.Close,
    increasing = list(line = list(color = "#00ff00")),
    decreasing = list(line = list(color = "white")),
    showlegend = FALSE,
    name = "candlestick"
  ) %>%
    layout(xaxis = list(rangeslider = list(visible = FALSE)))
}

askModalTrace <- function(currency_pair, index){
  d <- getAskBid(currency_pair, index, modal = TRUE)
  plot_ly(data = d, x = ~Date, y = ~Ask,
    type = "scatter", mode = "lines",
    showlegend = FALSE
  )
}

bidModalTrace <- function(currency_pair, index){
  d <- getAskBid(currency_pair, index, modal = TRUE)
  plot_ly(
    data = d, x = ~Date, y = ~Bid,
    type = "scatter", mode = "lines",
    showlegend = FALSE
  )
}

getModalFig <- function(currency_pair, index){
  ask <- askModalTrace(currency_pair, index)
  bid <- bidModalTrace(currency_pair, index)
  subplot(ask, bid, nrows = 2, shareX = TRUE, shareY = FALSE) %>%
    layout(
      margin = list(b = 0, r = 5, l = 10, t = 5),
      xaxis = list(showticklabels = FALSE, title = ""),
      paper_bgcolor = "#262a30",
      plot_bgcolor = "#262a30",
      autosize = TRUE
    )
}

getFig <- function(currency_pair, ask, bid,
                   type_trace, studies,
                   period = c("5Min", "15Min", "30Min")){
  d <- getOHLCData(currency_pair, period)
  subplotTraces <- list(
    "accumulationTrace",
    "cciTrace",
    "rocTrace",
    "stocTrace",
    "momTrace"
  )
  selectedSubplotStudies <- list()
  #selectedFirstRowStudies <- list()
  type_trace_fun <- match.fun(type_trace)
  fig <- type_trace_fun(d)
  totalRows <- 1
  if (is.list(studies) & length(studies) > 0){
    for (study in studies){
      if (study %in% subplotTraces){
        totalRows <- totalRows + 1
        selectedSubplotStudies <- c(selectedSubplotStudies, study)
      } else {
        s <- match.fun(study)
        fig <- add_trace(
          s(d, fig),
          type = "scatter",
          mode = "lines",
          inherit = FALSE
        )
      }
    }
  }
  plotList <- vector("list", totalRows)
  plotList[[1]] <- fig
  currentRow <- 1
  for (study in selectedSubplotStudies){
    currentRow <- currentRow + 1
    s <- match.fun(study)
    plotList[[currentRow]]  <- s(d)
  }
  subplot(plotList, nrows = totalRows, shareX = TRUE, shareY = TRUE) %>%
    layout(
      margin = list(b = "20%", r = "0", l = "50", t = "0"),
      xaxis = list(
        tickformat = "%H:%M",
        rangeslider = list(visible = FALSE)
      ),
      paper_bgcolor = "#262a30",
      plot_bgcolor = "#262a30"
    )
}

chartDiv <- function(pair){
  pairname <- pair$Symbol[[1]]
  htmlDiv(
    list(
      # Menu Div
      htmlDiv(
        list(
          # current menu tab:
          htmlDiv(
            "Studies",
            id = paste0(pairname, "menu_tab"),
            style = list(display = "none")
          ),
          htmlSpan(
            "x",
            id = paste0(pairname, "close_menu"),
            n_clicks = 0,
            style = list(
              fontsize = "16",
              float = "right",
              paddingRight = "5",
              verticalAlign = "textTop",
              cursor = "pointer"
            )
          ),
          htmlSpan(
            "Style",
            id = paste0(pairname, "style_header"),
            n_clicks_timestamp = 2,
            style = list(
              top = "0",
              float = "left",
              marginLeft = "5",
              marginRight = "5",
              textDecoration = "none",
              cursor = "pointer"
            )
          ),
          htmlSpan(
            "Studies",
            id = paste0(pairname, "studies_header"),
            n_clicks_timestamp = 1,
            style = list(
              float = "left",
              textDecoration = "none",
              cursor = "pointer"
            )
          ),
          htmlBr(),
          htmlDiv(
            htmlDiv(
              dccChecklist(
                id = paste0(pairname, "studies"),
                options = list(
                  list(
                    label = "Accumulation/D", value = "accumulationTrace"
                  ),
                  list(
                    label = "Bollinger bands", value = "bollingerBands"
                  ),
                  list(
                    label = "MA", value = "movingAverageTrace"
                  ),
                  list(
                    label = "EMA", value = "eMovingAverageTrace"
                  ),
                  list(
                    label = "CCI", value = "cciTrace"
                  ),
                  list(
                    label = "ROC", value = "rocTrace"
                  ),
                  list(
                    label = "Pivot points", value = "ppTrace"
                  ),
                  list(
                    label = "Stochastic oscillator", value = "stocTrace"
                  ),
                  list(
                    label = "Momentum indicator", value = "momTrace"
                  )
                ),
                value = list()
              ),
              style = list(marginTop = "30", textAlign = "left")
            ),
            id = paste0(pairname, "studies_tab"),
            style = list(display = "none")
          ),
          htmlDiv(
            dccRadioItems(
              id = paste0(pairname, "chart_type"),
              options = list(
                list(label = "candlestick", value = "candlestickTrace"),
                list(label = "line", value = "lineTrace"),
                list(label = "mountain", value = "areaTrace"),
                list(label = "bar", value = "barTrace"),
                list(label = "colored bar", value = "coloredBarTrace")
              ),
              value = "coloredBarTrace"
            ),
            id = paste0(pairname, "style_tab"),
            style = list(marginTop = "30", textAlign = "left")
          )
        ),
        id = paste0(pairname, "menu"),
        className = "not_visible",
        style = list(
          overflow = "auto",
          backgroundColor = "#32363f",
          zIndex = "20",
          width = "45%",
          height = "60%",
          position = "absolute",
          left = "10px",
          top = "50px",
          border = "1px solid #808080",
          borderRadius = "4px"
        )
      ),
      # Top Bar:
      htmlDiv(
        list(
          htmlDiv(
            list(
              htmlSpan(
                pairname,
                style = list(
                  margin = "0 1rem 0 2rem",
                  paddingLeft = "25",
                  marginRight = "15",
                  fontSize = 18
                )
              ),
              htmlSpan(
                "☰",
                id = paste0(pairname, "menu_button"),
                n_clicks = 0,
                style = list(
                  margin = "0 1rem",
                  cursor = "pointer",
                  fontSize = 18
                )
              )
            )
          ),
          htmlDiv(
            list(
              htmlDiv(
                dccDropdown(
                  id = paste0(pairname, "dropdown_period"),
                  options = list(
                    list(label = "5m", value = "5Min"),
                    list(label = "15m", value = "15Min"),
                    list(label = "30m", value = "30Min")
                  ),
                  value = "15Min",
                  clearable = FALSE,
                  searchable = FALSE
                ),
                style = list(
                  width = "100px",
                  paddingRight = "30"
                )
              ),
              htmlSpan(
                "x",
                id = paste0(pairname, "close"),
                n_clicks = 0,
                style = list(
                  margin = "0 2rem",
                  color = "#b4b4b4",
                  cursor = "pointer",
                  fontSize = 18,
                  paddingRight = "25"
                )
              )
            ),
            style = list(
              display = "flex",
              alignItems = "center"
            )
          )
        ),
        id = paste0(pairname, "toprow"),
        style = list(
          display = "flex",
          justifyContent = "space-between",
          alignItems = "center",
          height = "50px",
          color = "#b4b4b4",
          backgroundColor = "#32363f"
        )
      ),
      # GraphDiv:
      htmlDiv(
        dccGraph(
          id = paste0(pairname, "chart"),
          config = list(displayModeBar = FALSE, scrollZoom = TRUE),
          style = list(width = "100%", height = "100%")
        ),
        id = paste0(pairname, "graph"),
        style = list(width = "100%", height = "100%")
      )
    ),
    id = paste0(pairname, "graph_div"),
    className = "",
    style = list(display = "none")
  )
}

bottomPanel <- function(){
  htmlDiv(
    list(
      htmlDiv(
        list(
          dccLocation(id = "bottom_tab", refresh = FALSE),
          htmlDiv(
            dccDropdown(
              id = "open_close_dropdown",
              style = list(
                paddingLeft = "2rem",
                border = "none"
                ),
              options = list(
                list(label = "Open Positions", value = "/"),
                list(label = "Closed Positions", value = "/closed")
              ),
              value = "/",
              clearable = FALSE
            ),
            style = list(width = "21rem")
          ),
          htmlDiv(
            dccDropdown(
              id = "closable_orders",
              placeholder = "Close order",
              style = list(paddingLeft = "18px")
            ),
            id = "close_orders_div",
          )
        ),
        style = list(
          backgroundColor = "#2d2f36",
          display = "flex",
          height = "50px",
          justifyContent = "space-between",
          alignItems = "center",
          padding = "0px 25px"
        )
      ),
      htmlDiv(
        children = list(
          htmlTable(id = "orders_table")
        ),
        className = "row",
        style = list(
          padding = "3",
          textalign = "center"
        ),
        id = "bottom_content"
      )
    ),
    id = "bottom_panel",
    style = list(
      overflowY = "auto",
      marginLeft = "15px",
      marginRight = "auto",
      height = "20%",
      width = "calc(100% - 35px)",
      backgroundColor = "#262a30"
    )
  )
}

modal <- function(pair){
  pairname <- pair$Symbol[1]
  htmlDiv(
    htmlDiv(
      list(
        htmlDiv(
          list(
            htmlSpan(
              "x",
              id = paste0(pairname, "closeModal"),
              style = list(
                float = "right",
                cursor = "pointer",
                color = "#b6b6b6",
                fontSize = 18
              )
            ),
            htmlSpan(
              pairname,
              id = paste0("modal", pairname)
            ),
            htmlDiv(
              list(
                htmlDiv(
                  list(
                    dccGraph(
                      id = paste0(pairname, "modal_graph"),
                      config = list(displayModeBar = FALSE)
                    )
                  ),
                  className = "seven columns"
                ),
                htmlDiv(
                  list(
                    htmlDiv(
                      list(
                        htmlP(
                          "Volume",
                          style = list(marginBottom = "0")
                        ),
                        dccInput(
                          id = paste0(pairname, "volume"),
                          type = "number",
                          value = 0.1,
                          min = 0,
                          step = 0.1
                        )
                      ),
                      style = list(marginBottom = "5")
                    ),
                    htmlDiv(
                      list(
                        htmlP(
                          "Type", style = list(marginBottom = "0")
                        ),
                        dccRadioItems(
                          id = paste0(pairname, "trade_type"),
                          options = list(
                            list(label = "Buy", value = "buy"),
                            list(label = "Sell", value = "sell")
                          ),
                          value = "buy",
                          labelStyle = list(display = "inline-block")
                        )
                      ),
                      style = list(marginBottom = "5"),
                    ),
                    htmlDiv(
                      list(
                        htmlP(
                          "SL TPS",
                          style = list(marginBottom = "0")
                        ),
                        dccInput(
                          id = paste0(pairname, "SL"),
                          type = "number",
                          min = 0,
                          step = 1
                        )
                      ),
                      style = list(marginBottom = "5")
                    ),
                    htmlDiv(
                      list(
                        htmlP(
                          "TP TPS",
                          style = list(marginBottom = "0")
                        ),
                        dccInput(
                          id = paste0(pairname, "TP"),
                          type = "number",
                          min = 0,
                          step = 1
                        )
                      ),
                      style = list(marginBottom = "5")
                    )
                  ),
                  className = "four columns"
                )
              ),
              style = list(
                display = "flex",
                justifyContent = "space-between"
              )
            ),
            htmlDiv(
              htmlButton(
                "Order",
                id = paste0(pairname, "button_order"),
                n_clicks = 0
              ),
              style = list(textAlign = "center", marginTop = "12")
            )
          ),
          className = "modal-content"
        )
      ),
      className = "modal"
    ),
    id = paste0(pairname, "modal"),
    style = list(display = "none", color = "white")
  )
}

ordersDiv <- function(){
  list(
    htmlDiv(id = "EURUSDorders", style = list(display = "none")),
    htmlDiv(id = "USDCHForders", style = list(display = "none")),
    htmlDiv(id = "USDJPYorders", style = list(display = "none")),
    htmlDiv(id = "GBPUSDorders", style = list(display = "none"))
  )
}

ordersRows <- function(list_order, st){
  headers <- list(
    "Order ID",
    "Time",
    "Type",
    "Volume",
    "Symbol",
    "TP",
    "SL",
    "Price",
    "Profit",
    "Status"
  )
  if (st == "closed"){
    headers <- c(headers, list("Close Time", "Close Price"))
  }
  headerList <- htmlTr(lapply(headers, htmlTh))
  if (!is.null(list_order)){
    if (sum(list_order[, "status"] == st) == 0){
      rows <- list(
        htmlTr(
          htmlTd(
            sprintf("No %s positions data available", st),
            colSpan = length(headers),
            style = list(
              textAlign = "center",
              color = "#b4b4b4"
            )
          )
        )
      )
    } else {
      rows <- lapply(
        seq(1, nrow(list_order)),
        function(i){
          if (list_order[i, "status"] == st){
            htmlTr(
              lapply(
                as.character(list_order[i, 1:length(headers)]),
                htmlTd
              ),
              style = list(
                background = ifelse(
                  list_order[i, "profit"] < 0,
                  "rgba(255,65,54,0.5)",
                  "rgba(61,153,112,0.5)"
                ),
                borderBottom = "1px solid #262a30"
              )
            )
          }
        }
      )
    }
    return(c(list(headerList), rows))
  }
  rows <- list(
    htmlTr(
      htmlTd(
        "No orders have been placed yet",
        colSpan = length(headers),
        style = list(
          textAlign = "center",
          color = "#b4b4b4"
        )
      )
    )
  )
  return(c(list(headerList), rows))
}

updateOrders <- function(orders, current_bids, current_asks, id_to_close){
  currencies <- c("EURUSD", "USDCHF", "USDJPY", "GBPUSD")
  o <- as.data.table(orders)
	olist <- data.table()
  for (i in 1:nrow(o)){
    if (o[i, "status"] == "open"){
      current_bid <- current_bids[match(o[i, "symbol"], currencies)]
      current_ask <- current_asks[match(o[i, "symbol"], currencies)]
      if (o[i, "type"] == "buy"){
        o[, profit :=
          as.numeric(profit)][
        i, profit :=
          round(
            o[i, "volume"] *
              100000 *
              ( (current_bid - o[i, "price"]) / o[i, "price"]),
            2
          )
        ]
      } else {
        o[, profit :=
          as.numeric(profit)][
        i, profit :=
          round(
            o[i, "volume"] *
              100000 *
              ( (o[i, "price"] - current_ask) / o[i, "price"]),
            2
          )
        ]
      }
      if (o[i, "type"] == "buy"){
        price <- current_bid
      } else {
        price <- current_ask
      }
      if (!is.null(id_to_close)){
        if (o[i, "id"] == id_to_close){
          o[i, status := "closed"]
          o[i, "close Time" := strftime(
              Sys.time(), format = "%Y-%m-%d %H:%M:%S"
              )
          ]
          o[i, "close Price" := price]
        }
      }
      if (o[i, "tp"] != 0 & price >= o[i, "tp"]){
        o[i, status := "closed"]
        o[i, "close Time" :=
          strftime(Sys.time(), format = "%Y-%m-%d %H:%M:%S")
        ]
        o[i, "close Price" := price]
      }
      if (o[i, "sl"] != 0 & price >= o[i, "sl"]){
        o[i, status := "closed"]
        o[i, "close Time" :=
          strftime(Sys.time(), format = "%Y/%m/%d %H:%M:%S")
        ]
        o[i, "close Price" := price]
      }
    }
    olist <- rbindlist(list(olist, o[i, ]), fill = TRUE)
  }
  return(olist)
}
