library(readr)
library(stringr)
library(gdata)
library(tidyr)
library(dplyr)
library(dash)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dashTable)
library(plotly)

####################################################################################################
# INITIATE APP

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  Sys.setenv(
    DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
    DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix
  )
}

app <- Dash$new(name = 'DashR LAS Report')

####################################################################################################
# IMPORT DATA

las_dir <- "data/"
las_file <- "alcor2.las"
las_path <- paste(las_dir, las_file, sep = "")

####################################################################################################
# HANDLE & STORE DATA

convert_las <- function(las_path) {
  las_data <- list()
  las_string <- read_file(las_path)
  las_sections <- str_split(las_string, "~")
  for (section in las_sections[[1]]) {
    # store section containing version and wrap mode information
    if (startsWith(section, "V")) {
      linesV <- str_split(section, "\r\n")
      contentV <- linesV[[1]][-c(1, length(linesV[[1]]))]
      las_data$V <- data.frame(x = contentV) %>%
        separate(col = x, into = c("label", "x"), sep = "\\.", extra = "merge") %>%
        separate(col = x, into = c("value", "description"), sep = ":\\s+", extra = "merge")
    }
    # store section containing well identification
    if (startsWith(section, "W")) {
      linesW <- str_split(section, "\r\n")
      contentW <- linesW[[1]][-c(1, length(linesW[[1]]))]
      if ((length(contentW) >=  2) & startsWith(contentW, "#")) {
        contentW <- contentW[-c(1, 2)]
      }
      las_data$W <- data.frame(x = contentW) %>%
        separate(col = x, into = c("mnemonic", "x"), sep = "\\.", extra = "merge") %>%
        separate(col = x, into = c("unit", "x"), sep = "\\s+", extra = "merge") %>%
        separate(col = x, into = c("value", "description"), sep = ":\\s+", extra = "merge")
    }
    # store section containing curve information
    if (startsWith(section, "C")) {
      linesC <- str_split(section, "\r\n")
      contentC <- linesC[[1]][-c(1, length(linesC[[1]]))]
      if ((length(contentC) >=  2) & startsWith(contentC, "#")) {
        contentC <- contentC[-c(1, 2)]
      }
      las_data$C <- data.frame(x = contentC) %>%
        separate(col = x, into = c("mnemonic", "x"), sep = "\\.", extra = "merge") %>%
        separate(col = x, into = c("unit", "x"), sep = "\\s+", extra = "merge") %>%
        separate(col = x, into = c("api code", "curve description"), sep = ":\\s+", extra = "merge")
    }
    # section contains ASCII log data
    if (startsWith(section, "A")) {
      contentA <- str_sub(section, 2, -2)
      las_data$A <- read.table(text = contentA, header = TRUE, sep = "")
    }
  }
  return (las_data)
}

las_data <- convert_las(las_path)

####################################################################################################
# DEFINE GLOBAL VARIABLES

table_header <- list (
  "mnemonic" = list(name = "mnemonic", id = "mnemonic", width = "10%"),
  "description" = list(name = "description", id = "description", width = "40%"),
  "unit" = list(name = "unit", id = "unit", width = "10%"),
  "value" = list(name = "value", id = "value", width = "40%")
)

####################################################################################################
# DEFINE FUNCTIONS FOR APP LAYOUT

generate_frontpage <- function() {
  desc <- str_trim(las_data$V[[3]][1], "left")
  return(
    list(
      htmlDiv(
        id = "las-header",
        children = list(
          htmlImg(id = "las-logo", src = "assets/logo.png"),
          htmlDiv(
            id = "las-header-text",
            children = list(
              htmlH1("LAS Report"),
              htmlDiv(
                id = "las-file-info",
                children = list(
                  htmlSpan(id = "las-filename", children = las_file),
                  htmlSpan(sprintf(" (%s)", desc))
                )
              )
            )
          )
        )
      )
    )
  )
}

generate_table <- function() {
  return(
    dashDataTable(
      id = "table",
      sort_action = 'native',
      filter_action = 'native',
      row_deletable = TRUE,
      css = list(list(
        selector = '.dash-cell div.dash-cell-value',
        rule = 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'
      )),
      style_data = list(
        whiteSpace = "normal"
      ),
      style_cell = list(
        padding = "15px",
        minWidth = "0px",
        width = "25%",
        textAlign = "center",
        border = "white"
      ),
      style_cell_conditional = list(list(
        'if' = list(row_index = "even"),
        backgroundColor = "#f9f9f9"
      )),
      columns = unname(table_header),
      data = df_to_list(las_data$W[c("mnemonic", "description", "unit", "value")])
    )
  )
}

get_item <- function(df, name1, name2, key1) {
  key2 <- ""
  for (i in 1:nrow(df)) {
    if (str_trim(df[i,name1],"left") == key1) {
      key2 <- df[i,name2]
    }
  }
  return (key2)
}

generate_axis_title <- function(df, key) {
  desc <- get_item(df, "mnemonic", "curve description", key)
  unit <- get_item(df, "mnemonic", "unit", key)
  axisWords <- str_split(desc, " ")
  curLine <- ""
  lines <- list()
  for (word in axisWords[[1]]) {
    # create new line if needed
    if ((nchar(curLine) + nchar(word)) > 15) {
      lines[[length(lines)+1]] <- curLine
      curLine <- ""
    }
    # concatenate strings
    curLine <- ifelse (
      identical(curLine, ""),
      paste(curLine, word, sep = ""),
      paste(curLine, word, sep = " ")
    )
  }
  lines[[length(lines)+1]] <- curLine
  title <- paste(lines, collapse = "\n")
  title <- sprintf("%s\n(%s)", title, unit)
  return(title)
}

get_plot_list <- function(nplots) {
  lapply(seq_len(nplots), function(x) plot_ly())
}

generate_curves <- function(height = 950,
                            width = 800,
                            bg_color = "white",
                            font_size = 10,
                            tick_font_size = 8,
                            line_width = 0.5) {
  yval <- "DEPT"
  xvals <- list(
    list("BTVPVS", "DGRC"),
    list("R15P", "R09P", "R39P", "R27P", "EWXT"),
    list("ALCDLC", "ALDCLC"),
    list("TNPS"),
    list("BTCSS", "BTCS")
  )
  plots <- list()
  for (i in 1:length(xvals)) {
    line_style <- ifelse(i == 2, 'dashdot', 'solid')
    p <- plot_ly(las_data$A, y = las_data$A[[yval]])
    for (xval in xvals[[i]]) {
      p <- add_trace(p, x = las_data$A[[xval]], name = xval, mode = 'lines', type = "scatter",
                     line = list(dash = line_style, width = line_width))
    }
    plots[length(plots)+1] <- p
  }
  fig <- subplot(
    plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
    get_plot_list(5), nrows = 1, shareY = TRUE, margin = 0
  )
  fig$x$data[[2]][["xaxis"]] <- "x6"
  fig$x$data[[7]][["xaxis"]] <- "x7"
  fig$x$data[[9]][["xaxis"]] <- "x8"
  fig$x$data[[12]][["xaxis"]] <- "x9"
  for (i in 1:length(xvals)) {
    for (xval in xvals[[i]]) {
      xaxis <- ifelse(i == 1, "xaxis", sprintf("xaxis%i", i))
      fig$x$layout[[xaxis]]$type <- ifelse(i == 2, "log", "linear")
    }
  }
  for (i in c(6:9)) {
    xaxis <- sprintf("xaxis%i", i)
    fig$x$layout[[xaxis]]$overlaying <- ifelse(i != 9, sprintf("x%i", i-5), sprintf("x%i", i-4))
    fig$x$layout[[xaxis]]$anchor <- "y"
    fig$x$layout[[xaxis]]$side <- "top"
  }
  for (i in 1:9) {
    key <- ifelse(i %in% c(1:5), xvals[[i]][[1]], ifelse(i !=  9, xvals[[i-5]][[length(xvals[[i-5]])]], xvals[[i-4]][[length(xvals[[i-4]])]]))
    ifelse(i == 1, xaxis <- "xaxis", xaxis <- sprintf("xaxis%i", i))
    fig$x$layout[[xaxis]]$mirror <- "all"
    fig$x$layout[[xaxis]]$automargin <- TRUE
    fig$x$layout[[xaxis]]$showline <- TRUE
    fig$x$layout[[xaxis]]$title <- list(text = generate_axis_title(las_data$C, key), font = list(family = "Arial, sans-serif", size = font_size))
    fig$x$layout[[xaxis]]$tickfont <- list(family = "Arial, sans-serif", size = tick_font_size)
  }
  # y-axis
  fig$x$layout[["yaxis"]]$mirror <- "all"
  fig$x$layout[["yaxis"]]$automargin <- TRUE
  fig$x$layout[["yaxis"]]$showline <- TRUE
  fig$x$layout[["yaxis"]]$title <- list(text = generate_axis_title(las_data$C, "DEPT"), font = list(family = "Arial, sans-serif", size = font_size))
  fig$x$layout[["yaxis"]]$tickfont <- list(family = "Arial, sans-serif", size = tick_font_size)
  fig$x$layout[["yaxis"]]$autorange <- 'reversed'
  # layout
  fig$x$layout$height <- height
  fig$x$layout$width <- width
  fig$x$layout$plot_bgcolor <- bg_color
  fig$x$layout$paper_bgcolor <- bg_color
  fig$x$layout$hovermode <- "y"
  fig$x$layout$legend <- list("font" = list("size" = tick_font_size))
  fig$x$layout$margin <- list(r = 0.0, t = 100.0, b = 50.0, l = 80.0, autoexpand = FALSE)
  return(dccGraph(figure = fig))
}

####################################################################################################
# DEFINE APP LAYOUT

app$layout(
  htmlDiv(
    list(
      htmlDiv(
        id = "controls",
        children = list(
          htmlButton(
            "Print",
            id = "las-print"
          )
        )
      ),
      htmlDiv(
        id = "frontpage",
        className = "page",
        children = generate_frontpage()
      ),
      htmlDiv(
        className = "section",
        children = list(
          htmlDiv(
            className = "section-title",
            children = "LAS well"
          ),
          htmlDiv(
            className = "page",
            children = list(
              htmlDiv(
                id = "las-table",
                children = generate_table()
              ),
              htmlDiv(
                id = "las-table-print"
              )
            )
          )
        )
      ),
      htmlDiv(
        className = "section",
        children = list(
          htmlDiv(
            className = "section-title",
            children = "LAS curves"
          ),
          htmlDiv(
            className = "page",
            children = list(
              htmlDiv(
                id = "las-curves",
                children = generate_curves()
              )
            )
          )
        )
      )
    )
  )
)

####################################################################################################
# DEFINE APP CALLBACKS

app$callback(
  output = list(id = "las-table-print", property = "children"),
  params = list(input(id = "table", property = "data")),
  function (data) {
    num_rows <- 34
    num_tables <- as.integer(length(data)/num_rows)+1
    output <- list()
    Th <- lapply(
      names(data[[1]]),
      function(s) {return(htmlTh(str_to_title(s), style = list(width = table_header[[s]]$width)))}
    )
    header <- list(htmlTr(Th))
    for (i in 1:num_tables) {
      start <- (i-1)*num_rows+1
      end <- ifelse (i*num_rows > length(data), length(data), i*num_rows)
      rows <- lapply(
        data[start:end],
        function(L) {
          Tr <- lapply(
            unname(L),
            function(s) {return(htmlTd(s))}
          )
          return(htmlTr(Tr))
        }
      )
      output[[i]] <- htmlDiv(className = "tablepage", children = htmlTable(c(header, rows)))
    }
    return (output)
  }
)

####################################################################################################
# RUN SERVER

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()
}

