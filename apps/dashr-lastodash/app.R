library(readr)
library(stringr)
library(gdata)
library(tidyr)
library(dplyr)
library(glue)
library(ggplot2)
library(dashR)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dashTable)
library(data.table)
library(plotly)

########################################################################################################################

las_object <- setRefClass("las_object", fields = list(
  V="data.frame", 
  W="data.frame", 
  C="data.frame", 
  P="data.frame", 
  O="data.frame", 
  A="data.frame")
)

convert_las <- function(las_path) {
  las_data <- las_object()
  las_string <- read_file(las_path)
  las_sections <- str_split(las_string, "~")
  for (section in las_sections[[1]]) {
    # section contains version and wrap mode information
    if (startsWith(section, "V")) {
      print("~V")
      linesV <- str_split(section, "\r\n")
      contentV <- linesV[[1]][-c(1, length(linesV[[1]]))]
      las_data$V <- data.frame(x = contentV) %>% 
        separate(col=x, into=c("label", "x"), sep="\\.", extra="merge") %>% 
        separate(col=x, into=c("value", "description"), sep=":\\s+", extra="merge")
    }
    # section contains well identification
    if (startsWith(section, "W")) {
      print("~W")
      linesW <- str_split(section, "\r\n")
      contentW <- linesW[[1]][-c(1, length(linesW[[1]]))]
      if ((length(contentW) >= 2) & startsWith(contentW, "#")) {
        contentW <- contentW[-c(1, 2)]
      }
      #las_data$W <- read.table(text=contentW, header=FALSE, sep="", strip.white=TRUE) %>%  
      las_data$W <- data.frame(x = contentW) %>%
        separate(col=x, into=c("mnemonic", "x"), sep="\\.", extra="merge") %>% 
        separate(col=x, into=c("unit", "x"), sep="\\s+", extra="merge") %>% 
        separate(col=x, into=c("value", "description"), sep=":\\s+", extra="merge")
    }
    # section contains curve information
    if (startsWith(section, "C")) {
      print("~C")
      linesC <- str_split(section, "\r\n")
      contentC <- linesC[[1]][-c(1, length(linesC[[1]]))]
      if ((length(contentC) >= 2) & startsWith(contentC, "#")) {
        contentC <- contentC[-c(1, 2)]
      }
      las_data$C <- data.frame(x = contentC) %>%
        separate(col=x, into=c("mnemonic", "x"), sep="\\.", extra="merge") %>%
        separate(col=x, into=c("unit", "x"), sep="\\s+", extra="merge") %>%
        separate(col=x, into=c("api code", "curve description"), sep=":\\s+", extra="merge")
    }
    # section contains parameters or constants
    if (startsWith(section, "P")) {
      print("~P")
    }
    # section contains other information such as comments
    if (startsWith(section, "O")) {
      print("~O")
    }
    # section contains ASCII log data
    if (startsWith(section, "A")) {
      print("~A")
      #linesA <- str_split(section, "\r\n")
      #contentA <- linesA[[1]][-c(1, length(linesC[[1]]))]
      #tableA <- read.table(text=contentA, header=FALSE, sep="", col.names=as.vector(tableC['mnemonic']))
      contentA <- str_sub(section, 2, -2)
      las_data$A <- read.table(text=contentA, header=TRUE, sep="")
      #las_data$A <- as.data.frame(apply(dataA,2,function(x)str_trim(x,"left")))
      #las_data$A <- sub("[[:space:]]+$","",read.table(text=contentA, header=TRUE, sep=""))
    }
  }
  return (las_data)
}

########################################################################################################################

las_dir <- "data/"
las_file <- "alcor2.las"
las_path <- paste(las_dir, las_file, sep="")
las_data <- convert_las(las_path)

table_cols <- list (
  "mnemonic" = list(name="mnemonic", id="mnemonic", width="100px"),
  "description" = list(name="description", id="description", width="300px"),
  "unit" = list(name="unit", id="unit", width="25px"),
  "value" = list(name="value", id="value", width="300px")
)

########################################################################################################################

get_item <- function(df, name1, name2, key1) {
  key2 <- ""
  for (i in 1:nrow(df)) {
    if (str_trim(df[i,name1],"left")==key1) {
      key2 <- df[i,name2]
    }
  }
  return (key2)
}

get_plot_list <- function(nplots) {
  lapply(seq_len(nplots), function(x) plot_ly())
}

########################################################################################################################

generate_frontpage <- function() {
  desc <- str_trim(las_data$V[[3]][1], "left")
  return(
    list(
      htmlDiv(
        id="las-header",
        children=list(
          htmlImg(id="las-logo", src="assets/logo.png"),
          htmlDiv(
            id="las-header-text",
            children=list(
              htmlH1("LAS Report"),
              htmlDiv(
                id="las-file-info",
                children=list(
                  htmlSpan(id="las-filename", children=las_file),
                  htmlSpan(glue(" ({desc})"))
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
    dashTable::dashDataTable(
      id="table",
      sorting = TRUE,
      filtering = TRUE,
      row_deletable = TRUE,
      style_cell = list(
        padding="15px", 
        width="auto", 
        textAlign="center", 
        border="white"
      ),
      style_cell_conditional = list(list(
        'if'=list(row_index="even"), 
        backgroundColor="#f9f9f9"
      )),
      columns = unname(table_cols),
      data = dashTable::df_to_list(las_data$W[c("mnemonic", "description", "unit", "value")])
    )
  )
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
    ifelse (
      identical(curLine, ""), 
      curLine <- paste(curLine, word, sep=""), 
      curLine <- paste(curLine, word, sep=" ")
    ) 
  }
  lines[[length(lines)+1]] <- curLine
  title <- paste(lines, collapse="\n")
  title <- glue("{title}\n({unit})")
  return(title)
}

generate_curves <- function(height=950,
                            width=800,
                            bg_color="white",
                            font_size=10,
                            tick_font_size=8,
                            line_width=0.5) {
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
    ifelse(i==2, line_style <- 'dashdot', line_style <- 'solid')
    p <- plot_ly(las_data$A, y=las_data$A[[yval]])
    for (xval in xvals[[i]]) {
      p <- add_trace(p, x=las_data$A[[xval]], name=xval, mode='lines', type="scatter",
                     line=list(dash=line_style, width=line_width))
    }
    plots[length(plots)+1] <- p
  }
  fig <- subplot(
    plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], 
    get_plot_list(5), nrows=1, shareY=TRUE, margin=0
  )
  fig$x$data[[2]][["xaxis"]] = "x6"
  fig$x$data[[7]][["xaxis"]] = "x7"
  fig$x$data[[9]][["xaxis"]] = "x8"
  fig$x$data[[12]][["xaxis"]] = "x9"
  for (i in 1:length(xvals)) {
    for (xval in xvals[[i]]) {
      ifelse(i==1, xaxis <- "xaxis", xaxis <- glue("xaxis{i}"))
      ifelse(i==2, fig$x$layout[[xaxis]]$type <- "log", fig$x$layout[[xaxis]]$type <- "linear")
    }
  }
  for (i in c(6:9)) {
    xaxis <- glue("xaxis{i}")
    ifelse(i!=9, fig$x$layout[[xaxis]]$overlaying <- glue("x{i-5}"), fig$x$layout[[xaxis]]$overlaying <- glue("x{i-4}"))
    fig$x$layout[[xaxis]]$anchor <- "y"
    fig$x$layout[[xaxis]]$side <- "top"
    fig$x$layout[[xaxis]]$title <- glue("TITLE {i}")
  }
  for (i in 1:9) {
    ifelse(i %in% c(1:5), key <- xvals[[i]][[1]], ifelse(i!=9, key <- xvals[[i-5]][[length(xvals[[i-5]])]], key <- xvals[[i-4]][[length(xvals[[i-4]])]]))
    ifelse(i==1, xaxis <- "xaxis", xaxis <- glue("xaxis{i}"))
    fig$x$layout[[xaxis]]$mirror <- "all"
    fig$x$layout[[xaxis]]$automargin <- TRUE
    fig$x$layout[[xaxis]]$showline <- TRUE
    fig$x$layout[[xaxis]]$title <- list(text=generate_axis_title(las_data$C, key), font=list(family="Arial, sans-serif", size=font_size))
    fig$x$layout[[xaxis]]$tickfont <- list(family="Arial, sans-serif", size=tick_font_size)
  }
  # y-axis
  fig$x$layout[["yaxis"]]$mirror <- "all"
  fig$x$layout[["yaxis"]]$automargin <- TRUE
  fig$x$layout[["yaxis"]]$showline <- TRUE
  fig$x$layout[["yaxis"]]$title <- list(text=generate_axis_title(las_data$C, "DEPT"), font=list(family="Arial, sans-serif", size=font_size))
  fig$x$layout[["yaxis"]]$tickfont <- list(family="Arial, sans-serif", size=tick_font_size)
  fig$x$layout[["yaxis"]]$autorange <- 'reversed'
  # layout
  fig$x$layout$height <- height
  fig$x$layout$width <- width
  fig$x$layout$plot_bgcolor <- bg_color
  fig$x$layout$paper_bgcolor <- bg_color
  fig$x$layout$hovermode <- "y"
  #fig$x$layout$legend <- list("xpad"=0.0,"font"=list("size"=tick_font_size))
  fig$x$layout$legend <- list("x"=0.85,"xpad"=0.0,"font"=list("size"=tick_font_size))
  fig$x$layout$margin <- list(r=0, pad=0, autoexpand=FALSE)
  
  return(dccGraph(figure=fig))
}

########################################################################################################################

app <- Dash$new(name='DashR LAS Report')

app$layout(
  htmlDiv(
    list(
      htmlDiv(
        id="controls", 
        children=list(
          htmlButton(
            "Print", 
            id="las-print"
          )
        )
      ),
      htmlDiv(
        id="frontpage",
        className="page",
        children=generate_frontpage()
      ),
      htmlDiv(
        className="section",
        children=list(
          htmlDiv(
            className="section-title",
            children="LAS well"
          ),
          htmlDiv(
            className="page",
            children=list(
              htmlDiv(
                id="las-table", 
                children=generate_table()
              ),
              htmlDiv(
                id="las-table-print"
              )
            )
          )
        )
      ),
      htmlDiv(
        className="section",
        children=list(
          htmlDiv(
            className="section-title", 
            children="LAS curves"
          ),
          htmlDiv(
            className="page",
            children=list(
              htmlDiv(
                id="las-curves", 
                children=generate_curves()
              )
            )
          )
        )
      )
    )
  )
)

app$callback(
  output = list(id="las-table-print", property="children"),
  params = list(input(id="table", property="data")),
  function (data) {
    #print("------------------ PRINT DATA ------------------")
    #str(data[[2]]['mnemonic'])
    #str(data)
    table_list <- list()
    num_tables = as.integer(length(data)/34)+1 # 34 rows max per page
    for (i in 1:num_tables) {
      table_rows <- list()
      Th <- list()
      for (key in names(data[[1]])) {
        Th[[length(Th)+1]] <- htmlTh(str_to_title(key), style=list(width=table_cols[[key]]$width))
      }
      table_rows[[1]] <- htmlTr(Th)
      for (j in 1:34) {
        if ((i-1)*34+j > length(data)) {
          break
        }
        Td <- list()
        for (key in names(data[[1]])) {
          #print("------------------ PRINT NOW ------------------")
          #str(data[[(i-1)*34+j]])
          Td[[length(Td)+1]] <- htmlTd(data[[(i-1)*34+j]][key])
        }
        table_rows[[j+1]] <- htmlTr(Td)
      }
      #print("------------------ PRINT TABLE i ------------------")
      #str(table_rows[[2]])
    }
    table_list[[i]] <- htmlDiv(className="tablepage", children=htmlTable(table_rows))
    return (table_list)
  }
)

app$run_server(showcase=TRUE)
