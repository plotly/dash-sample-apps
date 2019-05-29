library(glue)
library(dashR)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dashTable)
library(data.table)
library(plotly)
source("lastodf.R")

app <- Dash$new(name='DashR LAS Report')

LASdir <- "data/"
LASfile <- "alcor1.las"
LASpath <- paste(LASdir, LASfile, sep="")
LASdata <- convertLAS(LASpath)

generate_frontpage <- function() {
  desc <- LASdata$V[[3]][1]
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
                  htmlB(id="las-filename", children=LASfile),
                  htmlSpan(glue("({desc})"))
                )
              )
            )
          )
        )
      )
    )
  )
}

generate_axis_title <- function(desc, unit) {
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

plotList <- function(nplots) {
  lapply(seq_len(nplots), function(x) plot_ly())
}

generate_curves <- function(height=950,
                            width=725,
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
    p <- plot_ly(LASdata$A,y=LASdata$A[[yval]])
    for (xval in xvals[[i]]) {
      p <- add_trace(p, x=LASdata$A[[xval]], name=xval, mode='lines', 
                     line=list(dash=line_style, width=line_width))
    }
    plots[length(plots)+1] <- p
  }
  fig <- subplot(
    plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], 
    plotList(5), nrows=1, shareY=TRUE, margin=0
  )
  #fig["data"][2]["xaxis"] = "x6"
  #fig["data"][7]["xaxis"] = "x7"
  #fig["data"][9]["xaxis"] = "x8"
  #fig["data"][12]["xaxis"] = "x9"
  fig <- layout(fig, 
    height = height, 
    width = width, 
    plot_bgcolor = bg_color, 
    paper_bgcolor = bg_color, 
    hovermode = "y", 
    legend = list(font = list(size = tick_font_size)), 
    margin = list(r=100),
    yaxis = list(
      title = list(text=generate_axis_title("Depth", "F"), font=list(family="Arial, sans-serif", size=font_size)), 
      autorange = "reversed",
      mirror = "all",
      automargin = TRUE,
      showline = TRUE,
      tickfont = list(family="Arial, sans-serif", size=tick_font_size)
    ),
    xaxis = list(
      title = list(text="TITLE 1", font=list(family="Arial, sans-serif", size=font_size)), 
      mirror = "all",
      automargin = TRUE,
      showline = TRUE,
      tickfont = list(family="Arial, sans-serif", size=tick_font_size)
    ),
    xaxis2 = list(
      title = list(text="TITLE 2", font=list(family="Arial, sans-serif", size=font_size)), 
      mirror = "all",
      automargin = TRUE,
      showline = TRUE,
      tickfont = list(family="Arial, sans-serif", size=tick_font_size),
      type="log"
    ),
    xaxis3 = list(
      title = list(text="TITLE 3", font=list(family="Arial, sans-serif", size=font_size)), 
      mirror = "all",
      automargin = TRUE,
      showline = TRUE,
      tickfont = list(family="Arial, sans-serif", size=tick_font_size)
    ),
    xaxis4 = list(
      title = list(text="TITLE 4", font=list(family="Arial, sans-serif", size=font_size)), 
      mirror = "all",
      automargin = TRUE,
      showline = TRUE,
      tickfont = list(family="Arial, sans-serif", size=tick_font_size)
    ),
    xaxis5 = list(
      title = list(text="TITLE 5", font=list(family="Arial, sans-serif", size=font_size)), 
      mirror = "all",
      automargin = TRUE,
      showline = TRUE,
      tickfont = list(family="Arial, sans-serif", size=tick_font_size)
    ),
    xaxis6 = list(
      title = list(text="TITLE 6", font=list(family="Arial, sans-serif", size=font_size)),
      overlaying="x1",
      anchor="y",
      side="top",
      mirror = "all",
      automargin = TRUE,
      showline = TRUE,
      tickfont = list(family="Arial, sans-serif", size=tick_font_size)
    ),
    xaxis7 = list(
      title = list(text="TITLE 7", font=list(family="Arial, sans-serif", size=font_size)),
      overlaying="x1",
      anchor="y",
      side="top",
      mirror = "all",
      automargin = TRUE,
      showline = TRUE,
      tickfont = list(family="Arial, sans-serif", size=tick_font_size)
    )
  )
  return(dccGraph(figure=fig))
}

generate_table <- function() {
  #columns <- list()
  #for (key in names(LASdata$W)) {
  #  columns[[length(columns)+1]] <- list(name = key, id = key)
  #}
  cols <- list (
    list(name = "mnemonic", id = "mnemonic"),
    list(name = "description", id = "description"),
    list(name = "unit", id = "unit"),
    list(name = "value", id = "value")
  )
  return(
    dashTable::dashDataTable(
      id = "table",
      sorting = TRUE,
      filtering = TRUE,
      row_deletable = TRUE,
      style_cell = list(padding="5px", width="auto", textAlign="left"),
      columns = cols,
      data = dashTable::df_to_list(LASdata$W)
    )
  )
}

app$layout(
  htmlDiv(
    list(
      htmlDiv(id="controls", children=list(htmlButton("Print", id="las-print"))),
      htmlDiv(id="frontpage", className="page", children=generate_frontpage()),
      htmlDiv(className="section-title", children="LAS well"),
      htmlDiv(
        className="page",
        children=list(
          htmlDiv(id="las-table", children=generate_table()),
          htmlDiv(id="las-table-print")
        ),
      ),
      htmlDiv(className="section-title", children="LAS curves"),
      htmlDiv(
        className="page",
        children=list(htmlDiv(id="las-curves", children=generate_curves()))
      )
    )
  )
)

app$callback(
  output = list(id="las-table-print", property="children"),
  params = list(input(id="table", property="data")),
  function (data) {
    colwidths = list(
      mnemonic = "100px",
      descr = "300px",
      unit = "25px",
      value = "300px"
    )
    table_list <- list()
    num_tables = length(data)/34+1 # 34 rows max per page
    for (i in 1:num_tables) {
      table_rows <- list()
      Th <- list()
      for (key in names(data[[1]])) {
        Th[[length(Th)+1]] <- htmlTh(key, style=list(width=colwidths[[key]]))
      }
      table_rows[[1]] <- htmlTr(Th)
      for (j in 1:34) {
        if (i*34+j <= length(data)) {
          break
        }
        Td <- list()
        for (key in names(data[[1]])) {
          Td[[length(Td)+1]] <- htmlTd(data[i*34+j][key])
        }
        table_rows[[length(table_rows)+1]] <- htmlTr(Td)
      }
    }
    table_list[[length(table_list)+1]] <- htmlDiv(className="tablepage", children=htmlTable(table_rows))
    return (table_list)
  }
)

#print("TEST FUNCTION (generate axis):")
#print(generate_axis_title("BAT comp slowness", "uspf"))

app$run_server(showcase=TRUE)
