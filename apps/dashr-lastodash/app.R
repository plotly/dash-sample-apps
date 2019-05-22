library(argparse)
library(dashR)
library(dashHtmlComponents)
library(dashCoreComponents)
library(data.table)

app <- Dash$new()

parse_args <- function() {
  
  return('')
}

generate_frontpage <- function() {
  return('')
}

generate_axis_title <- function() {
  return('')
}

generate_curves <- function() {
  return('')
}

generate_table <- function() {
  return('')
}

app$layout(
  htmlDiv(
    list(
      htmlDiv(id="controls", list(htmlButton("Print", id="las-print"))),
      htmlDiv(id="frontpage", className="page", list(generate_frontpage())),
      htmlDiv(className="section-title", list("LAS well")),
      htmlDiv(
        className="page",
        list(
          htmlDiv(id="las-table", list(generate_table()))
          #htmlDiv(id="las-table-print"),
        ),
      ),
      htmlDiv(className="section-title", list("LAS curves")),
      htmlDiv(
        className="page",
        list(htmlDiv(id="las-curves", list(generate_curves())))
      )
    )
  )
)

app$run_server(showcase = TRUE)
