library(dashR)
library(dashHtmlComponents)
library(dashCoreComponents)

app <- Dash$new()

build_banner <- function() {
  return (
    htmlDiv(
      id = 'banner',
      className = 'banner',
      list(
        htmlH5('Manufacturing SPC Dashboard - Process Control and Exception Reporting'),
        htmlButton(
          id = 'learn-more-button',
          list('LEARN MORE'),
          n_clicks = 0
        ),
        htmlImg(
          src = 'https://s3-us-west-1.amazonaws.com/plotly-tutorials/logo/new-branding/dash-logo-by-plotly-stripe-inverted.png'
        )
      )
    )
  )
}

build_tabs <- function() {
  return ('')
}

init_df <- function() {
  return ('')
}

populate_ooc <- function() {
  return ('')
}

init_value_setter_store <- function() {
  return ('')
}

build_tab_1 <- function() {
  return ('')
}

build_value_setter_line <- function() {
  return ('')
}

generate_modal <- function() {
  return ('')
}

app$layout(
  htmlDiv(list(
    build_banner(),
    build_tabs(),
    htmlDiv(
      id = 'app-content',
      className = 'container scalable'
    ),
    dccStore(
      id = 'value-setter-store',
      data = init_value_setter_store()
    ),
    generate_model()
  ))
)

app$run_server(showcase = TRUE)