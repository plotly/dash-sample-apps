library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("/app/apps/%s", appName))
}

app <- Dash$new()

app$layout(
  htmlDiv(
    list(
      htmlDiv("This is an app"),
      dccGraph(
        figure = list(
          data = list(
            list(
              x = c(1,2,3),
              y = c(4, 6, 2)
            )
          )
        )
      )
    )
  )
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
