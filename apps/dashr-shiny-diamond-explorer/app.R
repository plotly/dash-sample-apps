library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
#library(dplyr) #package for data manipulation

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf("/app/apps/%s", appName))
}

data(diamonds, package = "ggplot2")
nms <- names(diamonds)

#main application
external_ss = list('https://codepen.io/chriddyp/pen/bWLwgP.css')
app <- Dash$new(external_stylesheets=external_ss)
#myList = list()
#nms
lapply(nms, function(x){list(label = x, value = x)})

app$layout(
  htmlDiv(htmlImg(src='assets/dash-logo.png')),
  htmlDiv(htmlH1("Diamonds Explorer")),
  htmlDiv(htmlLabel("Sample Size")),
  htmlDiv(
    dccSlider(min = 1, max=nrow(diamonds), step = 500, 
      marks = as.list(
        setNames(
          round(seq(1, nrow(diamonds), length.out = 10),0), 
          round(seq(1,nrow(diamonds),length.out = 10), 0)
          )
        ), value = 1000)
      ),
  htmlBr(), htmlBr(),
  htmlDiv(htmlLabel("X")),
  htmlDiv(
    dccDropdown(options = lapply(
      nms, function(x) {list(label = x, value = x)}
                                 ),
                value = "carat")
    ),
  htmlBr(),
  htmlDiv(htmlLabel("Y")),
  htmlDiv(dccDropdown(options = lapply(nms, function(x) {list(label = x, value = x)}), value = "price")),
  htmlBr(),
  htmlDiv(htmlLabel("Color")), 
  htmlDiv(dccDropdown(options = lapply(nms, function(x) {list(label = x, value = x)}), value="clarity")),
  htmlBr(),
  htmlDiv(htmlLabel("Facet Row")),
  htmlDiv(dccDropdown(options = lapply(nms, function(x) {list(label = x, value = x)}), value="clarity")),
  htmlBr(),
  htmlDiv(htmlLabel("Facet Column")),
  dccDropdown(options = lapply(nms, function(x) {list(label = x, value = x)}), value="None"),
  htmlBr(),
  # include marks for R ---- For slider components
  htmlDiv(htmlLabel("Height of plot (in pixels)")),
  
  htmlDiv(
    list(
      htmlLabel("Height of plot (in pixels)"),
      dccSlider(min= 100, max = 2000, marks = c("100", "290", "480", "670", "860", "1,240", "1,620", "2000"), value = "1000")
    )
  ),
  
  htmlBr()
)
  
#app$callback(
  # build graph with ggplot syntax
   #p <- ggplot(dataset(), aes_string(x = input$x, y = input$y, color = input$color)) + 
            #geom_point()
   # if at least one facet column/row is specified, add it
        #facets <- paste(input$facet_row, '~', input$facet_col)
       # if (facets != '. ~ .') p <- p + facet_grid(facets)
        
        #ggplotly(p) %>% 
       #     layout(height = input$plotHeight, autosize=TRUE))

#adding conditional statement for running app on server and local

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server(showcase = TRUE)
}
