library(shiny)
library(plotly)

shinyUI(fluidPage(
    titlePanel("Movie Ratings!"),
    sidebarPanel(
        sliderInput("bins", "Number of bins:", min = 1, max = 50, value = 10)
    ),
    mainPanel(
        plotlyOutput("trendPlot")
    )
))