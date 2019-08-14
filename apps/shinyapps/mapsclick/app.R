library(shiny)
library(plotly)

ui <- fluidPage(
    plotlyOutput("plot"),
    verbatimTextOutput("click")
)

server <- function(input, output, session) {
    
    output$plot <- renderPlotly({
        # specify some map projection/options
        g <- list(
            scope = 'usa',
            projection = list(type = 'albers usa'),
            lakecolor = toRGB('white')
        )
        plot_ly(z = state.area, text = state.name, locations = state.abb,
                type = 'choropleth', locationmode = 'USA-states') %>%
            layout(geo = g)
    })
    
    output$click <- renderPrint({
        d <- event_data("plotly_click")
        if (is.null(d)) "Click on a state to view event data" else d
    })
    
}

shinyApp(ui, server)