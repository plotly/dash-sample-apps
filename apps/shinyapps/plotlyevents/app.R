library(plotly)
library(shiny)

ui <- fluidPage(
    radioButtons("plotType", "Plot Type:", choices = c("ggplotly", "plotly")),
    plotlyOutput("plot"),
    verbatimTextOutput("hover"),
    verbatimTextOutput("click"),
    verbatimTextOutput("brush"),
    verbatimTextOutput("zoom")
)

server <- function(input, output, session) {
    
    output$plot <- renderPlotly({
        # use the key aesthetic/argument to help uniquely identify selected observations
        key <- row.names(mtcars)
        if (identical(input$plotType, "ggplotly")) {
            p <- ggplot(mtcars, aes(x = mpg, y = wt, colour = factor(vs), key = key)) + 
                geom_point()
            ggplotly(p) %>% layout(dragmode = "select")
        } else {
            plot_ly(mtcars, x = ~mpg, y = ~wt, key = ~key) %>%
                layout(dragmode = "select")
        }
    })
    
    output$hover <- renderPrint({
        d <- event_data("plotly_hover")
        if (is.null(d)) "Hover events appear here (unhover to clear)" else d
    })
    
    output$click <- renderPrint({
        d <- event_data("plotly_click")
        if (is.null(d)) "Click events appear here (double-click to clear)" else d
    })
    
    output$brush <- renderPrint({
        d <- event_data("plotly_selected")
        if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
    })
    
    output$zoom <- renderPrint({
        d <- event_data("plotly_relayout")
        if (is.null(d)) "Relayout (i.e., zoom) events appear here" else d
    })
    
}

shinyApp(ui, server)
