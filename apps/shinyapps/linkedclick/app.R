library(plotly)
library(shiny)

# compute a correlation matrix
correlation <- round(cor(mtcars), 3)
nms <- names(mtcars)

ui <- fluidPage(
    mainPanel(
        plotlyOutput("heat"),
        plotlyOutput("scatterplot")
    ),
    verbatimTextOutput("selection")
)

server <- function(input, output, session) {
    output$heat <- renderPlotly({
        plot_ly(x = nms, y = nms, z = correlation, 
                key = correlation, type = "heatmap", source = "heatplot") %>%
            layout(xaxis = list(title = ""), 
                   yaxis = list(title = ""))
    })
    
    output$selection <- renderPrint({
        s <- event_data("plotly_click")
        if (length(s) == 0) {
            "Click on a cell in the heatmap to display a scatterplot"
        } else {
            cat("You selected: \n\n")
            as.list(s)
        }
    })
    
    output$scatterplot <- renderPlotly({
        s <- event_data("plotly_click", source = "heatplot")
        if (length(s)) {
            vars <- c(s[["x"]], s[["y"]])
            d <- setNames(mtcars[vars], c("x", "y"))
            yhat <- fitted(lm(y ~ x, data = d))
            plot_ly(d, x = ~x) %>%
                add_markers(y = ~y) %>%
                add_lines(y = ~yhat) %>%
                layout(xaxis = list(title = s[["x"]]), 
                       yaxis = list(title = s[["y"]]), 
                       showlegend = FALSE)
        } else {
            plotly_empty()
        }
    })
    
}

shinyApp(ui, server)