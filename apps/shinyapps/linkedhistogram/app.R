library(plotly)
library(shiny)

# user interface
ui <- fluidPage(
    titlePanel("Linked highlighting with plotly and shiny"),
    mainPanel(
        htmltools::div(style = "display:inline-block", plotlyOutput("x", width = 400, height = 250)),
        wellPanel(
            style = "display:inline-block; vertical-align:top;", 
            sliderInput("xbins", "Number of x bins", 
                        min = 1, max = 50, value = 20, width = 250),
            sliderInput("ybins", "Number of y bins", 
                        min = 1, max = 50, value = 20, width = 250)
        ),
        br(),
        htmltools::div(style = "display:inline-block", plotlyOutput("xy", width = 400, height = 400)),
        htmltools::div(style = "display:inline-block", plotlyOutput("y", width = 250, height = 400))
    )
)

# marker objects
m <- list(color = toRGB("black"))
m2 <- list(color = toRGB("black", 0.2))

server <- function(input, output, session) {
    
    # convenience function for computing xbin/ybin object given a number of bins
    compute_bins <- function(x, n) {
        list(
            start = min(x),
            end = max(x),
            size = (max(x) - min(x)) / n
        )
    }
    
    # the 'x' histogram
    output$x <- renderPlotly({
        x <- cars$speed
        xbins <- compute_bins(x, input$xbins)
        p <- plot_ly(x = x, type = "histogram", autobinx = F, 
                     xbins = xbins, marker = m2)
        # obtain plotlyjs selection
        s <- event_data("plotly_selected")
        # if points are selected, subset the data, and highlight
        if (length(s$x) > 0) {
            p <- add_trace(p, x = s$x, type = "histogram", autobinx = F, 
                           xbins = xbins, marker = m)
        }
        p %>%
            config(displayModeBar = F, showLink = F) %>%
            layout(showlegend = F, barmode = "overlay", yaxis = list(title = "count"),
                   xaxis = list(title = "", showticklabels = F))
    })
    
    # basically the same as 'x' histogram
    output$y <- renderPlotly({
        y <- cars$dist
        ybins <- compute_bins(y, input$ybins)
        p <- plot_ly(y = y, type = "histogram", autobiny = F, 
                     ybins = ybins, marker = m2)
        s <- event_data("plotly_selected")
        if (length(s$y) > 0) {
            p <- add_trace(p, y = s$y, type = "histogram", autobiny = F, 
                           ybins = ybins, marker = m)
        }
        p %>%
            config(displayModeBar = F, showLink = F) %>%
            layout(showlegend = F, barmode = "overlay", xaxis = list(title = "count"),
                   yaxis = list(title = "", showticklabels = F))
    })
    
    output$xy <- renderPlotly({
        cars %>% 
            plot_ly(x = ~speed, y = ~dist, 
                    mode = "markers", marker = m) %>%
            layout(dragmode = "select")
    })
    
}

shinyApp(ui, server)