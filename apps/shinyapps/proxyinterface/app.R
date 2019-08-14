library(shiny)
library(plotly)

ui <- fluidPage(
    div(
        h3("PLOTLYPROXY", align = 'center'), 
        style = "border-bottom: solid; border-width: thin;"
    ),
    br(),
    fluidRow(
        column(3,
               h4("addTraces", style = "text-decoration: underline;"),
               fluidRow(
                   column(6, textInput("xvalue", NULL, "1,2,3,4")),
                   column(6, textInput("yvalue", NULL, "3,3,3,3"))),
               actionButton("add", "Add", width = '210px')
        ),
        column(3,
               h4("restyle", style = "text-decoration: underline;"),
               textInput("traceNo", NULL, "Trace #"),
               selectInput("marker", NULL, 
                           schema()$traces$scatter$attributes$marker$symbol$values, 
                           selected = "circle")
        ),
        column(3,
               h4("relayout", style = "text-decoration: underline;"),
               textInput("text", NULL, "Enter New Title Here"),
               selectInput("color", NULL, colors(), selected = "black")
        ),
        column(3,
               h4("deleteTraces", style = "text-decoration: underline;"),
               textInput("traceNo2", NULL, "Trace #"),
               actionButton("delete", "Delete", width = '210px')
        )
    ),
    br(),
    div(
        plotlyOutput("plot"),
        style = "border: solid black 1px"
    )
)

server <- function(input, output, session) {
    
    output$plot <- renderPlotly({
        plot_ly(
            type = 'scatter',
            mode = 'markers',
            x = c(1,2,3,4), 
            y = c(2,4,2,4)
        ) %>%
            layout(
                title = 'Original Title',
                showlegend = T,
                margin = list(t=40)
            )
    })
    
    # plotly.addTraces
    observeEvent(input$add, {
        plotlyProxy("plot", session) %>%
            plotlyProxyInvoke("addTraces", list(x = as.integer(unlist(strsplit(input$xvalue,","))), 
                                                y = as.integer(unlist(strsplit(input$yvalue, ","))),
                                                type = 'scatter',
                                                mode = 'markers'))
    })
    
    # plotly.restyle
    observeEvent(input$marker, {
        plotlyProxy("plot", session) %>%
            plotlyProxyInvoke("restyle", list(marker = list(symbol = input$marker)), list(input$traceNo))
    })
    
    # plotly.relayout
    observeEvent(input$text, {
        plotlyProxy("plot", session) %>%
            plotlyProxyInvoke("relayout", list(title = input$text))
    })
    
    observeEvent(input$color, {
        plotlyProxy("plot", session) %>%
            plotlyProxyInvoke("relayout", list(plot_bgcolor = input$color, paper_bgcolor = input$color))
    })
    
    # plotly.deleteTraces
    observeEvent(input$delete, {
        plotlyProxy("plot", session) %>%
            plotlyProxyInvoke("deleteTraces", list(as.integer(input$traceNo2)))
    })
    
}

shinyApp(ui, server)