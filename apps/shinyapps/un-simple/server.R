shinyServer(function(input, output, session) {
    
    output$trendPlot <- renderPlotly({
        
        if (length(input$name) == 0) {
            print("Please select at least one country")
        } else {
            df_trend <- ideal[ideal$Name == input$name, ]
            ggplot(df_trend) +
                geom_line(aes(x = Year, y = Ideal.point, by = Name, color = Name)) +
                labs(x = "Year", y = "Ideology", title = "Ideal Points for Countries") +
                scale_colour_hue("clarity", l = 70, c = 150) + ggthemes::theme_few()
        }
        
    })
})
