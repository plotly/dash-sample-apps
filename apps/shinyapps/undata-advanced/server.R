#server script for United Nations Advanced Example
shinyServer(function(input, output, session) {
    
    output$trendPlot <- renderPlotly({
        if (length(input$name) == 0) {
            print("Please select at least one country")
        } else {
            df_trend <- ideal[ideal$Name == input$name, ]
            
            # Graph title
            if (length(input$name) > 2) {
                j_names_comma <- paste(input$name[-length(input$name)], collapse = ', ')
                j_names <- paste0(j_names_comma, ", and ", input$name[length(input$name)])
            } else {
                j_names <- paste(input$name, collapse = ' and ')
            }
            
            graph_title  <- paste("Ideal Points for ", j_names, sep="")
            
            ggideal_point <- ggplot(df_trend) +
                geom_line(aes(x = Year, y = Ideal.point, by = Name, color = Name)) +
                labs(x = "Year", y = "Ideology", title = graph_title) +
                scale_colour_hue("clarity", l = 70, c = 150) + ggthemes::theme_few() +
                theme(legend.direction = "horizontal", legend.position = "bottom")
            
            # Convert ggplot object to plotly
            gg <- plotly_build(ggideal_point)
            
            # Use Plotly syntax to further edit the plot:
            gg$layout$annotations <- NULL # Remove the existing annotations (the legend label)
            gg$layout$annotations <- list()
            
            # Add colored text annotations next to the end of each line
            # More about plotly annotations: https://plot.ly/r/reference/#annotation
            # Each key that we update is documented in that link above.
            for (i in 1:length(gg$data)) { # data is a list of the lines in the graph
                gg$layout$annotations[[i]] <- list(
                    text = gg$data[[i]]$name,  # The text label of the annotation, e.g. "Canada"
                    font = list(color = gg$data[[i]]$line$color), # Match the font color to the line color
                    showarrow = FALSE, # Don't show the annotation arrow
                    y = gg$data[[i]]$y[[length(gg$data[[i]]$y)]], # set the y position of the annotation to the last point of the line
                    yref = "y1", # the "y" coordinates above are with respect to the yaxis
                    x = 1, # set the x position of the graph to the right hand side of the graph
                    xref = "paper", # the x coordinates are with respect to the "paper", where 1 means the right hand side of the graph and 0 means the left hand side
                    xanchor = "left" # position the x coordinate with respect to the left of the text
                );
            }
            
            gg$layout$showlegend <- FALSE # remove the legend
            gg$layout$margin$r <- 170 # increase the size of the right margin to accommodate more room for the annotation labels
            gg
            
        }
    })
    
    output$termPlot <- renderPlot({
        df_term <- ideal  %>% filter(Name %in% input$name) %>%
            group_by(Name) %>% summarise(terms = n())
        
        trans_theme <- theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.background = element_rect(fill = NA),
            plot.background = element_rect(fill = NA)
        )
        
        ggplot(df_term, aes(x = reorder(Name, terms), y = terms))+
            geom_bar(stat = "identity", fill = "#2980b9") + coord_flip() +
            theme_bw() + trans_theme + labs(y = "Terms (in years)", x = "")
        
    }, bg="transparent")
})
