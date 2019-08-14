# server.R definition
server <- function(input, output){
    
    # Observes the second feature input for a change
    observeEvent(input$featureInput2,{
        
        # Create a convenience data.frame which can be used for charting
        plot.df <- data.frame(BreastCancer[,input$featureInput1],
                              BreastCancer[,input$featureInput2],
                              Class = BreastCancer$Class)
        
        # Add column names
        colnames(plot.df) <- c("x", "y", "Class")
        
        # Do a plotly contour plot to visualize the two featres with
        # the number of malignant cases as size
        # Note the use of 'source' argument
        output$Plot1 <- renderPlotly({
            plot_ly(plot.df, x = ~x, y = ~y, mode = "markers", type = "scatter", color = ~Class, source = "subset",
                    marker = list(size = 30)) %>%
                layout(title = paste(input$featureInput1, "vs ", input$featureInput2),
                       xaxis = list(title = input$featureInput1),
                       yaxis = list(title = input$featureInput2),
                       dragmode =  "select",
                       plot_bgcolor = "6A446F")
        })
        
        # Create a contour plot of the number of malignant cases
        output$Plot2 <- renderPlotly({
            
            plot.df %>%
                group_by(x, y, Class) %>%
                summarize(Count = n()) %>%
                filter(Class == "malignant") %>%
                plot_ly(x = ~x, y = ~y, z = ~Count, type = "contour") %>%
                layout(title = "Contour map of number of malignant cases",
                       xaxis = list(title = input$featureInput1),
                       yaxis = list(title = input$featureInput2))
            
        })
        
        # Assign to parent environment
        plot.df <<- plot.df
    })
    
    # Coupled event 1
    output$Plot3 <- renderPlotly({
        
        # Get subset based on selection
        event.data <- event_data("plotly_selected", source = "subset")
        
        # If NULL dont do anything
        if(is.null(event.data) == T) return(NULL)
        
        # Get number of malignant and benign cases from selection
        malig.class <- subset(plot.df, Class == "malignant")[subset(event.data, curveNumber == 0)$pointNumber + 1,]
        benign.class <- subset(plot.df, Class == "benign")[subset(event.data, curveNumber == 1)$pointNumber + 1,]
        
        # Combine
        plot.subset <- rbind(malig.class, benign.class)
        
        # Summarize
        plot.summ <- plot.subset %>%
            group_by(x, y, Class) %>%
            summarize(Count = n())
        
        # Assign to parent frame
        plot.summ <<- plot.summ
        
        # Plot
        plot_ly(plot.summ, x = ~Class, y = ~Count, type = "bar", source = "select", color = ~Class) %>%
            layout(title = "No. of Malignant and Benign cases <br> in Selection",
                   plot_bgcolor = "6A446F",
                   yaxis = list(domain = c(0, 0.9)))
    })
    
    # Coupled event 2
    output$Plot4 <- renderPlotly({
        
        # Get subset based on selection
        event.data <- event_data("plotly_click", source = "select")
        
        # If NULL dont do anything
        if(is.null(event.data) == T) return(NULL)
        
        # If Malignant
        if(event.data[3] == "malignant"){
            
            tab <- subset(plot.summ, Class == "malignant")
            
            p1 <- plot_ly(tab, x = x, y = Count, type = "box", showlegend = F) %>%
                layout(yaxis = list(title = "Count"),
                       xaxis = list(title = input$featureInput1))
            
            p2 <- plot_ly(tab, x = y, y = Count, type = "box", showlegend = F) %>%
                layout(title = "Box plot for Malignant cases",
                       yaxis = list(title = "Count"),
                       xaxis = list(title = input$featureInput2))
            
            subplot(p1, p2)
        }else{
            tab <- subset(plot.summ, Class == "benign")
            
            p1 <- plot_ly(tab, x = ~x, y = ~Count, type = "box", showlegend = F) %>%
                layout(yaxis = list(title = "Count"),
                       xaxis = list(title = input$featureInput1))
            
            p2 <- plot_ly(tab, x = ~y, y = ~Count, type = "box", showlegend = F) %>%
                layout(title = "Box plot for Benign cases",
                       yaxis = list(title = "Count"),
                       xaxis = list(title = input$featureInput2))
            
            subplot(p1, p2)
        }
        
    })
    
}