library(shiny)
library(plotly)

shinyUI(fluidPage(
    
    # Application title
    titlePanel("Ideal Points"),
    
    sidebarPanel(
        h3("Ideal Points Estimation"),
        # Select Justices name here
        selectizeInput("name",
                       label = "Country Name(s) of Interest",
                       choices = unique(ideal$Name),
                       multiple = T,
                       options = list(maxItems = 5, placeholder = 'Select a name'),
                       selected = "United States of America"),
        # Term plot
        plotOutput("termPlot", height = 200),
        helpText("Data: Bailey, Michael, Anton  Strezhnev and Erik Voeten. Forthcoming.  'Estimating Dynamic State Preferences from United Nations Voting Data.' Journal of Conflict Resolution. ")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
        plotlyOutput("trendPlot")
    )
)
)