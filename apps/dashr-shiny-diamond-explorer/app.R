appName <- Sys.getenv("DASH_APP_NAME")
pathPrefix <- sprintf("/%s/", appName)

Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
           DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)

################# LOAD DATA, FUNCTIONS, & CREATE GLOBAL OBJECTS ########

data(diamonds, package = "ggplot2")
nms <- names(diamonds)

#rendering dropdown lists
lapply(nms, function(x) {
  list(label = x, value = x)
})

######### CREATE GLOBAL OBJECTS ####

DropDownMenuOptions <- list(
  list(label = "carat", value = "carat"),
  list(label = "cut", value = "cut"),
  list(label = "color", value = "color"),
  list(label = "clarity", value = "clarity"),
  list(label = "depth", value = "depth"),
  list(label = "table", value = "table"),
  list(label = "price", value = "price"),
  list(label = "x", value = "x"),
  list(label = "y", value = "y"),
  list(label = "z", value = "z")
)


FacetOptions <- list(
  list(label = "None", value = "None"),
  list(label = "cut", value = "cut-labels"),
  list(label = "color", value = "color-labels"),
  list(label = "clarity", value = "clarity")
)

SampleSlider <- dccSlider(id = "sample-slider",
                          min = 1,
                          max = 53940,
                          marks = list(
                            "1" = "1", 
                            "10K" = "10K",
                            "20K" = "20K",
                            "30K" = "30K",
                            "40K" = "40K", 
                            "53940" = "53940")
                          ,
                          value = 1000
)

heightSlider <- dccSlider(id = "height-slider",
                          min = 100,
                          max = 2000,
                          marks = list(
                            "100" = "100",
                            "290" = "290",
                            "480" = "480",
                            "670" = "670",
                            "860" = "860",
                            "1,050" = "1,050",
                            "1,240" = "1,240",
                            "1,430" = "1,430",
                            "1,620" = "1,620",
                            "1,810" = "1,810",
                            "2,000" = "2,000"),
                          value = 1000
)


##### CREATE LAYOUT VARIABLES #######

pageTitle <- htmlH2("Diamonds Explorer ")

plotlyLogo <-
  htmlA(list(htmlImg(id = "banner-image", src = "assets/image.png")), className = "logo",
        href = "https://dashr-docs.herokuapp.com/")

xDropDown <- dccDropdown(id = "x-dropdown",
                         options = DropDownMenuOptions,
                         value = "carat",
                         clearable = FALSE)

yDropDown <- dccDropdown(id = "y-dropdown",
                         options = DropDownMenuOptions,
                         value = "cut",
                         clearable = FALSE)

FacetDropDown <- dccDropdown(id="facet-row-dropdown",
                             options = FacetOptions,
                             value = "None",
                             clearable = FALSE)


###### APP START #######

app = Dash$new()

########### CREATE LAYOUT ##############

app$layout(
  htmlDiv(
    list(
      htmlDiv(list(
        plotlyLogo
      ))
    )
  )
)


################ CALLBACKS #####################

app$callback(output = list(id = "histogram", property = "figure"),
             params = list(input(id = "movies-slider", property = "value")),
             
             function(input, output) {
               #add reactive data information. Dataset = built in diamonds data
               dataset <- reactive({
                 diamonds[sample(nrow(diamonds), sampleSize)]
                          
             # build graph with ggplot syntax
             p <- ggplot(dataset(), 
                         aes_string(x = x, 
                                    y = y, 
                                    color = color)) + geom_point()
               )
             # if at least one facet column/row is specified, add it
             facets <- paste(input$facet_row, '~', input$facet_col)
             if (facets != '. ~ .') p <- p + facet_grid(facets)
             
             ggplotly(p) %>% 
               layout(height = input$plotHeight, autosize=TRUE)
             
             }
             return(p)
             })
)
########CONDITIONAL STATEMENT FOR APP RUNNING ON CLOUD SERVER & LOCAL

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server(showcase = TRUE)
}
