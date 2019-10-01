appName <- Sys.getenv('DASH_APP_NAME')
if (appName != '') {
  pathPrefix <- sprintf('/%s/', appName)
  
  Sys.setenv(
    DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
    DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix
  )
  
  setwd(sprintf('/app/apps/%s', appName))
}

library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)

# Load & Prep Data ------------------------------

data(diamonds, package = 'ggplot2')

nms <- names(diamonds)

# Create Layout Variables ------------------------------

DropDownMenuOptions <- lapply(nms, function(x) {
  list(label = x, value = x)
})

SampleSlider <- dccSlider(
  id = 'sample-slider',
  min = 1,
  max = nrow(diamonds),
  marks = list(
    '1' = '1',
    '10000' = list('label' = '10K'),
    '20000' = list('label' = '20K'),
    '30000' = list('label' = '30K'),
    '40000' = list('label' = '40K'),
    '53940' = '53940'
  ),
  value = 1000
)

heightSlider <- dccSlider(
  id = 'height-slider',
  min = 100,
  max = 2000,
  marks = list(
    '100' = '100',
    '290' = '290',
    '480' = '480',
    '670' = '670',
    '860' = '860',
    '1050' = list('label' = '1,050'),
    '1240' = list('label' = '1,240'),
    '1430' = list('label' = '1,430'),
    '1620' = list('label' = '1,620'),
    '1810' = list('label' = '1,810'),
    '2000' = list('label' = '2,000')
  ),
  value = 1000
)

pageTitle <- htmlH2('Diamonds Explorer')

plotlyLogo <- htmlA(
  list(htmlImg(id = 'banner-image', src = 'assets/image.png')),
  className = 'logo',
  href = 'https://dashr.plot.ly')

xDropDown <- dccDropdown(
  id = 'x-dropdown',
  options = DropDownMenuOptions,
  value = 'carat',
  clearable = FALSE
)

yDropDown <- dccDropdown(
  id = 'y-dropdown',
  options = DropDownMenuOptions,
  value = 'price',
  clearable = FALSE
)

colorDropDown <- dccDropdown(
  id = 'color-dropdown',
  options = DropDownMenuOptions,
  value = 'clarity',
  clearable = FALSE
)

facetRowDropDown <- dccDropdown(
  id = 'facet-row-dropdown',
  options = DropDownMenuOptions,
  value = 'clarity',
  clearable = FALSE
)

facetColumnDropDown <- dccDropdown(
  id = 'facet-column-dropdown',
  options = DropDownMenuOptions,
  value = '.',
  clearable = FALSE
)


# App Start ------------------------------

# Initiate application
app = Dash$new()


# Create Layout ------------------------------

app$layout(
  htmlDiv(list(
    plotlyLogo,
    pageTitle
    )),
  htmlDiv(
    list(
      htmlB('Sample Size'),
      htmlDiv(SampleSlider, style = list('marginBottom' = 40)),
      htmlB('X'),
      xDropDown,
      htmlBr(),
      htmlB('Y'),
      yDropDown,
      htmlBr(),
      htmlB('Color'),
      colorDropDown,
      htmlBr(),
      htmlB('FacetRow'),
      facetRowDropDown,
      htmlBr(),
      htmlB('FacetColumn'),
      facetColumnDropDown,
      htmlBr(),
      htmlB('Height of plot (in pixels)'),
      heightSlider
    ),
    className = 'three columns'
  ),
  htmlDiv(list(dccGraph(id = 'facetgrid')),
          className = 'eight columns')
)


# Callbacks Start ------------------------------

app$callback(output = list(id = 'facetgrid', property = 'figure'),
             params = list(
               input(id = 'sample-slider', property = 'value'),
               input(id = 'height-slider', property = 'value'),
               input(id = 'x-dropdown', property = 'value'),
               input(id = 'y-dropdown', property = 'value'),
               input(id = 'color-dropdown', property = 'value'),
               input(id = 'facet-row-dropdown', property = 'value'),
               input(id = 'facet-column-dropdown', property = 'value')
             ),
             
  function(sampleSize, plotHeight, x, y, color, facet_row,facet_col) {

    #Adding Reactive data info, dataset is built in diamonds data
    dataset <- diamonds[sample(nrow(diamonds), sampleSize),]

    #FacetGrid ggplot2 object
    initial_plot <- ggplot(dataset, aes_string(x = x, y = y, color = color)) + geom_point()

    # Add it if least one facet column/row is specified
    facets <- paste(facet_row, '~', facet_col)


    if (facets != '. ~ .') { initial_plot <- initial_plot + facet_grid(facets)
      final_plot <- ggplotly(initial_plot, height = plotHeight) %>%
      layout(autosize = TRUE)
    }
    
    return(final_plot)
})

if (appName != '') {
  
  app$run_server(host = '0.0.0.0', port = Sys.getenv("PORT", 8050))
  
} else {
  
  app$run_server(showcase = TRUE)
}
