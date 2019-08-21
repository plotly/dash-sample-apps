appName <- Sys.getenv('DASH_APP_NAME')
if (appName != '') {
  pathPrefix <- sprintf('/%s/', appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf('/app/apps/%s', appName))
}

library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)

################# LOAD DATA, FUNCTIONS, & CREATE GLOBAL OBJECTS ########

data(diamonds, package = 'ggplot2')
nms <- names(diamonds)

#Adding Reactive data info, dataset is built in diamonds data

#rendering dropdown lists
lapply(nms, function(x) {
  list(label = x, value = x)
})

######### CREATE GLOBAL OBJECTS ####

DropDownMenuOptions <- list(
  list(label = 'None', value = 'None'),
  list(label = 'carat', value = 'carat'),
  list(label = 'cut', value = 'cut'),
  list(label = 'color', value = 'color'),
  list(label = 'clarity', value = 'clarity'),
  list(label = 'depth', value = 'depth'),
  list(label = 'table', value = 'table'),
  list(label = 'price', value = 'price'),
  list(label = 'x', value = 'x'),
  list(label = 'y', value = 'y'),
  list(label = 'z', value = 'z')
)

SampleSlider <- dccSlider(id = 'sample-slider',
                          min = 1,
                          max = 53940,
                          marks = list(
                            '1' = '1', 
                            '10000' = list('label' = '10K'),
                            '20000' = list('label' = '20K'),
                            '30000' = list('label' = '30K'),
                            '40000' = list('label' = '40K'),
                            '53940' = '53940'),
                          value = 1000
)
heightSlider <- dccSlider(id = 'height-slider',
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
                            '2000' = list('label' = '2,000')),
                          value = 1000
)
##### CREATE LAYOUT VARIABLES #######

pageTitle <- htmlH2('Diamonds Explorer')
plotlyLogo <-
  htmlA(list(htmlImg(id = 'banner-image', src = 'assets/image.png')), className = 'logo',
        href = 'https://dashr.plot.ly')
xDropDown <- dccDropdown(id = 'x-dropdown',
                         options = DropDownMenuOptions,
                         value = 'carat',
                         clearable = FALSE)
yDropDown <- dccDropdown(id = 'y-dropdown',
                         options = DropDownMenuOptions,
                         value = 'price',
                         clearable = FALSE)
ColorDropDown <- dccDropdown(id='color-dropdown',
                             options = DropDownMenuOptions,
                             value = 'clarity',
                             clearable = FALSE)

FacetRowDropDown <- dccDropdown(id='facet-row-dropdown',
                              options = DropDownMenuOptions,
                              value = 'clarity',
                              clearable = FALSE)
FacetColumnDropDown <- dccDropdown(id ='facet-column-dropdown',
                               options = DropDownMenuOptions,
                               value = 'None',
                               clearable = FALSE)
###### APP START #######
app = Dash$new()

########### CREATE LAYOUT ##############
app$layout(
  htmlDiv(list(
      plotlyLogo,
      pageTitle
      ), className = 'twelve columns'
  ),
  htmlDiv(list(
    htmlb('Sample Size'),
    htmlDiv(SampleSlider, style = list('marginBottom' = 40)
    ),
    htmlB('X'),
    xDropDown,
    htmlBr(),
    htmlB('Y'),
    yDropDown,
    htmlBr(),
    htmlB('Color'),
    ColorDropDown,
    htmlBr(),
    htmlB('FacetRow'),
    FacetRowDropDown,
    htmlBr(),
    htmlB('FacetColumn'),
    FacetColumnDropDown,
    htmlBr(),
    htmlB('Height of plot (in pixels)'),
    heightSlider
  ), className = 'three columns'
  ),
  htmlDiv(list(
    dccGraph(id = 'facetgrid')
  ),
  className = 'eight columns'
  )
)
################ CALLBACKS #####################


########CONDITIONAL STATEMENT FOR APP RUNNING ON CLOUD SERVER & LOCAL

if (appName != '') {
  app$run_server(host = '0.0.0.0', port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server(showcase = TRUE)
}
