library(dashR)
library(plotly)
library(dashHtmlComponents)
library(dashCoreComponents)

app <- Dash$new()

app$layout(
  htmlDiv(list(
    htmlDiv(list(
      htmlH3("Plotly.js importing test using Edgar Anderson's iris data")
    ),
    className = 'row'),
    htmlDiv(list(
      dccDropdown(id = 'species-selector',
                  value = list('setosa', 'versicolor', 'virginica'),
                  options = list(list(value = "setosa",
                                      label = "Setosa"),
                                 list(value = "versicolor",
                                      label = "Versicolor"),
                                 list(value = "virginica",
                                      label = "Virginica")),
                  multi = TRUE)),
      className = 'row'),
    htmlDiv(list(
      htmlH4('Sepal width (cm) vs. sepal length (cm)', className = 'six columns'),
      htmlH4('Petal length (cm) vs. sepal length (cm)', className = 'six columns')
    ),
    className = 'row'),
    htmlDiv(list(
      htmlDiv(
        id = 'div-sepal-sepal',
        className = 'six columns'),
      htmlDiv(
        id = 'div-petal-sepal',
        className = 'six columns')), 
      className='row')
  )
  )
)

app$callback(output=list(id='div-petal-sepal', property='children'),
             params=list(
               input(id='species-selector', property='value')),
             function(species)
             {
               if (!(length(species))) {
                 subdata <- iris
               } else {
                 subdata <- iris[which(iris$Species == species),]                 
               }
               
               dccGraph(id = 'graph-petal-sepal', figure = plot_ly(data = subdata, x = ~Sepal.Length, y = ~Petal.Length))
             }
)

app$callback(output=list(id='div-sepal-sepal', property='children'),
             params=list(
               input(id='species-selector', property='value')),
             function(species)
             {
               if (!(length(species))) {
                 subdata <- iris
               } else {
                 subdata <- iris[which(iris$Species == species),]                 
               }
               
               dccGraph(id = 'graph-sepal-sepal', figure = plot_ly(data = subdata, x = ~Sepal.Length, y = ~Sepal.Width))
             }
)

app$run_server()
