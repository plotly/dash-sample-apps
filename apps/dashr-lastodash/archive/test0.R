library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)

app <- Dash$new()

df <- read.csv('dashr/chapters/getting-started/examples/gdp-life-exp-2007.csv', header = TRUE, sep = ",")
continents <- c(unique(df$continent))
data_gdp_life_exp_2007 <- c()
continents <- unique(df$continent)

data_gdp_life_exp_2007 <- list()
for (i in 1:length(continents)){
  data_gdp_life_exp_2007[[i]] <- list(
    x = df$gdp_per_capita[df$continent==continents[i]],
    y = df$life_expectancy[df$continent==continents[i]],
    opacity=0.7,
    text = continents[i],
    mode = 'markers',
    name = continents[i],
    marker = list(size = 15,
                  line = list(width = 0.5, color = 'white'))
  )}

app$layout(htmlDiv(list(
  dccGraph(
    id = 'life-exp-vs-gdp',
    figure = list(
      data =  data_gdp_life_exp_2007,
      layout = list(
        xaxis = list('type' = 'log', 'title' = 'GDP Per Capita'),
        yaxis = list('title' = 'Life Expectancy'),
        margin = list('l' = 40, 'b' = 40, 't' = 10, 'r' = 10),
        legend = list('x' = 0, 'y' = 1),
        hovermode = 'closest'
      )
    )
  )
)))

#app$run_heroku()