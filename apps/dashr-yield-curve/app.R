library(plotly)
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(quantmod)
library(magrittr)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(glue)
library(jsonlite)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

#Helper functions

genMarksEmpty <- function(n){
  l <- list(glue(''))
  names(l) <- 'label'
  m <- list(l)
  names(m) <- glue('{n}')
  return(m)
}

safemerge <- function(v1, v2){
  #for columns for vectors use mapply(safemerge, v1, v2)
  if(is.na(v1)||is.null(v1)||rlang::is_empty(v1)){
    return(v2)
  } else{
    return(v1)
  }
  
}

safelabel <- function(v1){
  
  if(is.na(v1)|| is.null(v1)||rlang::is_empty(v1)){
    return("3-month")
  } else{
    return('1-month')
  }
  
}

last_back <- 0
last_next <- 0

DT_ <- fread('data/yield_curve.csv') %>% drop_na(., x,y)
dates <- DT_$y
Dates <- as.Date(dates, format='%m/%d/%Y')
labels <-(DT_$x)[1:11]
DT <- DT_[, 3:13]
colnames(DT) <- labels
rownames(DT) <- dates
DM <- DT %>% data.matrix()

UPS <- list('0'= list(x=0, y=0, z=1), '1'=list(x=0, y=0, z=1), '2'=list(x=0, y=0, z=1), '3'=list(x=0, y=0, z=1), '4'=list(x=0, y=0, z=1), '5'=list(x=0, y=0, z=1))

CENTERS <- list('0'= list(x=0.3, y=0.8, z=-0.5), 
              '1'=list(x=0, y=0, z=-0.37), 
              '2'=list(x=0, y=1.1, z=-1.3), 
              '3'=list(x=0, y=-0.7, z=0), 
              '4'=list(x=0, y=-0.2, z=0), 
              '5'=list(x=-0.11, y=-0.5, z=0))

EYES <-  list('0'= list(x=2.7, y=2.7, z=0.3), 
            '1'=list(x=0.01, y=3.8, z=-0.37), 
            '2'=list(x=1.3, y=3, z=0), 
            '3'=list(x=2.6, y=-1.6, z=0), 
            '4'=list(x=3, y=-0.2, z=0), 
            '5'=list(x=-0.1, y=-0.5, z=2.66))

TEXTS <- list(
  list(
    htmlDiv(dccMarkdown(className='six columns', children=" #### Yield curve 101 "), style=list(width='120%')), #htmlBr(),
    htmlDiv(dccMarkdown(className='six columns', children=" The yield curve shows how much it costs the federal government to borrow
                 money for a given amount of time, revealing the relationship between long-
                   and short-term interest rates.  It is, inherently, a forecast for what the economy holds in the future -
                 how much inflation there will be, for example, and how healthy growth will
                 be over the years ahead - all embodied in the price of money today,
                 tomorrow and many years from now. "), style=list(width='120%'))), 
  
  list(htmlDiv(dccMarkdown(className='six columns', children= "#### Where we stand"), style=list(width='120%')), #htmlBr(),
     htmlDiv( dccMarkdown(className='six columns', children="On Wednesday, both short-term and long-term rates were lower than they have
    been for most of history - a reflection of the continuing hangover
    from the financial crisis. The yield curve is fairly flat, which is a sign that investors expect
    mediocre growth in the years ahead."), style=list(width='120%'))), 
  list(htmlDiv(dccMarkdown(className='six columns', children=' #### Deep in the valley'),style=list(width='120%')), #htmlBr(),
     htmlDiv(dccMarkdown(className='six columns', children='In response to the last recession, the Federal Reserve has kept short-term
    rates very low - near zero - since 2008. (Lower interest rates stimulate
    the economy, by making it cheaper for people to borrow money, but also
    spark inflation.) Now, the Fed is getting ready to raise rates again, possibly as early as
    June.'), style=list('width'='120%'))), 
  list(htmlDiv(dccMarkdown(className='six columns', children=' #### Last time, a puzzle'), style=list(width='120%')), #htmlBr(),
     htmlDiv(dccMarkdown(className='six columns', children='The last time the Fed started raising rates was in 2004. From 2004 to 2006,
    short-term rates rose steadily. But long-term rates did not rise very much. The Federal Reserve chairman called this phenomenon a "conundrum," and it
    raised questions about the ability of the Fed to guide the economy.
    Part of the reason long-term rates failed to rise was because of strong
    foreign demand.') , style=list('width'='120%'))), 
  list(htmlDiv(dccMarkdown(className='six columns', children='#### Long-term rates are low now, too'), style=list(width='120%')), #htmlBr(),
     htmlDiv(dccMarkdown(className='six columns', children= 'Foreign buyers have helped keep long-term rates low recently, too - as have
    new rules encouraging banks to hold government debt and expectations that
    economic growth could be weak for a long time. The 10-year Treasury yield was as low as it has ever been in July 2012 and
    has risen only modestly since.Some economists refer to the economic pessimism as "the new normalist"'), style=list('width'='120%'))), 
  list(htmlDiv(dccMarkdown(className='six columns', children=' #### Long-term rates are low now, too'),style=list(width='120%') ), #htmlBr(),
     htmlDiv(dccMarkdown(className='six columns', children= 'Here is the same chart viewed from above.'), style=list(width='120%')))
)

##INITALIZE THE GRAPH
graph <- DM %>% .[,dim(.)[2]:1] %>% plot_ly(x=rev(labels), y=as.Date(dates, format='%m/%d/%Y'),z=., type='surface', 
                        hoverinfo='x+y+z', lighting=list('ambient'=0.95, 'diffuse'=0.99, 'fresnel'=0.01, 'roughness'=0.01, 'specular'=0.01 ),
                        colorscale=list(list(0, 'rgb(230,245,254)'), list(0.4,'rgb(123,171,203)'), list(0.8, 'rgb(40,119,174)'), list(1, 'rgb(37,61,81)')), zmax=9.18, zmin=0, scene='scene'
)
#graph


#Setup the app
app <- Dash$new()

app$layout(
  htmlDiv(list(
    htmlDiv(list(
              dccMarkdown(
              className='eight columns offset-by-two'
              ,
              children=" ### A View of a Chart That Predicts The Economic Future: The Yield Curve"
              )
              ,
              dccMarkdown(className='eight columns offset-by-two'
              , 
              children='This interactive report is a rendition of a
      [New York Times original](https://www.nytimes.com/interactive/2015/03/19/upshot/3d-yield-curve-economic-growth.html).')
              )
            ,
            className='row'
            ,
            style=list('text-align'='center', 'margin-bottom'= '15px')

            )
    ,

    htmlDiv(list(
      htmlDiv(list(
          dccSlider(min=0, max=5, value=0, marks=lapply(0:5, genMarksEmpty), id='slider')
          ), className='row',style=list('margin-bottom'= '10px'))
      ,
      htmlDiv(list(
          htmlDiv(list(
              htmlButton('Back', id='back', style=list('display'='inline-block'),n_clicks=0)#, className='learn-more-button')
              ,
              htmlButton('Next', id='next', style=list('display'='inline-block'),n_clicks=0)#,  className='learn-more-button')
              ),className='two columns offset-by-two')
          ,
          htmlDiv(id='text', className='six columns', children='')

          ), className='row', style=list('margin-bottom'= '10px'))
      ,
      dccGraph(
        id='graph',
        style=list('height'= '60vh'),
        figure=graph
      )
      ,
      htmlDiv(id='bucket', style=list(display='none'))
      ,
      htmlDiv(id='dummy')
    ), id='page')
  ))
)

## BUTTON CONTROLS

####HELPER function

app$callback(
  output=list(id='bucket', property='children')
  ,
  params = list(
    input(id='back', property='n_clicks'), 
    input(id='next', property='n_clicks'))
  ,
  function(beta, nu)
  {
    if(is.na(beta)||is.null(beta)||rlang::is_empty(beta)){
      
      b = 0
      
    } else
    {
      b <- as.numeric(beta[[1]])    
    }
    
    if(is.na(nu)||is.null(nu)||rlang::is_empty(nu))
    {
      n = 1
    } else
    {
      n <- as.numeric(nu[[1]])    
    }
    
    V <- (n-b) %% 6
    
    return(toJSON(V, null='null'))
  }
)


app$callback(
  output = list(id='slider', property='value')
  ,
  params = list(input(id='bucket', property='children'), state(id='slider', property='value'))
  ,
  function(V, S){
    value <- (fromJSON(V))
    if(is.na(value)||is.null(value)||length(value)==0||rlang::is_empty(value)){return(0)} else {
    value <- value[[1]] %>% as.numeric
      if(value <0){
        return(0)
      } else if (value >5){
        return(5)
      } else{
        return(value)
      }
    }
  }
)

app$callback(
  output=list(id='dummy', property='children')
  ,
  params=list(input(id='bucket', property='children'))
  ,
  function(n){
    return(n)
  }
)

app$callback(
  output=list(id='text', property='children')
  ,
  params=list(input(id='slider', property='value'))
  ,
  function(value){
    if(is.null(value)||is.na(value)){value=0} else{value}
    return(TEXTS[[(value+1)]])
  }
) 

app$callback(
  output=list(id='graph', property='figure')
  ,
  params=list(input(id='slider', property='value'))
  ,
  function(value)
  {
    if(is.null(value)||is.na(value)||rlang::is_empty(value)){
      v <- 0
    } else{
      v <- value
    }

    if(v %in% c(0,2,3)){

      z_snd <- mapply(safemerge,DT_[['z[0]']], DT_[['z[1]']])
      x_snd <- sapply(DT_[['z[0]']], safelabel)
      y_snd <- Dates 
      opacity <- 0.7

    } else if(v==1){

      x_snd <- rev(labels) 
      y_snd <- rep(Dates[length(Dates)], 11)
      z_snd <- tail(DM %>% .[,dim(.)[2]:1] , 1)  %>% as.numeric
      opacity <- 0.7

    } else if(v==4){
      z_snd <- DT_[['z[8]']]
      x_snd <- c("10-year")
      y_snd <- Dates
      opacity=0.25
      
    }
    
    if(v %in% 0:4){
      p <- plotly::add_trace(
        graph, 
        type='scatter3d', 
        mode='lines', 
        x=x_snd, 
        y=y_snd, 
        z=z_snd, 
        hoverinfo='x+y+z', 
        line=list(color='#444444'))
     
      
    } else{
 
      p <- plotly::plot_ly(type='contour', 
                           x=Dates, 
                           y=rev(labels), 
                           z=t(DM),  
                           hoverinfo='x+y+z', 
                           lighting=list('ambient'=0.95, 
                                         'diffuse'=0.99, 
                                         'fresnel'=0.01, 
                                         'roughness'=0.01, 
                                         'specular'=0.01 ), 
                           colorscale=list(list(0, 'rgb(230,245,254)'), list(0.4,'rgb(123,171,203)'), list(0.8, 'rgb(40,119,174)'), list(1, 'rgb(37,61,81)')), 
                           zmax=9.18, 
                           zmin=0, 
                           showscale=F,  
                           line=list(smoothing=1, color='rgba(40,40,40,0.15)'), 
                           contours=list(coloring='heatmap'))
    
    }
    
    p <- p %>% plotly::layout(.,autosize=T, 
                              font=list(size=12, color='#CCCCCC'), 
                              margin=list(t=5, l=5, b=5, r=5), 
                              showlegend=F, hovermode='closest',
                              scene=list(
                                aspectmode='manual', 
                                aspectratio=list(x=2,y=5,z=1.5), 
                                camera=list(up=UPS[[v+1]], 
                                center=CENTERS[[v+1]], 
                                eye=EYES[[v+1]]),   
                                annotations=list(list(
                                showarrow=F,
                                y="2015-03-18",
                                x="1-month",
                                z=0.046,
                                text="Point 1",
                                xanchor="left",
                                xshift=10,
                                opacity=0.7
                              ), list(
                                y="2015-03-18",
                                x="3-month",
                                z=0.048,
                                text="Point 2",
                                textangle=0,
                                ax=0,
                                ay=-75,
                                font=list(
                                  color="black",
                                  size=12
                                ),
                                arrowcolor="black",
                                arrowsize=3,
                                arrowwidth=1,
                                arrowhead=1
                              )),
                             
                              yaxis=list(
                                "showgrid"= T,
                                "title"= "",
                                "type"= "date",
                                "zeroline"= F
                              )))
    return(p)
  }

)

app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
