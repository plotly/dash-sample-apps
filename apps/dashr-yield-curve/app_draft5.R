library(plotly)
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(quantmod)
library(magrittr)
library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)
library(glue)
library("aws.s3")
library(jsonlite)
library(rlang)

l. <- list

#Helper functions

genMarksEmpty <- function(n){
  l <- list(glue(''))
  names(l) <- 'label'
  m <- list(l)
  names(m) <- glue('{n}')
  return(m)
}



###############INTERNAL LOGIC######################

last_back <- 0
last_next <- 0

DT_ <- fread('data/yield_curve.csv') %>% drop_na(., x,y)
dates <- DT_$y
labels <-(DT_$x)[1:11]
DT <- DT_[, 3:13]
colnames(DT) <- labels
rownames(DT) <- dates
DM <- DT %>% data.matrix()

UPS <- l.('0'= l.(x=0, y=0, z=1), '1'=l.(x=0, y=0, z=1), '2'=l.(x=0, y=0, z=1), '3'=l.(x=0, y=0, z=1), '4'=l.(x=0, y=0, z=1), '5'=l.(x=0, y=0, z=1))

CENTERS <- l.('0'= l.(x=0.3, y=0.8, z=-0.5), 
              '1'=l.(x=0, y=0, z=-0.37), 
              '2'=l.(x=0, y=1.1, z=-1.3), 
              '3'=l.(x=0, y=-0.7, z=0), 
              '4'=l.(x=0, y=-0.2, z=0), 
              '5'=l.(x=-0.11, y=-0.5, z=0))

EYES <-  l.('0'= l.(x=2.7, y=2.7, z=0.3), 
            '1'=l.(x=0.01, y=3.8, z=-0.37), 
            '2'=l.(x=1.3, y=3, z=0), 
            '3'=l.(x=2.6, y=-1.6, z=0), 
            '4'=l.(x=3, y=-0.2, z=0), 
            '5'=l.(x=-0.1, y=-0.5, z=2.66))

TEXTS <- l.(
  l.(
    htmlDiv(dccMarkdown(className='six columns', children=" #### Yield curve 101 "), style=l.(width='120%')), #htmlBr(),
    htmlDiv(dccMarkdown(className='six columns', children=" The yield curve shows how much it costs the federal government to borrow
                 money for a given amount of time, revealing the relationship between long-
                   and short-term interest rates.  It is, inherently, a forecast for what the economy holds in the future -
                 how much inflation there will be, for example, and how healthy growth will
                 be over the years ahead - all embodied in the price of money today,
                 tomorrow and many years from now. "), style=l.(width='100%'))), 
  
  l.(htmlDiv(dccMarkdown(className='six columns', children= "#### Where we stand"), style=l.(width='120%')), #htmlBr(),
     htmlDiv( dccMarkdown(className='six columns', children="On Wednesday, both short-term and long-term rates were lower than they have
    been for most of history - a reflection of the continuing hangover
    from the financial crisis. The yield curve is fairly flat, which is a sign that investors expect
    mediocre growth in the years ahead."), style=l.(width='100%'))), 
  l.(htmlDiv(dccMarkdown(className='six columns', children=' #### Deep in the valley'),style=l.(width='120%')), #htmlBr(),
     htmlDiv(dccMarkdown(className='six columns', children='In response to the last recession, the Federal Reserve has kept short-term
    rates very low - near zero - since 2008. (Lower interest rates stimulate
    the economy, by making it cheaper for people to borrow money, but also
    spark inflation.) Now, the Fed is getting ready to raise rates again, possibly as early as
    June.'), style=l.('width'='100%'))), 
  l.(htmlDiv(dccMarkdown(className='six columns', children=' #### Last time, a puzzle'), style=l.(width='120%')), #htmlBr(),
     htmlDiv(dccMarkdown(className='six columns', children='The last time the Fed started raising rates was in 2004. From 2004 to 2006,
    short-term rates rose steadily. But long-term rates did not rise very much. The Federal Reserve chairman called this phenomenon a "conundrum," and it
    raised questions about the ability of the Fed to guide the economy.
    Part of the reason long-term rates failed to rise was because of strong
    foreign demand.') , style=l.('width'='100%'))), 
  l.(htmlDiv(dccMarkdown(className='six columns', children='#### Long-term rates are low now, too'), style=l.(width='120%')), #htmlBr(),
     htmlDiv(dccMarkdown(className='six columns', children= 'Foreign buyers have helped keep long-term rates low recently, too - as have
    new rules encouraging banks to hold government debt and expectations that
    economic growth could be weak for a long time. The 10-year Treasury yield was as low as it has ever been in July 2012 and
    has risen only modestly since.Some economists refer to the economic pessimism as "the new normal."'), style=l.('width'='110%'))), 
  l.(htmlDiv(dccMarkdown(className='six columns', children=' #### Long-term rates are low now, too'),style=l.(width='120%') ), #htmlBr(),
     htmlDiv(dccMarkdown(className='six columns', children= 'Here is the same chart viewed from above.'), style=l.(width='70vh')))
)

##INITALIZE THE GRAPH
graph <- DM %>% plot_ly(x=labels, y=base::as.Date(dates, format='%m/%d/%Y'),z=., type='surface', 
                        hoverinfo='x+y+z', lighting=l.('ambient'=0.95, 'diffuse'=0.99, 'fresnel'=0.01, 'roughness'=0.01, 'specular'=0.01 ),
                        colorscale=l.(l.(0, 'rgb(230,245,254)'), l.(0.4,'rgb(123,171,203)'), l.(0.8, 'rgb(40,119,174)', l.(1, 'rgb(37,61,81)')), zmax=9.18, zmin=0, scene='scene')
                        )
#graph


#Setup the app
app <- Dash$new()

###################################################
#LAYOUT BEGINS
###################################################

app$layout(
  
  htmlDiv(l.(
    #
    htmlDiv(l.(
      ##
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
              
      ##
            ,
            className='row'
            ,
            style=l.('text-align'='center', 'margin-bottom'= '15px')

            )
    #
    ,

    htmlDiv(l.(
      ##
      htmlDiv(l.(
        ###
          dccSlider(min=0, max=5, value=0, marks=lapply(0:5, genMarksEmpty), id='slider')
        ###
          ), className='row',style=l.('margin-bottom'= '10px'))
      ##
      ,
      ##
      htmlDiv(l.(
        ###
          htmlDiv(l.(
            ####
              htmlButton('Back', id='back', style=l.('display'='inline-block'),n_clicks=0)#, className='learn-more-button')
            ####
              ,
            ####
              htmlButton('Next', id='next', style=l.('display'='inline-block'),n_clicks=0)#,  className='learn-more-button')
            ####
              ),className='two columns offset-by-two')
        ###
          ,
        ###
          htmlDiv(id='text', className='six columns', children='')
        ###

          ), className='row', style=l.('margin-bottom'= '10px'))
      ##
      ,
      ##
      dccGraph(
        id='graph',
        style=l.('height'= '60vh'),
        figure=graph
      )
      ##
      ,
      ##
      htmlDiv(id='bucket', style=l.(display='none'))
      ,
      htmlDiv(id='dummy')
     
      
      
      ##
    ), id='page')
    #
     
    # #


  ))
  
)
###################################################
#LAYOUT ENDS
###################################################



###################################################

###################################################
#CALLBACK BEGINS
###################################################

## BUTTON CONTROLS

####HELPER function

soft_bd  <- function(n){
  if(n<0){return(0)} else if(n >5){return(5)}else{return(n)}
}

app$callback(
  output=l.(id='bucket', property='children')
  ,
  params = l.(input(id='back', property='n_clicks'), input(id='next', property='n_clicks'))
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
    
    return(V %>% toJSON(., null='null'))
  }
)


# 
app$callback(

  output = l.(id='slider', property='value')
  ,
  params = l.(input(id='bucket', property='children'), state(id='slider', property='value'))
  ,
  function(V, S){
    value <- (V %>%fromJSON)
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
  
  output=l.(id='dummy', property='children')
  ,
  params=l.(input(id='bucket', property='children'))
  ,
  function(n){
    return(n)
  }
)

## MAKE ANNOTATIONS

app$callback(
  
  output=l.(id='text', property='children')
  ,
  params=l.(input(id='slider', property='value'))
  ,
  function(value){
    
    if(is.null(value)||is.na(value)){value=0} else{value}
    
    return(TEXTS[[(value+1)]])
    
  }
) 

# ##LOAD GRAPH
# app$callback(
#   
#   output=l.(id='graph', property='figure')
#   ,
#   params=l.(input(id='slider', property='value'))
#   ,
#   function()
#   {
#     
#   }
#   
# )




#

###################################################
#CALLBACK ENDS
###################################################
 
app$run_server()

