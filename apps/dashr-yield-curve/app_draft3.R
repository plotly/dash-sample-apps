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

l. <- list

#Helper functions

genMarksEmpty <- function(n){
  l <- list(glue(''))
  names(l) <- 'label'
  m <- list(l)
  names(m) <- glue('{n}')
  return(m)
}


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
              htmlButton('Back', id='back', style=l.('display'='inline-block'))#, className='learn-more-button')
            ####
              ,
            ####
              htmlButton('Next', id='next', style=l.('display'='inline-block'))#,  className='learn-more-button')
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
        style=l.('height'= '60vh')
      )
      ##
    ), id='page')
    #
  ))
  
)
###################################################
#LAYOUT ENDS
###################################################

###############INTERNAL LOGIC######################

last_back <- 0
last_next <- 0

DT_ <- fread('data/yield_curve.csv')

DT <- drop_na(DT_, x,y)


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

###################################################

###################################################
#CALLBACK BEGINS
###################################################

## MAKE ANNOTATIONS

app$callback(
  
  output=l.(id='text', property='children')
  ,
  params=l.(input(id='slider', property='value'))
  ,
  function(value){
    
    if(is.null(value)||is.na(value)){value=0}
    
    return(TEXTS[[(value+1)]])
    
  }
  
)

###################################################
#CALLBACK ENDS
###################################################
  app$run_server()
  

