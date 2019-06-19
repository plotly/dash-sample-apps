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
          dccMarkdown(id='text', className='six columns')
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

###################################################

###################################################
#CALLBACK BEGINS
###################################################

###################################################
#CALLBACK ENDS
###################################################
  app$run_server()
  

