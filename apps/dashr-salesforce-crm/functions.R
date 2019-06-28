library(dash)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dplyr)
library(sqldf)
library(plotly)
library(salesforcer)
library(rjson)

millnames = list("", " K", " M", " B", " T")

app = Dash$new()

millify <- function(n){
  n = n
  millidx = max(0, min(length(millnames)-1, integer(floor(
    if (n==0){
      0
    } else {
      log10(abs(n))/3
    }
  ))))
  
  sprintf('%f', (n / 10 ** (3 * millidx)))
}


indicator <- function(color, text, id_value){
  return(htmlDiv(
    list(htmlP(
        text,
        className="twelve columns indicator_text"
      ),
      htmlP(
        id = id_value,
        className="indicator_value"
      )),
    className="four columns indicator"  
  ))}