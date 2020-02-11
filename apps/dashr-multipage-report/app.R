library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashTable)
library(stringi)
library(generator)
library(plotly)
library(random)
library(data.table)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

supplyDemand <- fread("Data/supplyDemand.csv")
actualSeasonal <- fread("Data/actualSeasonal.csv")
industrailProd <- fread("Data/industrailProd.csv")
globalMarket <- fread("Data/globalMarket.csv")
oecdCommersial <- fread("data/oecdCommersial.csv")
wtiPrices <- fread("Data/wtiPrices.csv")
epxEquity <- fread("Data/epxEquity.csv")
chinaSpr <- fread("Data/chinaSpr.csv")
oecdIndustry <- fread("Data/oecdIndustry.csv")
wtiOilprices <- fread("Data/wtiOilprices.csv")
productionCost <- fread("Data/productionCost.csv")
production2015 <- fread("Data/production2015.csv")
energyShare <- fread("Data/energyShare.csv")

set.seed(7685)
front_phone <- r_phone_numbers(5, use_hyphens = T)
front_email <- r_email_addresses(5)

# Generates and lorem ipsum script takes input number of paragraphs (n)
# and a logical input (l) TRUE of FALSE  for start point of lorem TRUE being from 
# the beginining
LipP <- function(n, l){
  stri_rand_lipsum(nparagraphs = n, start_lipsum = l)
}

color_1 <- "#003399"
color_2 <- "#00ffff"
color_3 <- "#002277"
color_b <- "#F8F8FF"

app <- Dash$new(external_stylesheets = list("https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",
                                            "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
                                            "https://codepen.io/bcd/pen/KQrXdb.css", 
                                            "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css"))
app$layout(htmlDiv(list(
  ## Page 1
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv(list(
          htmlDiv(list(
            htmlDiv(
              htmlImg(
                src  =  "assets/dash-logo-new.png",
                className = "", style = list('height' = "30px")
              )
            ),
            htmlDiv(list(
              htmlH6('Suscipit nibh'),
              htmlH5('LOREM IPSUM DOLOR'),
              htmlH6("Blandit pretium dui")
            ), style = list('color' = "white", 'margin-left' = "40px"))
          ), style = list('display' = "flex"))
        ), style = list('background-color' = color_1, 'margin-top' = "-38px", 'width' = "150%",
                        'margin-left' = "-38px", 'height' = "160px", 'padding' = "30px", 
                        'display' = "flex")),
        htmlDiv(list(
          htmlH1(list(htmlSpan('03', style = list('opacity' = '0.8')), htmlSpan('19'))),
          htmlH6('Suscipit nibh vita')
        ), style = list('color' = "white", 'background-color' = color_2,
                        'margin-top' = "-38px", 'margin-right' = "-38px", 'width' = "50%", 
                        'height' = "160px", 'padding' = "30px"))
      ), className="header", style = list('display' = "flex")),
      htmlDiv(list(
        htmlDiv(list(
          htmlH6("Felecia Conroy", style = list('color' = "#003399")),
          htmlP(front_phone[1]),
          htmlP(front_email[1])
        ), style = list('margin-top' = "30px")),
        htmlDiv(list(
          htmlH6("Olin Dach", style = list('color' = "#003399")),
          htmlP(front_phone[2]),
          htmlP(front_email[2])
        ), style = list('margin-top' = "30px")),
        htmlDiv(list(
          htmlH6("Dominique Durgan", style = list('color' = "#003399")),
          htmlP(front_phone[3]),
          htmlP(front_email[3])
        ), style = list('margin-top' = "30px")), 
        htmlDiv(list(
          htmlH6("Abraham Lemke", style = list('color' = "#003399")),
          htmlP(front_phone[4]),
          htmlP(front_email[4])
        ), style = list('margin-top' = "30px")),
        htmlDiv(list(
          htmlH6("Abraham Lemke", style = list('color' = "#003399")),
          htmlP(front_phone[5]),
          htmlP(front_email[5])
        ), style = list('margin-top' = "30px"))
      ), style = list('display' = "block", 'background-color' = color_b,
                      'width' = "250px", 'float' = "left", 'text-align' = "justify", 'padding-bottom' = "30px",
                      'horizontal-align' = "middle", 'padding-left' = "40px", 'padding-top' = "20px", 
                      'margin-left' = "15px", 'margin-top' = "30px")),
      htmlDiv(list(
        htmlDiv(list(
          htmlH6('Viverra, imperdiet, praesent pellentesque', style = list('color' = "#003399")),
          htmlP(LipP(1, F))
        ), style = list('border-left' = "5px", 'border-left-style' = "solid", 
                        "border-left-color" = color_1,
                        "border-left-width" = "7px", 'padding-left' = "20px", 'margin-top' = "30px")),
        htmlDiv(list(
          htmlH6('Facilisis mauris parturient, eget vitae', style = list('color' = "#003399")),
          htmlP(LipP(1, F))
        ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                        "border-left-color" = color_2,
                        "border-left-width" = "7px", 'padding-left' = "20px", 'margin-top' = "20px")),
        htmlDiv(list(
          htmlH6('A suspendisse mauris aliquam tincidunt hac', style = list('color' = "#003399")),
          htmlP(LipP(1, F))
        ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                        "border-left-color" = color_1,
                        "border-left-width" = "7px", 'padding-left' = "20px", 'margin-top' = "20px")),
        htmlDiv(list(
          htmlH6('A elementum lorem dolor aliquam nisi diam', style = list('color' = "#003399")),
          htmlP(LipP(1, F))
        ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                        "border-left-color" = color_2,
                        "border-left-width" = "7px", 'padding-left' = "20px", 'margin-top' = "20px"))
        
      ), style = list('display' = "block", 'float' = "right", 'width' = "400px", 'margin-left' = "10px",
                      'margin-right' = "15px"))
    ), className = "subpage")
  ), className = "page"),
  
  ## Page 2
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlH1('LOREM IPSUM')
      ), style = list('background-color' = color_1, 'color' = "white", 'width' = "200px",
                      'height' = "160px", 'padding' = "20px", 'margin-left' = "50px",
                      'margin-top' = "-38px")),  #'float' = "left"
      htmlDiv(list(
        htmlP(LipP(2, F), style = list('margin-top' = "40px")),
        htmlP(LipP(1, F), style = list('margin-top' = "20px")),
        htmlP(LipP(1, F), style = list('margin-top' = "20px"))
      ), style = list('margin-left' = "80px", 'margin-right' = "40px")), 
      htmlDiv(list(
        htmlP(LipP(1, F), style = list('margin-top' = "40px")),
        htmlP(LipP(2, F), style = list('margin-top' = "20px"))
      ), style = list('margin-left' = "80px", 'margin-right' = "40px"))
    ), className = "subpage")
  ), className = "page"),
  
  ## Page 3
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlH1('LOREM IPSUM')
      ), style = list('background-color' = color_2, 'color' = "white", 'width' = "200px",
                      'height' = "160px", 'padding' = "20px", 'margin-left' = "50px",
                      'margin-top' = "-38px")),
      htmlDiv(list(
        htmlDiv(list(
          htmlH6('Mauris feugiat quis lobortis nisl sed', style = list('color' = color_1, 'margin-top' = "40px")),
          htmlP(LipP(1, T), style = list('margin-top' = "15px"))
        )),
        htmlDiv(list(
          htmlDiv(list(
            htmlP(LipP(1, F), style = list('margin-left' = "20px"))
          ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                          "border-left-color" = color_1, 
                          "border-left-width" = "7px")),
          htmlDiv(list(
            htmlP(LipP(1, F), style = list('margin-left' = "20px"))
          ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                          "border-left-color" = color_2, 'margin-top' = "20px",
                          "border-left-width" = "7px")),
          htmlDiv(list(
            htmlP(LipP(1, F), style = list('margin-left' = "20px"))
          ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                          "border-left-color" = color_1, 'margin-top' = "20px",
                          "border-left-width" = "7px"))
        ), style = list('background-color' = "#F8F8FF", 'padding-top' = "10px", 
                        'padding-bottom' = "10px", 'padding-right' = "20px")),
        htmlDiv(list(
          htmlP(LipP(1, T), style = list('margin-top' = "20px"))
        ))
      ), className = "seven columns", style = list('margin-left' = "250px", 'margin-bottom' = "15px"))
    ), className = "subpage")
  ), className = "page"),
  
  ## Page 4
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv(list(
          htmlStrong("Ultricies fusce vel, ad ultricies enim, at, egestas", 
                     style = list('color' = color_1)),
          htmlP("Quis mauris dolor amet cubilia mattis, finibus magnis lacus",
                style = list('font-size' = "11"))
        ), className = "title six columns"),
        htmlDiv(list(
          htmlStrong("Feugiat justo, aliquam feugiat justo suspendisse leo blandit", 
                     style = list('color' = color_1)),
          htmlP("Praesent, morbi, rhoncus habitant at maximus mauris",
                style = list('font-size' = "11"))
        ), className = "title six columns")
      ), className = "thirdPage first row", style = list('position' = 'relative', 'margin-top' = '40px', 
                                                         'margin-left' = '30px')),
      htmlDiv(list(
        htmlDiv(list(
          dccGraph(
            figure = list(
              data=list(
                list(
                  x = supplyDemand$`Demand, x`,
                  y = supplyDemand$`Demand, y`,
                  hoverinfo = "y", 
                  line = list(
                    color = color_1, 
                    width = 1.5
                  ), 
                  name = "Demand", 
                  type = "scatter", 
                  id = "df2965"
                ),
                list(
                  x = supplyDemand$`Supply, x; Trace 2, x`,
                  y = supplyDemand$`Supply, y; Trace 2, y`,
                  hoverinfo = "y", 
                  line = list(
                    color = color_2, 
                    width = 1.5
                  ), 
                  name = "Supply", 
                  type = "scatter", 
                  id = "409f17"
                )
              ),
              layout = list(
                height = "250",
                'hovermode' = "closest",
                'margin' = list(
                  r = 10, 
                  t = 5, 
                  b = -1, 
                  l = 40, 
                  pad = 2
                ),
                'legend' = list(
                  x = 0.5, 
                  y = -0.4, 
                  font = list(size = 9), 
                  orientation = "h", 
                  xanchor = "center", 
                  yanchor = "bottom"
                ),
                xaxis = list(
                  range = c(1988, 2015), 
                  showgrid = FALSE, 
                  showticklabels = TRUE, 
                  tickangle = -90, 
                  tickcolor = "#b0b1b2", 
                  tickfont = list(
                    family = "Arial", 
                    size = 9
                  ), 
                  tickmode = "linear", 
                  tickprefix = "1Q", 
                  ticks = "", 
                  type = "linear", 
                  zeroline = TRUE, 
                  zerolinecolor = "#FFFFFF"
                ),
                yaxis = list(
                  autorange = FALSE, 
                  linecolor = "#b0b1b2", 
                  nticks = 9, 
                  range = c(-3000, 5000), 
                  showgrid = FALSE, 
                  showline = TRUE, 
                  tickcolor = "#b0b1b2", 
                  tickfont = list(
                    family = "Arial", 
                    size = 9
                  ), 
                  ticks = "outside", 
                  ticksuffix = " ", 
                  type = "linear", 
                  zerolinecolor = "#b0b1b2"
                )
              )
            ))
        ), className = "six columns", style = list('height' = "250px")),
        htmlDiv(list(
          dccGraph(
            figure = list(
              data = list(
                list(
                  x = actualSeasonal$`Actual, x; Crude ex. US SPR, x; Main Products, x`,
                  y = actualSeasonal$`Actual, y`,
                  hoverinfo = "y", 
                  line = list(
                    color = "#e41f23", 
                    width = 2
                  ), 
                  marker = list(
                    maxdisplayed = 0, 
                    opacity = 0
                  ), 
                  name = "Actual", 
                  type = "scatter"
                ),
                list(
                  x = actualSeasonal$`Seasonal*, x`,
                  y = actualSeasonal$`Seasonal*, y`,
                  hoverinfo = "y", 
                  line = list(
                    color = color_3, 
                    dash = "dot", 
                    width = 1.5
                  ), 
                  mode = "lines", 
                  name = "Seasonal*", 
                  type = "scatter"
                ),
                list(
                  x = actualSeasonal$`Actual, x; Crude ex. US SPR, x; Main Products, x`, 
                  y = actualSeasonal$`Crude ex. US SPR, y`, 
                  marker = list(color = color_2), 
                  name = "Crude ex. US SPR", 
                  type = "bar"
                ),
                list(
                  x = actualSeasonal$`Actual, x; Crude ex. US SPR, x; Main Products, x`, 
                  y = actualSeasonal$`Main Products, y`, 
                  marker = list(color = color_1), 
                  name = "Main Products", 
                  type = "bar"
                )
              ),
              layout = list(
                barmode =  'relative', 
                dragmode = "pan", 
                height = 250, 
                width = 310,
                hovermode = "closest", 
                legend = list(
                  x = -0.06413301662707839, 
                  y = -0.05555227415846632, 
                  bgcolor = "rgba(255, 255, 255, 0)", 
                  borderwidth = -1, 
                  font = list(size = 9), 
                  orientation = "h", 
                  traceorder = "reversed"
                ), 
                margin = list(
                  r = 10, 
                  t = 5, 
                  b = -1, 
                  l = 40, 
                  pad = 2
                ), 
                showlegend = TRUE, 
                titlefont = list(size = 16), 
                width = "550", 
                xaxis = list(
                  autorange = TRUE, 
                  range = c(0.5, 8.5), 
                  showgrid = FALSE, 
                  showline = FALSE, 
                  tickcolor = "#b0b1b2", 
                  tickfont = list(
                    family = "Arial", 
                    size = 9
                  ), 
                  tickmode = "array", 
                  ticks = "", 
                  ticktext = c("Jan-15", "Feb-15", "Mar-15", "Apr-15", "May-15", "Jun-15", "Jul-15", "Aug-15"), 
                  tickvals = c(1, 2, 3, 4, 5, 6, 7, 8), 
                  titlefont = list(size = 8), 
                  type = "linear", 
                  zeroline = TRUE, 
                  zerolinecolor = "#FFFFFF"
                ), 
                xaxis2 = list(
                  autorange = FALSE, 
                  fixedrange = TRUE, 
                  overlaying = "x", 
                  position = "0.38", 
                  range = c(0.5, 8.5), 
                  showgrid = FALSE, 
                  showticklabels = FALSE, 
                  ticks = "", 
                  ticktext = c("Jan-15", "Feb-15", "Mar-15", "Apr-15", "May-15", "Jun-15", "Jul-15", "Aug-15"), 
                  tickvals = c(1, 2, 3, 4, 5, 6, 7, 8)
                ), 
                yaxis = list(
                  autorange = FALSE, 
                  linecolor = "#b0b1b2", 
                  nticks = 8, 
                  range = c(-20, 50), 
                  showgrid = FALSE, 
                  showline = FALSE, 
                  tickcolor = "#b0b1b2", 
                  tickfont = list(
                    family = "Arial", 
                    size = 9
                  ), 
                  ticks = "outside", 
                  ticksuffix = " ", 
                  type = "linear", 
                  zerolinecolor = "#b0b1b2"
                )
              )
            )
          )
        ), className = "two columns", style = list('height' = "250px"))
      ), className="thirdPage row", style = list("padding-top" = "20px", 'margin-left' = '30px')),
      htmlDiv(list(
        htmlP("Bibendum tellus phasellus turpis sapien:"),
        htmlP(LipP(1, F), style = list('border-left' = "5px", 'border-left-style' = "solid", 'padding' = "30px",
                                       "border-left-color" = color_1, 'padding-left' = "20px", 
                                       "border-left-width" = "7px", 'background-color' = color_b))
      ), style = list('float' = "left", 'margin-top' = "20px", 'margin-left' = "30px"), className = "eleven columns"),
      htmlDiv(list(
        htmlDiv(list(
          htmlStrong("Ultricies fusce vel, ad ultricies enim, at, egestas", 
                     style = list('color' = color_1, "padding-top" = "100px")),
          htmlP("Quis mauris dolor amet cubilia mattis, finibus magnis lacus",
                style = list('font-size' = "11"))
        ), className = "title six columns"),
        htmlDiv(list(
          htmlStrong("Feugiat justo, aliquam feugiat justo suspendisse leo blandit", 
                     style = list('color' = color_1)),
          htmlP("Praesent, morbi, rhoncus habitant at maximus mauris",
                style = list('font-size' = "11"))
        ), className = "title six columns")
      ), className = "thirdPage first row", style = list('position' = 'relative', 'top' = '20px', 'margin-left' = '30px')),
      htmlDiv(list(
        htmlDiv(list(
          dccGraph(
            figure = list(
              data = list(
                trace1 <- list(
                  x = industrailProd$`Industrial Production, x`,
                  y = industrailProd$`Industrial Production, y`,
                  line = list(color = color_2), 
                  mode = "lines", 
                  name = "Industrial Production", 
                  type = "scatter", 
                  visible = TRUE
                ),
                trace2 <- list(
                  x = industrailProd$`Price (rhs), x`,
                  y = industrailProd$`Price (rhs), y`,
                  line = list(color = color_1), 
                  mode = "lines", 
                  name = "Price (rhs)", 
                  type = "scatter",  
                  visible = TRUE, 
                  yaxis = "y2"
                )
              ),
              layout = list(
                annotations = list(
                  list(
                    x = 0.95, 
                    y = -0.15, 
                    arrowhead = 7, 
                    ax = 0, 
                    ay = -40, 
                    font = list(size = 8), 
                    showarrow = FALSE, 
                    text = "months after shock", 
                    xref = "paper", 
                    yref = "paper"
                  )
                ), 
                autosize = TRUE, 
                dragmode = "pan", 
                height = "300",
                width = "300",
                hovermode = "closest", 
                legend = list(
                  x = 0.00805970149254, 
                  y = 1.03496503497, 
                  bgcolor = "rgba(255, 255, 255, 0)", 
                  font = list(size = 9)
                ), 
                margin = list(
                  r = 40, 
                  t = 5, 
                  b = -1, 
                  l = 20, 
                  pad = 0
                ), 
                paper_bgcolor = "rgba(0, 0, 0, 0)", 
                plot_bgcolor = "rgba(0, 0, 0, 0)", 
                showlegend = TRUE, 
                xaxis = list(
                  autorange = FALSE, 
                  nticks = 19, 
                  range = c(0.5, 18), 
                  showgrid = FALSE, 
                  tickfont = list(
                    color = "rgb(68, 68, 68)", 
                    size = "9"
                  ), 
                  ticks = "", 
                  type = "linear", 
                  zeroline = FALSE
                ), 
                yaxis = list(
                  autorange = FALSE, 
                  linecolor = "rgb(190, 191, 192)", 
                  mirror = TRUE, 
                  nticks = "9", 
                  range = c(-0.4, 1.2), 
                  showgrid = FALSE, 
                  showline = TRUE, 
                  side = "left", 
                  tickfont = list(
                    color = "rgb(68, 68, 68)", 
                    size = 9
                  ), 
                  ticks = "outside", 
                  ticksuffix = " ", 
                  type = "linear", 
                  zeroline = FALSE
                ), 
                yaxis2 = list(
                  anchor = "x", 
                  autorange = FALSE, 
                  exponentformat = "e", 
                  linecolor = "rgb(190, 191, 192)", 
                  nticks = 9, 
                  overlaying = "y", 
                  range = c(-0.1, 0.3), 
                  showgrid = FALSE, 
                  side = "right", 
                  tickfont = list(size = 9), 
                  tickprefix = " ", 
                  ticks = "outside", 
                  type = "linear", 
                  zerolinecolor = "rgb(190, 191, 192)"
                )
              )
            )
          )
        ), className = "six columns", style = list('height' = "250px")),
        htmlImg(src = "assets/IJsVT9P.png", 
                className = "exhibit six columns", 
                style = list('margin-top' = "20px"))
      ), className="", style = list("padding-top" = "20px", 'margin-left' = '30px'))
    ), className = "subpage")
  ), className = "page"),
  
  # Page 5
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv(list(
          htmlP(LipP(1, F))
        ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                        "border-left-color" = color_1, 'padding' = "20px",
                        "border-left-width" = "7px")),
        htmlDiv(list(
          htmlP(LipP(1, F))
        ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                        "border-left-color" = color_2, 'padding' = "20px",
                        "border-left-width" = "7px", 'margin-top' = "20px")),
        htmlDiv(list(
          htmlP(LipP(1, F))
        ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                        "border-left-color" = color_1, 'padding' = "20px",
                        "border-left-width" = "7px", 'margin-top' = "20px"))
      ), style = list('background-color' = color_b, 'margin-left' = "25px", 'margin-top' = "40px"
      ), className = "eleven columns row"),
      htmlP(LipP(1, F), style = list('margin-left' = "5px", 'margin-top' = "20px", 'margin-bottom' = "20px"), 
            className = "twelve columns row"),
      htmlDiv(list(
        htmlDiv(list(
          htmlStrong("Ultricies fusce vel, ad ultricies enim, at, egestas", 
                     style = list('color' = color_1)),
          htmlP("Quis mauris dolor amet cubilia mattis, finibus magnis lacus",
                style = list('font-size' = "11"))
        ), className = "title six columns"),
        htmlDiv(list(
          htmlStrong("Feugiat justo, aliquam feugiat justo suspendisse leo blandit", 
                     style = list('color' = color_1)),
          htmlP("Praesent, morbi, rhoncus habitant at maximus mauris",
                style = list('font-size' = "11"))
        ), className = "title six columns")
      ), className = "thirdPage first row", style = list('position' = 'relative', 'top' = '0px',
                                                         'margin-bottom' = "40px")),
      
      htmlDiv(list(
        htmlDiv(list(
          dccGraph(
            figure = list(
              data = list(
                list(
                  x = globalMarket$x,
                  y = globalMarket$y,
                  marker = list(color = color_1), 
                  name = "Global market imbalance", 
                  type = "bar"
                )
              ),
              layout = list( 
                dragmode = "pan", 
                bargap = 0.63,
                height = "250",
                width = "320",
                hovermode = "closest", 
                legend = list(
                  x = 0.0006061953460797935, 
                  y = -0.31665440684852813, 
                  bgcolor = "rgba(255, 255, 255, 0)", 
                  borderwidth = -1, 
                  font = list(size = 9),
                  orientation = "h"
                ), 
                margin = list(
                  r = 40, 
                  t = 5, 
                  b = -1, 
                  l = 20, 
                  pad = 0
                ),
                showlegend = TRUE, 
                title = "Click to enter Plot title", 
                xaxis = list(
                  autorange = TRUE, 
                  exponentformat = "none", 
                  nticks = 18, 
                  range = c(-0.5, 15.5), 
                  showline = FALSE, 
                  tickangle = "auto", 
                  tickfont = list(size = 9), 
                  tickmode = "linear", 
                  ticks = "", 
                  title = "Click to enter X axis title", 
                  type = "category"
                ), 
                yaxis = list(
                  autorange = TRUE, 
                  linecolor = "rgb(176, 177, 178)", 
                  nticks = 10, 
                  range = c(-1283.8982436029166, 3012.5614936594166), 
                  showgrid = FALSE, 
                  showline = TRUE, 
                  tickangle = "auto", 
                  tickfont = list(size = 9), 
                  ticks = "outside", 
                  title = "", 
                  type = "linear", 
                  zeroline = TRUE, 
                  zerolinecolor = "rgb(176, 177, 178)"
                )
              )
            )
          )
        ), className = "six columns", style = list('height' = "250px")),
        htmlDiv(list(
          dccGraph(
            figure = list(
              data = list(
                list(
                  x = oecdCommersial$`OECD commercial ex. US NGL & other, x`,
                  y = oecdCommersial$`OECD commercial ex. US NGL & other, y`,
                  line = list(color = color_1), 
                  mode = "lines", 
                  name = "OECD commercial ex. US NGL & other", 
                  type = "scatter"
                ),
                list(
                  x = oecdCommersial$`Seasonal (2000-2014), x`,
                  y = oecdCommersial$`Seasonal (2000-2014), y`,
                  line = list(color = color_2), 
                  mode = "lines", 
                  name = "Seasonal (2000-2014)", 
                  type = "scatter"
                )
              ),
              layout = list(
                height = "250", 
                width = "320",
                annotations = list(
                  list(
                    x = 0.9122213607288229, 
                    y = 1.0056975042758463, 
                    showarrow = FALSE, 
                    text = "GS forecast", 
                    xref = "paper", 
                    yref = "paper"
                  )
                ),
                hovermode = "closest", 
                legend = list(
                  x = 0.22370614460166696, 
                  y = -0.21665440684852813, 
                  font = list(size = 9), 
                  orientation = "h"
                ), 
                margin = list(
                  r = 40, 
                  t = 5, 
                  b = -1, 
                  l = 20, 
                  pad = 0
                ),
                paper_bgcolor = "rgba(255, 255, 255, 0)", 
                plot_bgcolor = "rgba(255, 255, 255, 0)", 
                shapes = list(
                  list(
                    line = list(
                      color = "rgb(0, 0, 0)", 
                      dash = "dot", 
                      width = 1
                    ), 
                    type = "line", 
                    x0 = 0.6208955223880597, 
                    x1 = 0.6208955223880597, 
                    xref = "paper", 
                    y0 = 0, 
                    y1 = 1, 
                    yref = "paper"
                  )
                ), 
                showlegend = TRUE, 
                xaxis = list(
                  autorange = FALSE, 
                  linecolor = "rgb(190, 191, 192)", 
                  nticks = 17, 
                  range = c(0.5, 16), 
                  showgrid = FALSE, 
                  showline = FALSE, 
                  tickfont = list(size = 9), 
                  ticks = "", 
                  ticksuffix = " ", 
                  title = "", 
                  type = "category", 
                  zeroline = FALSE, 
                  zerolinecolor = "rgb(190, 191, 192)"
                ), 
                yaxis = list(
                  autorange = FALSE, 
                  linecolor = "rgb(190, 191, 192)", 
                  nticks = 10, 
                  range = c(-800, 1000), 
                  showgrid = FALSE,
                  showline = TRUE,
                  tickfont = list(size = 10), 
                  ticks = "outside", 
                  ticksuffix = " ", 
                  title = "", 
                  type = "linear",
                  zeroline = TRUE,
                  zerolinecolor = "rgb(190, 191, 192)"
                )
              )
            )
          )
        ), className = "six columns", style = list('height' = "250px"))
      ), className="thirdPage row", style = list('margin-top' = "30px"))
    ), className = "subpage")
  ), className = "page"),
  
  # Page 6
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlP(LipP(1, F), style = list('float' = "left"), 
              className = "six columns"),
        htmlP(LipP(1, F), style = list('float' = "left", 'margin-left' = "5px"), 
              className = "six columns")
      ), className = "thirdPage row", style = list('position' = 'relative', 'margin-top' = '40px')),
      
      htmlDiv(list(
        htmlStrong("Vehicula elementum congue penatibus massa, eu sed", className = "eleven columns",
                   style = list('color' = color_1)),
        htmlDiv(
          htmlImg(src = "assets/DBkxRT2.png", 
                  style = list('margin-top' = "20px", 'height' = "250px"))
        )
      ), className = "eleven columns"),
      htmlDiv(list(
        htmlStrong("At velit pharetra ac fusce sit dictum pellentesque", className = "eleven columns", 
                   style = list('color' = color_1)),
        htmlDiv(list(
          dccGraph(
            figure = list(
              data = list(
                list(
                  x = wtiPrices$`WTI Prices, x`,
                  y = wtiPrices$`WTI Prices, y`,
                  line = list(
                    color = color_1, 
                    dash = "solid"
                  ), 
                  mode = "lines", 
                  name = "WTI Prices", 
                  type = "scatter", 
                  uid = "1380d4"
                ),
                list(
                  x = wtiPrices$`Sep-15 forecast, x`,
                  y = wtiPrices$`Sep-15 forecast, y`,
                  line = list(color = "rgb(228, 31, 35)"), 
                  mode = "lines", 
                  name = "Sep-15 forecast", 
                  type = "scatter", 
                  uid = "766493"
                ),
                list(
                  x = wtiPrices$`Forward, x`,
                  y = wtiPrices$`Forward, y`,
                  line = list(
                    color = color_2, 
                    dash = "solid"
                  ), 
                  mode = "lines", 
                  name = "Forward", 
                  type = "scatter", 
                  uid = "4e8248"
                ),
                list(
                  x = wtiPrices$`May-15 forecast, x`,
                  y = wtiPrices$`May-15 forecast, y`,
                  line = list(
                    color = color_3,
                    dash = "solid"
                  ), 
                  mode = "lines", 
                  name = "Forward", 
                  type = "scatter"
                )
              ),
              layout = list(
                height = 250, 
                hovermode = "closest", 
                legend = list(
                  x = 0.16039179104479998, 
                  y = 1, 
                  bgcolor = "rgba(255, 255, 255, 0)", 
                  bordercolor = "rgba(68, 68, 68, 0)", 
                  font = list(
                    color = "rgb(68, 68, 68)", 
                    size = 10
                  ), 
                  orientation = "h", 
                  traceorder = "normal"
                ), 
                margin = list(
                  r = 40, 
                  t = 5, 
                  b = -1, 
                  l = 20, 
                  pad = 0
                ), 
                showlegend = TRUE, 
                xaxis = list(
                  autorange = FALSE, 
                  gridwidth = 0, 
                  linecolor = "rgb(130, 132, 134)", 
                  mirror = FALSE, 
                  nticks = 14, 
                  range = c(0, 14), 
                  showgrid = FALSE, 
                  showline = TRUE, 
                  side = "bottom", 
                  tickfont = list(
                    color = "rgb(68, 68, 68)", 
                    size = 9
                  ), 
                  ticks = "outside", 
                  ticktext = c("Sep-14", "Nov-14", "Jan-15", "Mar-15", "May-15", "Jul-15", "Sep-15", 
                               "Nov-15", "Jan-16", "Mar-16", "May-16", "Jul-16", "Sept-16", "Nov-16"), 
                  tickvals = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), 
                  title = "", 
                  type = "linear", 
                  zeroline = FALSE, 
                  zerolinecolor = "rgb(130, 132, 134)"
                ), 
                yaxis = list(
                  autorange = FALSE, 
                  linecolor = "rgb(130, 132, 134)", 
                  nticks = 8, 
                  range = c(30, 100), 
                  showline = TRUE, 
                  tickfont = list(
                    color = "rgb(68, 68, 68)", 
                    size = 9
                  ), 
                  ticks = "outside", 
                  ticksuffix = " ", 
                  title = "", 
                  type = "linear", 
                  zerolinecolor = "rgb(130, 132, 134)"
                )
              )
              
            )
          )
        ))
      ), className = "eleven columns", style = list('margin-top' = "20px"))
    ), className = "subpage")
  ), className = "page"),
  
  # Page 7
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv(list(
          htmlP(LipP(1, F), style = list('float' = "left")),
          htmlP(LipP(1, F), style = list('float' = "left")),
          htmlP(LipP(1, F), style = list('float' = "left")),
          htmlP(LipP(1, F), style = list('float' = "left"))
        ), className = "seven columns", style = list('position' = 'relative', 'margin-top' = '40px',
                                                     'margin-left' = "5px")),
        htmlDiv(list(
          htmlDiv(list(
            htmlStrong("Vehicula elementum congue penatibus massa, eu sed sed dolor", style = list('color' = color_1)),
            htmlDiv(list(
              dccGraph(figure=list(
                data=list(
                  list(
                    x = c("AAA", "AA", "A", "BBB", "BB", "B", "CCC"), 
                    y = c("1497", "976", "1016", "1739", "993", "545", "31"),
                    type='bar',
                    name="y",
                    marker = list(color = color_1)
                  )
                ),
                layout = list(
                  height = "300",
                  autosize = TRUE, 
                  bargap = 0.75, 
                  legend = list(
                    x = 0.16039179104479998, 
                    y = -0.2720578174979476, 
                    bgcolor = "rgba(255, 255, 255, 0)", 
                    bordercolor = "rgba(68, 68, 68, 0)", 
                    font = list(
                      color = "rgb(68, 68, 68)", 
                      size = 10
                    ), 
                    orientation = "h", 
                    traceorder = "normal"
                  ),
                  hovermode = "closest", 
                  margin = list(
                    r = 0, 
                    t = 10, 
                    b = 30, 
                    l = 60
                  ), 
                  xaxis = list(
                    autorange = FALSE, 
                    nticks = 10, 
                    range = c(-0.5, 6.5), 
                    tickfont = list(size = 9), 
                    ticks = "", 
                    title = "", 
                    type = "category"
                  ), 
                  yaxis = list(
                    autorange = FALSE, 
                    dtick = "250", 
                    nticks = 9, 
                    range = c(0, 2250), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickfont = list(size = 9), 
                    ticks = "outside", 
                    title = "2015E production by rating (mboe)<br><br>", 
                    titlefont = list(size = 9), 
                    type = "linear", 
                    zeroline = TRUE
                  )
                )
              ))
            ))
          ), style = list('margin-top' = "5px", height="250")),
          htmlDiv(list(
            htmlStrong("At velit pharetra ac fusce sit dictum pellentesque, dictumst", style = list('color' = color_1)),
            htmlDiv(
              dccGraph(
                figure = list(
                  data = list(
                    list(
                      x = epxEquity$`EPX equity sector, x`,
                      y = epxEquity$`EPX equity sector, y`,
                      line = list(
                        color = color_1, 
                        width = 2
                      ), 
                      mode = "lines", 
                      name = "EPX equity sector", 
                      type = "scatter", 
                      uid = "1305fa", 
                      visible = TRUE
                    ),
                    list(
                      x = epxEquity$`WTI 2-yr swap, x`,
                      y = epxEquity$`WTI 2-yr swap, y`,
                      line = list(
                        color = color_2, 
                        width = 2
                      ), 
                      mode = "lines", 
                      name = "WTI 2-yr swap", 
                      type = "scatter", 
                      uid = "2ce7c1", 
                      visible = TRUE
                    ),
                    list(
                      x = epxEquity$`HY energy spread ratio (rhs, inverted), x`,
                      y = epxEquity$`HY energy spread ratio (rhs, inverted), y`,
                      line = list(
                        color = "red", 
                        width = 2
                      ), 
                      mode = "lines", 
                      name = "HY energy spread ratio (rhs, inverted)", 
                      type = "scatter", 
                      uid = "400a20", 
                      visible = TRUE, 
                      yaxis = "y2"
                    )
                  ),
                  layout = list(
                    height = "300",
                    autosize = TRUE, 
                    hovermode = "closest", 
                    legend = list(
                      x = 0.008033242860512229, 
                      y = -0.3007047167087806, 
                      bgcolor = "rgba(255, 255, 255, 0)", 
                      font = list(
                        color = "rgb(68, 68, 68)", 
                        size = 9
                      ), 
                      orientation = "h"
                    ), 
                    margin = list(
                      r = 30, 
                      t = 10, 
                      b = 20, 
                      l = 30
                    ), 
                    showlegend = TRUE, 
                    xaxis = list(
                      autorange = FALSE, 
                      linecolor = "rgb(130, 132, 134)", 
                      linewidth = 1, 
                      range = c(0, 12), 
                      showgrid = FALSE, 
                      showline = TRUE, 
                      tickfont = list(size = 9), 
                      ticks = "outside", 
                      ticktext = c("Sep-14", "Oct-14", "Nov-14", "Dec-14", "Jan-15", "Feb-15", "Mar-15", "Apr-15", "May-15", "Jun-15", "July-15", "Aug-15"), 
                      tickvals = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), 
                      title = "", 
                      type = "linear", 
                      zeroline = FALSE
                    ), 
                    yaxis = list(
                      autorange = FALSE, 
                      linecolor = "rgb(130, 132, 134)", 
                      nticks = 8, 
                      range = c(30, 100), 
                      showgrid = FALSE, 
                      showline = TRUE, 
                      tickfont = list(size = 9), 
                      ticks = "outside", 
                      title = "", 
                      type = "linear", 
                      zeroline = FALSE
                    ), 
                    yaxis2 = list(
                      anchor = "x", 
                      linecolor = "rgb(130, 132, 134)", 
                      nticks = 10, 
                      overlaying = "y", 
                      range = c(1.8, 0.9), 
                      showgrid = FALSE, 
                      showline = TRUE, 
                      side = "right", 
                      tickfont = list(size = 9), 
                      ticks = "outside", 
                      title = "Click to enter Y axis title", 
                      type = "linear", 
                      zeroline = FALSE
                    )
                  )
                ))
            )
          ), style = list('margin-top' = "30px"))
        ), className="eight columns", style = list('padding' = "30px", 'margin-bottom' = "30px", 
                                                   'position' = 'relative', 'margin-top' = '10px', 
                                                   'margin-left' = "30px", 'border-left' = "10px",
                                                   'border-left-style' = "solid", 'border-left-color' = color_b,
                                                   'border-left-width' = "2px"))
      ), style = list('display' = "flex"))
    ), className = "subpage")
  ), className = "page"),
  
  # Page 8
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv(list(
          htmlH6('Aliquet ut mauris nostra habitant egestas, massa vulputate. Magnis nullam leo eget ullamcorper lacus congue laoreet ex sed', 
                 style = list('color' = color_1, 'margin-top' = "40px")),
          htmlP(LipP(1, F), style = list('margin-top' = "15px")),
          htmlP(LipP(1, F), style = list('margin-top' = "15px", 'color' = color_1))
        )),
        htmlDiv(list(
          htmlDiv(list(
            htmlP(LipP(1, F), style = list('margin-left' = "20px"))
          ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                          "border-left-color" = color_1, 
                          "border-left-width" = "7px")),
          htmlDiv(list(
            htmlP(LipP(1, F), style = list('margin-left' = "20px"))
          ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                          "border-left-color" = color_2, 'margin-top' = "20px",
                          "border-left-width" = "7px")),
          htmlDiv(list(
            htmlP(LipP(1, F), style = list('margin-left' = "20px"))
          ), style = list('border-left' = "5px", 'border-left-style' = "solid",
                          "border-left-color" = color_1, 'margin-top' = "20px",
                          "border-left-width" = "7px"))
        ), style = list('background-color' = "#F8F8FF", 'padding-top' = "10px", 
                        'padding-bottom' = "10px", 'padding-right' = "20px")),
        htmlDiv(list(
          htmlP(LipP(1, F), style = list('margin-top' = "20px"))
        ))
      ), className = "nine columns", style = list('margin-left' = "100px"))
    ), className = "subpage")
  ), style = list('margin-top' = "50px"), className = "page"),
  
  # Page 9
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv(list(
          htmlP("Aenean felis et libero nullam pretium quis est in sit. Commodo nec ante aenean a. Commodo at facilisis vestibulum cursus elementum nascetur et, placerat class aliquam convallis porttitor accumsan. Ultricies sed laoreet eleifend maximus venenatis", 
                style = list('color' = color_1)),
          htmlStrong("Congue nisl iaculis interdum cubilia maximus"),
          htmlImg(src = "assets/wX5mQYn.png", 
                  className = "exhibit eleven columns", style = list('height' = "250", 'margin-top' = "20px"))
        ), style = list('float' = "left")),
        htmlDiv(list(
          htmlP("Id nulla sollicitudin taciti ac tempus amet ligula accumsan. Elementum, nullam dui ligula ut. Adipiscing sed ultricies ut vitae augue etiam nostra nibh.", 
                style = list(color = color_1)),
          htmlStrong("Convallis et eu habitant leo leo luctus venenatis"),
          htmlDiv(list(
            dccGraph(
              figure = list(
                data = list(
                  list(
                    x = chinaSpr$`OECD commercial ex. US NGL & other, x`,
                    y = chinaSpr$`OECD commercial ex. US NGL & other, y`,
                    line = list(
                      color = color_1, 
                      width = 2
                    ), 
                    mode = "lines", 
                    name = "OECD commercial ex. US NGL & other", 
                    type = "scatter", 
                    uid = "57b5cc", 
                    visible = TRUE
                  ),
                  list(
                    x = chinaSpr$`Non-OECD stocks ex. China SPR, x`,
                    y = chinaSpr$`Non-OECD stocks ex. China SPR, y`,
                    line = list(
                      color = color_2, 
                      width = 2
                    ), 
                    mode = "lines", 
                    name = "Non-OECD stocks ex. China SPR", 
                    type = "scatter", 
                    uid = "bcd3b2"
                  )
                ),
                layout = list(
                  annotations = list(
                    list(
                      x = 12.0815219907062, 
                      y = 948.201438849, 
                      font = list(size = 9), 
                      showarrow = FALSE, 
                      text = "GS forecast", 
                      xref = "x", 
                      yref = "y"
                    )
                  ), 
                  autosize = TRUE,
                  height = "300",
                  dragmode = "zoom", 
                  hovermode = "closest", 
                  legend = list(
                    x = 0.0913178294574, 
                    y = -0.167832167832, 
                    bgcolor = "rgba(255, 255, 255, 0)", 
                    font = list(size = 9), 
                    orientation = "h"
                  ), 
                  margin = list(
                    r = 10, 
                    t = 10, 
                    b = 0, 
                    l = 40, 
                    pad = 0
                  ), 
                  shapes = list(
                    list(
                      line = list(
                        color = "rgb(68, 68, 68)", 
                        dash = "dot", 
                        width = 1
                      ), 
                      type = "line", 
                      x0 = 0.6541331802525385, 
                      x1 = 0.6541331802525385, 
                      xref = "paper", 
                      y0 = 0, 
                      y1 = 1, 
                      yref = "paper"
                    )
                  ), 
                  showlegend = TRUE, 
                  xaxis = list(
                    autorange = FALSE, 
                    nticks = 10, 
                    range = c(-0.25, 15.5), 
                    showgrid = FALSE, 
                    showline = FALSE, 
                    tickfont = list(size = 9), 
                    ticktext = c("1Q13", "2Q13", "3Q13", "4Q13", "1Q14", "2Q14", "3Q14", "4Q14", "1Q15", "2Q15", "3Q15E", "4Q15E", "1Q16E", "2Q16E", "3Q16E", "4Q16E"), 
                    tickvals = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), 
                    title = "", 
                    type = "linear", 
                    zeroline = FALSE, 
                    zerolinecolor = "rgb(130, 132, 134)", 
                    zerolinewidth = 1
                  ), 
                  yaxis = list(
                    autorange = FALSE, 
                    nticks = 10, 
                    range = c(-800, 1000), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickfont = list(
                      color = "rgb(68, 68, 68)", 
                      size = 9
                    ), 
                    ticks = "outside", 
                    title = "", 
                    type = "linear", 
                    zeroline = TRUE
                  )
                )
              ))
          ), style = list('margin-top' = "30px"))
        ), style = list('margin-top' = "50px", 'float' = "left"))
      ), style = list('margin-left' = "30px", 'margin-top' = "40px", 'margin-left' = "30px"), className = "exibit six columns"),
      htmlDiv(list(
        htmlP(LipP(1, F)),
        htmlP(LipP(1, F)),
        htmlP(LipP(1, F)),
        htmlP(LipP(1, F))
      ), className = "five columns", style = list('margin-top' = "40px"))
    ), style =  list('display' = "flex", 'flex-wrop' = "wrap"),  className = "subpage")
  ), className = "page"),
  
  # Page 10
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv(list(
          htmlStrong("Nulla diam conubia nec lacus urna in ligula nec ut egestas sed. Diam inceptos nec venenatis",
                     style = list('color' = color_1)),
          htmlP("Nulla diam conubia nec lacus urna in ligula nec ut egestas sed", style = list("font-size" = "11")),
          htmlDiv(list(
            dccGraph(
              figure = list(
                data = list(
                  list(
                    x = oecdIndustry$`OECD industry stock changes , x`,
                    y = oecdIndustry$`OECD industry stock changes , y`,
                    line = list(color = color_1), 
                    mode = "lines", 
                    name = "OECD industry stock changes ", 
                    type = "scatter", 
                    uid = "fd77c7"
                  ),
                  list(
                    x = oecdIndustry$`IEA miscellaneous to balance (rhs), x`,
                    y = oecdIndustry$`IEA miscellaneous to balance (rhs), y`,
                    line = list(color = color_2), 
                    mode = "lines", 
                    name = "IEA miscellaneous to balance (rhs)", 
                    type = "scatter", 
                    uid = "562a09", 
                    yaxis = "y2"
                  )
                ),
                layout = list(
                  height = "250",
                  autosize = TRUE, 
                  hovermode = "closest", 
                  legend = list(
                    x = 0.13880597014925372, 
                    y = -0.25129840919314606, 
                    bgcolor = "rgba(255, 255, 255, 0)", 
                    orientation = "h"
                  ), 
                  margin = list(
                    r = 30, 
                    t = 10, 
                    b = 0, 
                    l = 30
                  ), 
                  shapes = list(
                    list(
                      fillcolor = "rgba(31, 119, 180, 0)", 
                      line = list(
                        color = "rgb(255, 0, 0)", 
                        dash = "dash", 
                        width = 1
                      ), 
                      opacity = 1.1, 
                      type = "square", 
                      x0 = 1997.25, 
                      x1 = 1998.75, 
                      xref = "x", 
                      y0 = -1713.7349397590363, 
                      y1 = 2391.5662650602408, 
                      yref = "y"
                    ), 
                    list(
                      fillcolor = "rgba(31, 119, 180, 0)", 
                      layer = "above", 
                      line = list(
                        color = "rgb(255, 0, 0)", 
                        dash = "dash", 
                        width = 1
                      ), 
                      opacity = 1.1, 
                      type = "square", 
                      x0 = 2013.25, 
                      x1 = 2014.75, 
                      xref = "x", 
                      y0 = -1674.2105263157894, 
                      y1 = 2286.315789473684, 
                      yref = "y"
                    )
                  ), 
                  showlegend = TRUE, 
                  xaxis = list(
                    autorange = FALSE, 
                    nticks = 30, 
                    range = c(1986, 2015), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickangle = "auto", 
                    tickfont = list(size = 8), 
                    tickprefix = "1Q", 
                    ticks = "outside", 
                    type = "linear", 
                    zeroline = TRUE
                  ), 
                  yaxis = list(
                    autorange = FALSE, 
                    nticks = 10, 
                    range = c(-2000, 2500), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickfont = list(size = 8), 
                    ticks = "outside", 
                    type = "linear"
                  ), 
                  yaxis2 = list(
                    anchor = "x", 
                    autorange = FALSE, 
                    nticks = 12, 
                    overlaying = "y", 
                    range = c(-2500, 2500), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    side = "right", 
                    tickfont = list(size = 8), 
                    ticks = "outside", 
                    type = "linear", 
                    zeroline = FALSE
                  )
                )
              )
            )
          ))
        ), className = "thirdPage first row", style = list('margin-top' = '0px')),
        htmlDiv(list(
          htmlStrong("Risus amet quam, eget, lacus, orci, dui facilisis dolor sodales arcu facilisi consectetur",
                     style = list('color' = color_1)),
          htmlP("Diam, maximus ultricies neque adipiscing tellus eros proin", style = list("font-size" = "11")),
          htmlDiv(list(
            dccGraph(
              figure = list(
                data = list(
                  list(
                    x = wtiOilprices$x,
                    y = wtiOilprices$y,
                    line = list(color = color_1), 
                    name = "WTI oil prices (S/bbl, 2015 $)    ", 
                    type = "scatter"
                  )
                ),
                layout = list(
                  height = "250",
                  autosize = TRUE, 
                  legend = list(
                    x = 0.16818221960553428, 
                    y = -0.30969810073003856, 
                    bgcolor = "rgba(255, 255, 255, 0)", 
                    font = list(size = 9)
                  ), 
                  margin = list(
                    r = 10, 
                    t = 10, 
                    b = 40, 
                    l = 30
                  ), 
                  shapes = list(
                    list(
                      fillcolor = "rgba(31, 119, 180, 0)", 
                      line = list(
                        color = "rgb(255, 0, 0)",
                        dash = "dash",
                        width = 1
                      ), 
                      opacity = 1.1, 
                      type = "square", 
                      x0 = 1985.6994029850746, 
                      x1 = 1987.4587313432835, 
                      xref = "x", 
                      y0 = 10, 
                      y1 = 85, 
                      yref = "y"
                    ), 
                    list(
                      fillcolor = "rgba(31, 119, 180, 0)", 
                      layer = "above", 
                      line = list(
                        color = "rgb(253, 0, 0)", 
                        dash = "dash",
                        width = 1
                      ), 
                      opacity = 1.1, 
                      type = "square", 
                      x0 = 1998.1650746268656, 
                      x1 = 1999.989328358209, 
                      xref = "x", 
                      y0 = 5, 
                      y1 = 70, 
                      yref = "y"
                    )
                  ), 
                  showlegend = TRUE, 
                  xaxis = list(
                    autorange = FALSE, 
                    nticks = 24, 
                    range = c(1972, 2015.5), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickfont = list(size = 9), 
                    ticks = "outside", 
                    titlefont = list(color = "rgb(92, 53, 143)"), 
                    type = "linear", 
                    zeroline = FALSE
                  ), 
                  yaxis = list(
                    autorange = FALSE, 
                    nticks = 11, 
                    range = c(0, 180), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickfont = list(size = 9), 
                    ticks = "outside", 
                    type = "linear"
                  )
                )
              )
            )
          ))
        ), className = "thirdPage first row", style = list('margin-top' = '20px'))
      ), style = list('margin-top' = '40px')),
      htmlDiv(list(
        htmlDiv(list(
          htmlStrong("Porttitor felis eget nibh quam duis et at a massa varius.",
                     style = list('color' = color_1)),
          htmlP("Risus amet quam, eget, lacus, orci, dui facilisis ", style = list("font-size" = "11")),
          htmlDiv(list(
            dccGraph(
              figure = list(
                data = list(
                  list(
                    x = productionCost$x,
                    y = productionCost$y,
                    line = list(color = color_1), 
                    type = "scatter"
                  )
                ),
                layout = list(
                  height = "200",
                  margin = list(
                    r = 20, 
                    t = 10, 
                    b = 50, 
                    l = 40
                  ), 
                  xaxis = list(
                    autorange = FALSE, 
                    exponentformat = "none", 
                    linecolor = "rgb(171, 172, 173)", 
                    nticks = 5, 
                    range = c(0, 40000), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickfont = list(size = 9), 
                    ticks = "outside", 
                    title = "Cumulative peak oil production (kb/d)", 
                    titlefont = list(size = 9), 
                    type = "linear", 
                    zeroline = FALSE
                  ), 
                  yaxis = list(
                    autorange = FALSE, 
                    linecolor = "rgb(171, 172, 173)", 
                    nticks = 10, 
                    range = c(0, 45), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickfont = list(size = 9), 
                    ticks = "outside", 
                    title = "Production cost (US$/bbl)", 
                    titlefont = list(size = 9), 
                    type = "linear", 
                    zeroline = FALSE
                  )
                )
              )
            )
          ))
        ), className = "six columns"),
        htmlDiv(list(
          htmlStrong("Arcu aenean litora quam dignissim penatibus sem ultrices",
                     style = list('color' = color_1)),
          htmlP("Aenean ipsum nostra magna ut sagittis venenatis", style = list("font-size" = "11")),
          htmlDiv(list(
            dccGraph(
              figure = list(
                data = list(
                  list(
                    x = production2015$`Canadian Producers, x`,
                    y = production2015$`Canadian Producers, y`,
                    marker = list(
                      color = "rgb(255, 0, 0)", 
                      symbol = "diamond"
                    ), 
                    mode = "markers", 
                    name = "Canadian Producers", 
                    type = "scatter", 
                    visible = TRUE
                  ),
                  list(
                    x = production2015$`US E&Ps and Integrated, x`,
                    y = production2015$`US E&Ps and Integrated, y`,
                    marker = list(
                      color = color_2, 
                      symbol = "diamond"
                    ), 
                    mode = "markers", 
                    name = "US E&Ps and Integrated", 
                    type = "scatter"
                  ),
                  list(
                    x = production2015$`Others, x`,
                    y = production2015$`Others, y`,
                    marker = list(
                      color = color_1, 
                      symbol = "diamond"
                    ), 
                    mode = "markers", 
                    name = "Others", 
                    type = "scatter"
                  )
                ),
                layout = list(
                  height = "200",
                  autosize = TRUE, 
                  hovermode = "closest", 
                  legend = list(
                    x = -0.06, 
                    y = -0.36, 
                    font = list(size = 9), 
                    orientation = "h"
                  ), 
                  margin = list(
                    r = 10, 
                    t = 10, 
                    b = 0, 
                    l = 40
                  ), 
                  showlegend = TRUE, 
                  xaxis = list(
                    autorange = FALSE, 
                    nticks = 8, 
                    range = c(0, 100), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickfont = list(size = 9), 
                    ticks = "outside", 
                    ticksuffix = "%", 
                    title = "2015 Net Debt / Capital Employed", 
                    titlefont = list(size = 9), 
                    type = "linear", 
                    zeroline = FALSE
                  ), 
                  yaxis = list(
                    autorange = FALSE, 
                    nticks = 12, 
                    range = c(0, 45), 
                    showgrid = FALSE, 
                    showline = TRUE, 
                    tickfont = list(size = 9), 
                    ticks = "outside", 
                    title = "2015 Production Cost $/bbl", 
                    titlefont = list(size = 9), 
                    type = "linear", 
                    zeroline = FALSE
                  )
                )
              )
            )
          ))
        ), className = "six columns")
      ), className = "thirdPage first row", style = list('margin-top' = '20px'))
      
    ), className = "subpage")
  ), className = "page"),
  
  # Page 11
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlH6('In tempor mauris non, maximus non odio. Lacus mi arcu, ut parturient ac sed curae \
               sed litora amet quam, massa purus condimentum', 
               style = list('color' = color_1, 'margin-top' = "40px")),
        htmlP(LipP(1, F), style = list('margin-top' = "10px")),
        htmlP(LipP(1, F), style = list('margin-top' = "10px", 'color' = color_1)),
        htmlP(LipP(1, F), style = list('margin-top' = "10px"))
      ), className = "eleven columns"),
      htmlDiv(list(
        htmlP('Non amet tempor pellentesque facilisis velit, dui nulla hendrerit sociosqu fusce', 
              style = list('color' = color_1, 'margin-top' = "10px")),
        dccGraph(
          figure = list(
            data = list(
              list(
                x = energyShare$`Energy share of HY Issuance, x`,
                y = energyShare$`Energy share of HY Issuance, y`,
                line = list(color = color_1), 
                mode = "lines", 
                name = "Energy share of HY Issuance", 
                type = "scatter",
                uid = "203858", 
                visible = TRUE
              ),
              list(
                x = energyShare$`US oil rig count  (monthly change), x`,
                y = energyShare$`US oil rig count  (monthly change), y`,
                line = list(color = color_2), 
                mode = "lines", 
                name = "US oil rig count  (monthly change)", 
                type = "scatter", 
                uid = "e17c2f",
                yaxis = "y2"
                
              )
            ),
            layout = list(
              height = "300px",
              autosize = TRUE, 
              hovermode = "closest", 
              legend = list(
                x = 0.39727646537238737, 
                y = -0.12197967025477964, 
                bgcolor = "rgba(255, 255, 255, 0)", 
                font = list(
                  color = "rgb(68, 68, 68)", 
                  size = 9
                ), 
                orientation = "h", 
                traceorder = "reversed"
              ), 
              margin = list(
                r = 30, 
                t = 10, 
                b = 0, 
                l = 30
              ), 
              showlegend = TRUE, 
              xaxis = list(
                autorange = TRUE, 
                nticks = 10, 
                range = c(-0.007132542, 8.1854778101), 
                showgrid = FALSE, 
                tickfont = list(size = 9), 
                ticks = "", 
                ticktext = c(" Jan-13", "  May-13", "  Sep-13", "  Jan-14", "  May-14", "  Sep-14", "  Jan-15", "  May-15"), 
                tickvals = c(0, 1, 2, 3, 4, 5, 6, 7), 
                title = "", 
                type = "linear", 
                zeroline = TRUE, 
                zerolinecolor = "rgb(171, 172, 173)", 
                zerolinewidth = 1
              ), 
              yaxis = list(
                autorange = FALSE, 
                linecolor = "rgb(136, 137, 140)", 
                nticks = 10, 
                range = c(-300, 150), 
                showgrid = FALSE, 
                showline = FALSE, 
                tickfont = list(size = 9), 
                ticks = "outside", 
                title = "", 
                type = "linear", 
                zeroline = TRUE, 
                zerolinecolor = "rgb(171, 172, 173)", 
                zerolinewidth = 1
              ), 
              yaxis2 = list(
                anchor = "x", 
                autorange = FALSE, 
                linecolor = "rgb(136, 137, 140)", 
                nticks = 8, 
                overlaying = "y", 
                range = c(0, 35), 
                showgrid = FALSE, 
                showline = TRUE, 
                side = "right", 
                tickfont = list(size = 9), 
                ticks = "outside", 
                ticksuffix = " %", 
                type = "linear", 
                zeroline = FALSE, 
                zerolinecolor = "rgb(171, 172, 173)", 
                zerolinewidth = 1
              )
            )
          )
        )
      ), className = "eleven columns")
    ), className = "subpage")
  ), style = list('margin-top' = "50px"), className = "page"),
  
  # Page 12
  htmlDiv(list(
    htmlDiv(list(
      htmlDiv(list(
        htmlH6("Erat cras porta inceptos nibh sociis justo. Natoque mauris nunc etiam, dis quam, tempor consectetur ac \
               Pulvinar nunc vitae dui elit hac ante, facilisi, primis nascetur. Non nostra torquent ipsum ac amet", 
               style = list('color' = color_1, 'margin-top' = "40px")),
        htmlP(LipP(1, F), style = list('margin-top' = "30px")),
        htmlP(LipP(2, F), style = list('margin-top' = "30px")),
        htmlH6("Ultrices phasellus dignissim, accumsan platea volutpat, sapien mi enim. Pharetra ipsum netus in turpis, \
               lorem tempus et. Eget sed. Eu porta cum tempor convallis sed nostra, pellentesque eros.",
               style = list('color' = color_1, 'margin-top' = "30px")),
        htmlImg(src = "assets/c7PM25P.png", style = list('margin-top' = "30px"), className = "twelve columns")
      ), className = "twelve columns")
    ), className = "subpage")
  ), style = list('margin-top' = "50px"), className = "page")
)))

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}



