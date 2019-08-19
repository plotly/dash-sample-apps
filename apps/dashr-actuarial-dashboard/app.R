library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashTable)
library(data.table)
library(lifecontingencies)
library(demography)
library(forecast)
library(sqldf)
library(plotly)
library(caret)
library(ISLR)
library(readxl)
library(StMoMo)
library(quantmod)
library(ragtop)

data("Smarket")

source('risk.R')
source('mortality.R')
source('market.R')
source('functions.R')

setwd('/Users/kevinphan/desktop/dash-sample-apps/apps/dashr-actuarial-dashboard')

app = Dash$new()

States = list('AL', 'AK', 'AZ', 'AR', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'HI', 'ID', 'IL', 'IN', 'IA', 'KS', 'KY', 'LA',
              'ME', 'MD', 'MA', 'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 'NV', 'NH', 'NJ', 'NM', 'NY', 'NC', 'ND', 'OH', 'OK',
              'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 'WA', 'WV', 'WI', 'WY')

tickersymbols = list('MSFT', 'AAPL', 'AMZN', 'FB', 'GOOG',
                     'JPM', 'XOM', 'V', 'BAC', 'INTC', '
                     CSCO', 'DIS', 'HD', 'VZ', 'CVX', 'MA')

blues = list('#00688B', '#009ACD', '#0099CC', '#00B2EE', 
             '#00BFFF', '#BFEFFF', '#33A1C9', '#507786', 
             '#87CEEB','#38B0DE', '#0BB5FF', '#6996AD', 
             '#236B8E', '#325C74', '#517693', '#4682B4', 
             '#525C65', '#62B1F6', '#0D4F8B') 

# SP = read.csv('data/sp.csv')
# dowjones = read.csv('data/^DJI.csv')
# nasdaq = read.csv('data/^IXIC.csv')
# gdp = read.csv('data/GDP.csv')
# betas = read_xls('data/PerformanceGraphExport.xls')
# 
# mortality = read.csv('NCHS_-_Drug_Poisoning_Mortality_by_State__United_States.csv')
# USAMortalityData<-hmd.mx(country="USA", username="YOUR USERNAME",
#                          password="YOUR PASSWORD", label="U.S.A.")
# 
# premiums = read_excel('data/premiums.xlsx',sheet = 1)

indexcolor = list(
  'SP'= list('#119DFF'),
  "^IXIC" = list('black'),
  "^DJI" = list('red'))

tab_style = list(
  'borderBottom'= '1px solid #d6d6d6',
  'padding'= '6px',
  'fontWeight'= 'bold'
)

tab_selected_style = list(
  'backgroundColor'= '#119DFF',
  'color'= 'white',
  'padding'= '6px'
)


app$layout(htmlDiv(list(
  htmlImg(
    src="/assets/logo.png",
    height = '70px'
  ),
  htmlH1('ACTUARIAL RISK DASHBOARD'),
  dccTabs(id='Market/Economic Data', value='Market/Economic Data', children=list(
    dccTab(label='Market/Economic Data', value='Market/Economic Data', style=tab_style, selected_style=tab_selected_style),
    dccTab(label='Mortality', value='Mortality', style=tab_style, selected_style=tab_selected_style),
    dccTab(label='Risk', value='Risk', style=tab_style, selected_style=tab_selected_style)
  )),
  htmlDiv(id='tabs-content-example')
)))

app$callback(
  output = list(id='tabs-content-example', property = 'children'),
  params = list(input(id = 'Market/Economic Data', property = 'value')),
  function(tab){
    if(tab == 'Market/Economic Data'){
      return(marketlayout)
    }
    
    else if(tab == 'Mortality'){
      return(dccLoading(mortalitylayout))
    }
    else if(tab == 'Risk'){
      return(layout)
    }
  }
)

app$callback(output = list(id='output', property = 'children'),
             params = list(input('button-2', 'n_clicks'),
                           state('input-1', 'value'),
                           state('input-2', 'value')),
             function(n_clicks, input1, input2){
               #browser()
               newdata = data.frame(Lag1 = 0.0038344, Lag2 = 0.0039192, Lag3 = 0.001716, 
                                    Volume = as.numeric(input1), Year = as.numeric(input2))
               return(sprintf('The probability of the S&P 500 going up is %s percent',
                              100 * round(predict(glm.fit, newdata, 
                                                  type = 'response'), digits =4)))
             })
app$callback(
  output = list(id='Arima', property = 'figure'),
  params = list(input(id =  'stockpicker', property = 'value')),
  function(stock){
    layout = list(
      title = 'Log Returns Prediction of Stocks Using Arima',
      xaxis=list(
        title= 'Time Lag'
      ),
      yaxis=list(
        title='Logarithmic Returns'
      )
    )
    #browser()
    aggregation = list()
    stock = as.vector(unlist(stock))
    j = 1
    for (row in stock) {
      aggregation[[j]] <- list(
        x = c(1:10),
        y = c(unlist(arimapred(row))),
        text = row,
        mode = 'lines+markers',
        name = row,
        marker = list(color = blues[[j]])
      )
      j = j + 1
    }
    return(list(
      'data' = aggregation,
      'layout' = layout))
  }
)


app$callback(output = list(id='BSOutput', property = 'children'),
             params = list(input("bsbutton", 'n_clicks'),
                           state('S0', 'value'),
                           state('K', 'value'),
                           state('r', 'value'),
                           state('time', 'value'),
                           state('vola', 'value'),
                           state('type', 'value')
             ),
             function(button,S0,K,r,time,vola, type){
               variables = list(S0[[1]],K[[1]],r[[1]],time[[1]],vola[[1]], type[[1]])
               
               if(TRUE %in% c(unlist(lapply(variables, is.null)))){
                 return(sprintf('Error! Missing Parameter!'))
               } else if("" %in% variables){
                 return(sprintf('Error! Missing Parameter!'))
               } else{
                 if(type == 'Call'){
                   type = 1
                 } else if(type == 'Put'){
                   type= -1
                 }
                 blackSc = blackscholes(callput=as.numeric(type), S0=as.numeric(S0), 
                                        K=as.numeric(K), r=as.numeric(r), 
                                        time=as.numeric(time),
                                        vola=as.numeric(vola), default_intensity=0.07)
                 return(sprintf('Delta: %g | \n 
                                Vega: %s | \n 
                                Price: %s',
                                round(bc$Delta, 2), round(bc$Vega, 2), round(bc$Price,2)))}
             })

app$callback(
  output = list(id='LC', property = 'figure'),
  params = list(input(id = 'genderpicker', property = 'value')),
  function(gender){
    if(gender == 'male'){
      return(deathrategraph('male'))
    }
    
    else if(gender == 'female'){
      return(deathrategraph('female'))
    }
  }
)   



app$callback(
  output = list(id='LeeCarter', property = 'figure'),
  params = list(input(id = 'projections', property = 'value')),
  function(projection){
    return(leeCarter(projection))
  }
)   

app$callback(
  output = list(id='stockgraph', property='figure'),
  params = list(input(id='stockpicker', property='value')),
  function(selected_value){
    dat = read.csv(paste0(selected_value,'.csv'))
    data = list(
      list(
        type='Scatter GL',
        x=dat$Date,
        y=dat$Open,
        marker=list(
          color=indexcolor[[selected_value]][[1]]
        )
      )
    )
    return(list(
      'data' = data,
      'layout' = list(
        title = 'Stock Trend',
        xaxis = list(title = "Date",
                     showgrid = FALSE,
                     showline = FALSE,
                     zeroline = FALSE),
        yaxis = list(title = "Index",
                     showgrid = FALSE,
                     showline = FALSE,
                     ticks = 'outside',
                     zeroline = FALSE)),
      
      hovermode = 'closest',
      gridcolor = 'rgb(255,255,255)'
    ))
  })

app$callback(
  output = list(id="memory-output", property = 'data'),
  params = list(input(id = "stockpicker", property = 'value'),
                input(id = "datepicker", property = 'start_date'),
                input(id = "datepicker", property = 'end_date')
  ),
  
  function(selected_value, start_date, end_date){
    df = read.csv(paste0(selected_value,'.csv'))
    df$Date = as.Date(df$Date)
    
    if(grepl(',', start_date) == TRUE){
      start = as.Date(start_date, format = "%Y,%m,%d")}
    else{
      start = start_date
    }
    if(grepl(',', end_date) == TRUE){
      end = as.Date(end_date, format = "%Y,%m,%d")}
    else{
      end = end_date
    }
    
    df = df[df$Date >= start & df$Date <= end,]
    
    if(length(selected_value) < 1){
      return(df_to_list(df))
    }
    return(df_to_list(df))
  })

app$callback(
  output = list(id="memory-table", property = 'data'),
  params = list(input(id = "memory-output", property = 'data')),
  function(data){
    if(is.null(data) == TRUE){
      return()
    }
    return(data)
  })

app$run_server()   
