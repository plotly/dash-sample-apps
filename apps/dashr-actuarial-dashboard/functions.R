leeCarter = function(Years){

USALcaM<-lca(USAMortalityData, series="male", max.age=100)
USALcaF <-lca(USAMortalityData, series="female", max.age=100)
USALcaT <-lca(USAMortalityData, series="total", max.age=100)

FM<-forecast(USALcaM,h=Years)
FF<-forecast(USALcaF,h=Years)
FT<-forecast(USALcaT,h=Years)

ratesM<-cbind(USAMortalityData$rate$male[1:100,],FM$rate$male[1:100,])
ratesF<-cbind(USAMortalityData$rate$female[1:100,],FF$rate$female[1:100,])
ratesT<-cbind(USAMortalityData$rate$total[1:100,],FT$rate$total[1:100,])

plot_ly(
  x = seq(min(USAMortalityData$year),max(USAMortalityData$year)+Years),
  y = ratesM[65,],
  type = 'scatter',
  mode = 'lines',
  name = 'Male',
  fill = 'tozeroy',
  fillcolor = 'rgba(168, 216, 234, 0.5)',
  line = list(color = unlist(blues[[9]]), width = 0.5)
  
) %>%
  add_trace(  x = seq(min(USAMortalityData$year),max(USAMortalityData$year)+Years),
              y = ratesF[65,],
              type = 'scatter',
              mode = 'lines',
              name = 'Female',
              fill = 'tozeroy',
              fillcolor = 'rgba(168, 216, 234, 0.5)'
              ) %>%
  
  add_trace(  x = seq(min(USAMortalityData$year),max(USAMortalityData$year)+Years),
              y = ratesT[65,],
              type = 'scatter',
              mode = 'lines',
              name = 'Total') %>%
  layout(
    title = 'Lee-Carter Mortality Projections',
    xaxis=list(
      title= 'Time'
    ),
    yaxis=list(
      title='Delta'
    ),
    plot_bgcolor = '#F0F8FF',
    paper_bgcolor = '#F0F8FF' 
  )
}

choropleth_map <- function(){
  df = sqldf('SELECT AVG("CrudeDeathRate") AS DeathRate, State FROM mortality GROUP BY State')
  df = df[-c(9,45),]
  df$State = unlist(States)
  l <- list(color = 'white', width = 2)
  # specify some map projection/options
  g <- list(
    scope = 'usa',
    projection = list(type = 'albers usa'),
    showlakes = TRUE,
    lakecolor = 'white'
  )
  
  p <- plot_geo(df, locationmode = 'USA-states') %>%
    add_trace(
      z = ~df$DeathRate, locations = ~df$State,
      color = ~df$DeathRate, colors = 'Blues'
    ) %>%
    colorbar(title = "USD") %>%
    layout(
      title = 'Aggregate Crude Death Rate Per State',
      geo = g,
      bgcolor = '#F0F8FF'
    )
  #return(p)
}

deathrategraph = function(gender){
  data = list()
  layout = list(
    title = 'Logarithmic Death Rates by Gender',
    xaxis=list(
      title= 'Year'
    ),
    yaxis=list(
      title='Logarithmic Death Rates'
    ),
    plot_bgcolor = '#F0F8FF',
    paper_bgcolor = '#F0F8FF'
  )
  
  df = as.data.frame(USAMortalityData[["rate"]][[gender]])
  
  for (i in 1:nrow(df)) {
    new_trace = list(
      x = as.vector(colnames(df)),
      y = as.vector(unlist(unname(df[i,]))),
      name = as.character(USAMortalityData$age[i]),
      mode = 'lines+markers',
      hoverinfo = 'x+y+name',
      marker = list(symbol = "hexagram-open",
                  color = 'Blues',
                  line = list("width" = "0.5")),
      line = list(
        shape = 'spline'
      ),
      showlegend = TRUE
    )
    data = append(data, list(new_trace))
  }
  return(list('data' = data[-c(1)], 'layout' = layout))
  
}


betas$delta = c(3571.18,tail(betas$`S&P 500 High Beta Index`, -1) 
                - head(betas$`S&P 500 High Beta Index`,-1))

betagraph = function(){
  plot_ly(
    x = betas$`Effective date`[-c(1)],
    y = betas$delta[-c(1)],
    type = 'area',
    mode = 'lines',
    name = 'Male') %>%
    layout(
      margin=list(t=50),
      title = 'Changes in Market Beta Values',
      xaxis=list(
        title= 'Time'
      ),
      yaxis=list(
        title='Delta'
      ),
      plot_bgcolor = '#F0F8FF',
      paper_bgcolor = '#F0F8FF' 
    )
}

arimapred = function(SYMBOLS){

  getSymbols(SYMBOLS, from='2018-01-01', to='2019-08-30')
  # Select the relevant close price series
  stock_prices = get(SYMBOLS)[,4]
  # Compute the log returns for the stock
  stock = diff(log(stock_prices),lag=1)
  stock = stock[!is.na(stock)]
  #adf = adf.test(stock)
  #p value < 0.05. so its stationary for most. We would need to apply
  #differencing if no stationary. (stationary for this particular app assumed.)
  
  breakpoint = floor(nrow(stock)*(2.9/3))
  #acf.stock = acf(stock[c(1:breakpoint),], main='ACF Plot', lag.max=100)
  #pacf.stock = pacf(stock[c(1:breakpoint),], main='PACF Plot', lag.max=100)
  
  #Model Testing
  #model checking
  #estimate models and store for later use
  ar2 <- arima(stock,order=c(2,0,2))
  
  ar2e <-residuals(ar2)
  
  #acf(ar2e)
  #pacf(ar2e)
  ljungbox = Box.test(ar2$residuals,lag=12,type='Ljung')
  pv=1-pchisq(14.7136,11)
  
  #N-Step ahead Forecasting 
  #Forecast 10 step ahead
  return(c(predict(ar2,10)$pred))
}



