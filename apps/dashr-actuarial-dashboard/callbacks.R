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
        name = row
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
                 return(dccLoading(children = blackSc$Vega))}
             })
