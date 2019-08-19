layout = app$layout(
  htmlDiv(
    list(
      htmlH1('Market Beta Visualization'),
      dccGraph(
        figure = betagraph()
      )
    )),
  htmlH3('Generalized Logistic Regression on SP500'),
  htmlDiv(list(
    htmlLabel('Volume '),
    dccInput(id='input-1'),
    
    htmlLabel('Year '),
    dccInput(id='input-2'),
    
    htmlButton("Predict", id="button-2"),
    
    dccLoading(id = 'hi2', children = htmlDiv(id='output'))
  )),
  htmlBr(),
  htmlDiv(list(
    dccDropdown(
      id = 'stockpicker',
      options=lapply(tickersymbols, function(symbol){list(label = symbol, value = symbol)}),
      value = list("AMZN", 'AAPL', 'GOOG', 'FB'),
      multi = TRUE,
      style = list('width' = 'auto')
    ),
    
    htmlDiv(list(
      dccLoading(id = 'hi',
                 children = htmlDiv(list(dccGraph(id = 'Arima',
                 style = list('float' = 'left',
                            'width' = '1000px',
                            'height' = '600px'))
               ))),
      htmlDiv(list(
        htmlH3('Black-Scholes Price Model',
               style = list('color' = '#003366')),
        htmlLabel('Stock Price (S0)'),
        htmlBr(),
        
        dccInput(id='S0'),
        htmlBr(),
        
        htmlLabel('Strike'),
        htmlBr(),
        
        dccInput(id='K'),
        htmlBr(),
        
        htmlLabel('Risk Free Interest Rate(r)'),
        htmlBr(),
        
        dccInput(id='r'),
        htmlBr(),
        
        htmlLabel('Time'),
        htmlBr(),
        
        dccInput(id='time'),
        htmlBr(),
        
        htmlLabel('Volatility'),
        htmlBr(),
        
        dccInput(id='vola'),
        htmlBr(),
        
        dccRadioItems(
          id = 'type',
          options=list(
            list("label" = "Call", "value" = "Call"),
            list("label" = "Put", "value" = "Put")
          ),
          value = "Call",
          labelStyle = list("display" = "inline-block")
        ),
        htmlBr(),
        
        
        htmlButton("Predict", id="bsbutton"),
        htmlBr(),
        dccLoading(id= 'hi3', children = htmlDiv(id='BSOutput'))
      ), style = list("border" = "2px black solid",
                      'background-color' = '#E5E5E5',
                      'width' = '30vw',
                      'float' = 'right',
                      'padding' = '30px'))
      
      ), style = list('display' = 'flex'))
  )
  ))

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



Smarket$Direction = as.factor(Smarket$Direction)
names(Smarket)
glm.fit <- glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + 
                 Lag5 + Volume, data = Smarket, 
               family = binomial)


train = Smarket[Smarket$Year < 2005,]
test = Smarket[!(Smarket$Year < 2005),]

glm.fit <- glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume, 
               data = Smarket, 
               family = binomial, subset = Smarket$Year<2005)

summary(glm.fit)

glm.probs <- predict(glm.fit, 
                     newdata = test, 
                     type = "response")

glm.pred <- ifelse(glm.probs > 0.5, "Up", "Down")

Direction.2005 = test$Direction

table(glm.pred, Direction.2005)

mean(glm.pred == Direction.2005)


#we dont have a good accuracy. lets improve by using a smaller
#model since lag 4 and 5 dont have good p values.


glm.fit = glm(formula = Direction ~ Lag1 + Lag2 + Lag3 + Volume + Year, family = binomial, 
              data = Smarket)

summary(glm.fit)

# 
# glm.probs <- predict(glm.fit, 
#                      newdata = test, 
#                      type = "response")
# 
# glm.pred <- ifelse(glm.probs > 0.5, "Up", "Down")
# 
# Direction.2005 = test$Direction
# 
# table(glm.pred, Direction.2005)
# 
# mean(glm.pred == Direction.2005)
# 

newdata = data.frame(Lag1 = 0.381, Lag2 = -0.192, Lag3 = -2.624, 
                     Volume = 2, Year = 2020)

predict(glm.fit, newdata, type = 'response')




