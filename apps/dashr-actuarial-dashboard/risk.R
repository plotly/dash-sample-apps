layout = app$layout(
  htmlDiv(
    list(
      htmlH1('Market Beta Visualization'),
      dccGraph(
        figure = betagraph(),
        className = 'container'
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
    htmlHr(),
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
                            'height' = '600px'),
                 className = 'container')
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
      ), style = list("border" = "2px #F0F8FF solid",
                      'background-color' = '#F0F8FF',
                      'width' = '30vw',
                      'float' = 'right',
                      'padding' = '30px'),
      className = 'container')
      
      ), style = list('display' = 'flex'))
  )
  ))

