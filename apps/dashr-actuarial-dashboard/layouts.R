marketlayout = app$layout(
  htmlDiv(
    list(
      dccStore(id='memory-output'),
      dccDropdown(
        id = 'stockpicker',
        options=list(
          list(label = "S&P 500", value = "sp"),
          list(label = "Nasdaq", value = "^IXIC"),
          list(label = "Dow Jones Industrial Average", value = "^DJI")
        ),
        value="sp"
      ),
      dccGraph(
        id = 'stockgraph'
      ),
      htmlBr(),
      dccDatePickerRange(
        id = 'datepicker',
        display_format = "MMMM Y, DD",
        start_date_placeholder_text = "MMMM Y, DD",
        start_date = format(as.Date("2019-01-01"), format  =  "%Y,%m,%d"),
        end_date = format(as.Date("2019-01-31"), format  =  "%Y,%m,%d")
      ),
      htmlHr(),
      htmlDiv(list(
        dashDataTable(
          id='memory-table',
          columns= lapply(colnames(SP), function(x){list('name' = x, 'id' = x)})
          ),
        dccGraph(
          figure=list(
            data=list(
              list(
                x=gdp$Index,
                y=gdp$GDPC1,
                type='bar',
                name='SF',
                marker=list(
                  color= '#119DFF'
                )
              ))))
        ),
        style = list(display = "flex")
      )
    )
  )
)

mortalitylayout = app$layout(
  htmlDiv(
    list(
      htmlHr(),
dccGraph(
  id="map",
  figure = choropleth_map(),
  style=list("height"= "90%", "width"= "98%"),
  config=list(displayModeBar=FALSE)
),
dccRadioItems(
  id = 'genderpicker',
  options=list(
    list("label" = "Male", "value" = "male"),
    list("label" = "Female", "value" = "female")
  ),
  value = "female",
  labelStyle = list("display" = "inline-block")
),
dccGraph(
  id = 'LC'
),
htmlH4('Projection Horizon'),
dccInput(
  id = 'projections',
  placeholder = "Enter Projection Period...",
  type = "number",
  value = 5
),
dccGraph(
  id = 'LeeCarter',
  figure = leeCarter(20)
)
)))

risklayout = app$layout(
  htmlDiv(
    list(
      htmlH1('Market Beta Visualization'),
      dccGraph(
        figure = betagraph()
      )
    )),
  htmlDiv(list(
    htmlLabel('Volume '),
    dccInput(id='input-1'),
    
    htmlLabel('Year '),
    dccInput(id='input-2'),
    
    htmlButton("Predict", id="button-2"),
    
    htmlDiv(id='output')
  )),
  htmlBr(),
  htmlDiv(list(
    
    htmlDiv(list(
                 dccDropdown(
                   id = 'stockpicker',
                   options=lapply(tickersymbols, function(symbol){list(label = symbol, value = symbol)}),
                   value = list("AMZN", 'AAPL', 'GOOG', 'FB'),
                   multi = TRUE,
                   style = list('width' = 'auto')
                 ),
      dccGraph(id = 'Arima',
               style = list('float' = 'left',
                            'width' = '1100px',
                            'height' = '600px')))),
    
  
    
    
    htmlDiv(list(
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
      htmlDiv(id='BSOutput')
    ), style = list("border" = "2px black solid",
                    'background-color' = '#E5E5E5',
                    'width' = '30vw',
                    'float' = 'right',
                    'padding' = '50px'))
  )
))




