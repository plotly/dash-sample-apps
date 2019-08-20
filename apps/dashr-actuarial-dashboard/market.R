marketlayout = app$layout(
  htmlDiv(
    list(
      dccStore(id='memory-output'),
      dccDropdown(
        id = 'stockpicker',
        options=list(
          list(label = "S&P 500", value = "data/SP"),
          list(label = "Nasdaq", value = "data/^IXIC"),
          list(label = "Dow Jones Industrial Average", value = "data/^DJI")
        ),
        value="data/SP"
      ),
      dccLoading(id = 'hiii',children = htmlDiv(list(dccGraph(
        id = 'stockgraph',
        className = 'container'
      )))),
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
          columns = lapply(colnames(SP), function(x){list('name' = x, 'id' = x)}),
          style = list('float' = "left")
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
              )),
            'layout' = list(
              'title' = 'North American GDP 1947-2019',
              xaxis=list(
                title= 'Date'
              ),
              yaxis=list(
                title='GDP'
              ),
              plot_bgcolor = '#F0F8FF' ,
              paper_bgcolor = '#F0F8FF'
            )),
          className = 'container'),
        style = list('float' = "right")
      )
      )
    )
  )
)