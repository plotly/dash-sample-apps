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
      )
      )
    )
  ,
  htmlH4('Projection Horizon'),
  dccInput(
    id = 'projections',
    placeholder = "Enter Projection Period...",
    type = "number",
    value = 5
  ),
  htmlDiv(list(
    dccLoading(id = 'LCC', 
               children = htmlDiv(list(dccGraph(
                 id = 'LC'))), type="default"
    ),
    dccLoading(id = 'testt', children = htmlDiv(list(dccGraph(
      id = 'LeeCarter',
      figure = leeCarter(20)
    ))))
  ),
  style = list('display' = 'flex'))
  )
