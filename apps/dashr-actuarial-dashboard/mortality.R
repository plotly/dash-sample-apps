mortalitylayout = app$layout(
  htmlDiv(
    list(
      htmlHr(),
      dccGraph(
        id="map",
        figure = choropleth_map(),
        style=list("height"= "90%", "width"= "98%"),
        config=list(displayModeBar=FALSE),
        className = 'container'
      ),
      htmlH3('lee-Carter Projection Time'),
      htmlDiv(dccInput(
        id = 'projections',
        placeholder = "Enter Projection Period...",
        type = "number",
        value = 5
      )
      ))
    ),
  htmlDiv(list(
    dccLoading(id = 'testt', children = htmlDiv(list(dccGraph(
      id = 'LeeCarter',
      figure = leeCarter(20),
      #style = list('width' = '800px', 'height' = '600px'),
      className = 'container'
    )))),
    dccRadioItems(
      id = 'genderpicker',
      options=list(
        list("label" = "Male", "value" = "male"),
        list("label" = "Female", "value" = "female")
      ),
      value = "female"
    ),
    dccLoading(id = 'LCC', 
               children = htmlDiv(list(dccGraph(
                 id = 'LC',
                 #style = list('width' = '800px', 'height' = '600px'),
                 className = 'container'))), type="default"
    )
  ),
  style = list('display' = 'flex'))
  )
