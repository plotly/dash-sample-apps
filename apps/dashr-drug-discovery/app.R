library(dash)
library(plotly)
library(dashCoreComponents)
library(dashHtmlComponents)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("./app/apps/%s", appName))
}

setwd(sprintf("./app/apps/%s", appName))

df <- read.csv("data/small_molecule_drugbank.csv")


### Graph & Layout Objects ###
colorscale <-   list(
    list(0, "rgb(54,50,153)"),
    list(0.3, "rgb(17,123,215)"),
    list(0.4, "rgb(37,180,167)"),
    list(0.5, "rgb(134,191,118)"),
    list(0.65, "rgb(249,210,41)"),
    list(1, "rgb(244,236,21")
)

appTitle <- htmlH3("dash for drug discovery", className = "uppercase title")
subtitle1Bold <- htmlSpan("Hover ", className = "uppercase bold")
subtitle1 <-  htmlSpan("over a drug in the graph to see its structure.")
subtitle2Bold <- htmlSpan("Select ", className = "uppercase bold")
subtitle2 <- htmlSpan(
  "a drug in the dropdown to add it to the drug candidates at the bottom."
)

threeD <- plot_ly(
  df,
  x = ~ PKA,
  y = ~ SOL,
  z = ~ LOGP,
  text = df$NAME,
  mode = "markers",
  marker = list(
    color = ~ MW,
    sizeref = 55,
    size = ~ MW,
    colorscale = colorscale,
    colorbar = list(title = "Molecular<br>Weight"),
    showscale = TRUE,
    symbol = "circle",
    sizemode = "diameter",
    line = list(color = "#444")
  )
) %>% layout(
  scene = list(
    xaxis = list(title = "pkA"),
    yaxis = list(title = "LogP", type = "log"),
    zaxis = list(title = "Solubility (mg/ml)")),
  plot_bgcolor = "#d8d8d8"
  )

histo <- plot_ly(df, x =  ~ PKA, y = ~ LOGP) %>% add_histogram2d(
  colorscale = "Greys",
  showscale = FALSE) %>% add_trace(
    df,
    x =  ~ PKA,
    y =  ~ LOGP,
    text = df$NAME,
    type = "scatter",
    mode = "markers",
    marker = list(
      color = ~ MW,
      sizeref = 55,
      size = ~ MW,
      colorscale = colorscale,
      colorbar = list(title = "Molecular<br>Weight"),
      symbol = "circle",
      opacity = 0.7,
      sizemode = "diameter",
      line = list(color = "#444")
    )
  ) %>% layout(
    xaxis = list(
      tickcolor = "rgb(250, 250, 250)",
      title = "pkA",
      titlefont = list(color = "rgb(0,0,0)"),
      zeroline = FALSE,
      showticklabels = FALSE,
      showline = FALSE,
      showgrid = FALSE
    ),
    yaxis = list(
      tickcolor = "rgb(250, 250, 250)",
      title = "LogP",
      titlefont = list(color = "rgb(0,0,0)"),
      zeroline = FALSE,
      showticklabels = FALSE,
      showline = FALSE,
      showgrid = FALSE
    ),
    plot_bgcolor = "rgb(0, 0, 0)",
    showlegend = FALSE
  )

scatterGraph <- plot_ly(
  df,
  x =  ~ PKA,
  y =  ~ LOGP,
  text = df$NAME,
  type = "scatter",
  mode = "markers",
  marker = list(
    color = ~ MW,
    sizeref = 55,
    size = ~ MW,
    colorscale = colorscale,
    colorbar = list(title = "Molecular<br>Weight"),
    symbol = "circle",
    opacity = 0.7,
    sizemode = "diameter",
    line = list(color = "#444")
  )
) %>% layout(
    xaxis = list(
      zeroline = FALSE,
      gridcolor = "rgb(255, 255, 255)",
      showticklabels = FALSE,
      title = "pkA"
    ),
    yaxis = list(
      zeroline = FALSE,
      gridcolor = "rgb(255, 255, 255)",
      showticklabels = FALSE,
      title = "LogP"
    ),
    plot_bgcolor = "#d8d8d8"
  )

# Options for dccDropDown
drugOptions <-  lapply(unique(df$NAME),
  function(drugName) {
    list(label = drugName, value = drugName)
})

startingDrug <- "Levobupivacaine"
# Get idescription from csv, initalized with starting drug
drugDescription <- df$DESC[df$NAME == startingDrug]
# Get img url from csv, initalized with starting drug
drugImg <- df$IMG_URL[df$NAME == startingDrug]

### Helper Function(s) ###
MakeDashTable <- function(selection){
  # Subset df to selected drug name
  dfSubset <- df[is.element(df$NAME, selection), ]
  # Wrap each cell then row content with html tags
  table <- lapply(1:length(selection),
                  function(i) {
                    htmlTr(children = list(
                      htmlTd(toString(dfSubset$NAME[[i]])),
                      htmlTd(toString(dfSubset$FORM[[i]])),
                      htmlTd(htmlTd(htmlImg(
                        src = toString(dfSubset$IMG_URL[[i]])
                      ))),
                      htmlTd(htmlA(
                        href = toString(dfSubset$PAGE[[i]]),
                        children = "Datasheet"
                      ))
                    ))
                  })
  return(table)
}

### App Start ###
app <- Dash$new()

app$layout(
  htmlDiv(list(
    htmlDiv(list(
      #Dash Icon
      htmlImg(src = "https://dash.plot.ly/assets/images/logo.png")
    ), className = "app__banner"),

    htmlDiv(list(
      htmlDiv(list(
        htmlDiv(list(
          # Title + Subtitles
          appTitle,
          subtitle1Bold,
          subtitle1,
          htmlBr(),
          subtitle2Bold,
          subtitle2
        ))
      ), className = "app__header"),
      htmlDiv(list(
        dccDropdown(
          id = "chem_dropdown",
          multi = TRUE,
          value = startingDrug,
          options = drugOptions
        )
      ), className = "app__dropdown"),
      htmlDiv(list(
        htmlDiv(list(
          dccRadioItems(
            id = "charts_radio",
            options = list(list(label = "3D Scatter", value = "scatter3d"),
                           list(label = "2D Scatter", value = "scatter"),
                           list(label = "2D Histogram", value = "histogram2d")),
            labelStyle = list(display = "inline"),
            labelClassName = "radio__labels",
            inputClassName = "radio__input",
            value = "scatter3d",
            className = "radio__group"
          ),
          dccGraph(
            id = "clickable-graph",
            hoverData = list(range = NULL),
            figure = threeD
          )
        ), className = "eight columns"),
        htmlDiv(list(
          htmlDiv(
            htmlImg(
              src = drugImg, id = "chem_img"),
            style = list("text-align" = "center")
            ),
          htmlA(
            children = startingDrug,
            id = "chem_name",
            href = "https://www.drugbank.ca/drugs/DB01002",
            target = "_blank"
          ),
          htmlP(
            children = drugDescription, id = "chem_desc"
          )
        ), className = "chem__desc__container four columns")
      ), className = "container card app__content bg-white row"),
      htmlDiv(list(
        htmlTable(
          id = "table-element",
          children = MakeDashTable(startingDrug),
          className = "table__container"
        )
      ), className = "container bg-white p-0")
    ), className = "app__container")
  ))
)

### Callbacks ###

# Change graph based on radio button
app$callback(
  output = list(id = "clickable-graph", property = "figure"),
  params = list(input(id = "charts_radio", property = "value")),
  function(plotType) {
    switch(
      plotType,
      "scatter3d" = return(threeD),
      "scatter" = return(scatterGraph),
      "histogram2d" = return(histo),
      list()
    )
  }
)

# Description
app$callback(
  output = list(id = "chem_desc", property = "children"),
  params = list(input(id = "clickable-graph", property = "hoverData")),
  function(hover) {
    if (is.null(hover$points[[1]])) {
      return(df$DESC[df$NAME == "Levobupivacaine"])
    } else if (length(hover) == 2) {
      return(drugDescription)
    } else {
      info <- hover$points[[1]]$text
      description <- df$DESC[df$NAME == info]
      if (nchar(as.character(description)) > 0){
        return(description)
      }
    }
  }
)

# Change image on hover
app$callback(
  output = list(id = "chem_img", property = "src"),
  params = list(input(id = "clickable-graph", property = "hoverData")),
  function(hover){
    if (is.null(hover$points[[1]])) {
      return(df$IMG_URL[df$NAME == "Levobupivacaine"])
    } else if (length(hover) == 2) {
      return(drugImg)
    } else {
      info <- hover$points[[1]]$text
      img <- df$IMG_URL[df$NAME == info]
      if (nchar(as.character(img)) > 0){
        return(img)
      }
    }
  }
)

# Chemical link on hover
app$callback(
  output = list(id = "chem_name", property = "href"),
  params = list(input(id = "clickable-graph", property = "hoverData")),
  function(hover) {
    if (is.null(hover$points[[1]])) {
      return(df$PAGE[df$NAME == "Levobupivacaine"])
    } else if (length(hover) == 2) {
      return("https://www.drugbank.ca/drugs/DB01002")
    } else {
      info <- hover$points[[1]]$text
      pageUrl <- df$PAGE[df$NAME == info]
      if (nchar(as.character(pageUrl)) > 0){
        return(pageUrl)
      }
    }
  }
)

# Chemical name on hover
app$callback(
  output = list(id = "chem_name", property = "children"),
  params = list(input(id = "clickable-graph", property = "hoverData")),
  function(hover) {
    if (is.null(hover$points[[1]])) {
      return("Levobupivacaine")
    } else if (length(hover) == 2) {
      return(startingDrug)
    } else {
      info <- hover$points[[1]]$text
      name <- df$NAME[df$NAME == info]
      if (nchar(as.character(name)) > 0){
      return(name)
      }
    }
  }
)

# Append to table
app$callback(
  output = list(id = "table-element", property = "children"),
  params = list(input(id = "chem_dropdown", property = "value")),
  function(chemDropdownValue){
    return(MakeDashTable(chemDropdownValue))
  }
)


if (!appName == "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv("PORT", 8050))
} else {
  app$run_server(debug = TRUE)
}
