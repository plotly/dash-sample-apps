library(dashR)
library(plotly)
library(dashCoreComponents)
library(dashHtmlComponents)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

app <- Dash$new()

df <- read.csv('./data/small_molecule_drugbank.csv', header = TRUE, sep = ",")


###GRAPH PLOTLY OBJECTS###
colorscale <- list(list(0, 'rgb(54,50,153)'), list(0.3, 'rgb(17,123,215)'), list(0.4, 'rgb(37,180,167)'),
                   list(0.5, 'rgb(134,191,118)'), list(0.65, 'rgb(249,210,41)'), list(1, 'rgb(244,236,21'))


title <- htmlH3("dash for drug discovery", className = "uppercase title")
subtitle1_bold <- htmlSpan("Hover ", className = 'uppercase bold')
subtitle1 <-  htmlSpan(
  "over a drug in the graph to see its structure."
)

subtitle2_bold <-htmlSpan("Select ", className = 'uppercase bold')

subtitle2 <- htmlSpan(
  "a drug in the dropdown to add it to the drug candidates at the bottom."
)

three_d <- plot_ly(df, x = ~PKA, z = ~LOGP, y = ~SOL, text = df$NAME, mode = 'markers',
                   marker = list(color = ~MW, 
                                 sizeref = 55,
                                 size = ~MW,
                                 colorscale = colorscale, 
                                 colorbar = list(title = 'Molecular<br>Weight'),
                                 showscale = TRUE, 
                                 symbol = 'circle', 
                                 sizemode='diameter',
                                 line = list( color = '#444')
                   )) %>%
  layout(scene = list(xaxis = list(title = 'pkA'),
                      yaxis = list(title = 'LogP', type = 'log'),
                      zaxis = list(title = 'Solubility (mg/ml)')),
         plot_bgcolor = '#d8d8d8')

histo <- plot_ly(df, x=~PKA, y = ~LOGP) %>% add_histogram2d(colorscale = 'Greys', showscale = FALSE)%>%
  add_trace(df, x=~PKA, y=~LOGP, text = df$NAME,
            type = 'scatter', mode = 'markers',
            marker= list(color = ~MW,
                         sizeref= 55,
                         size = ~MW,
                         colorscale = colorscale,
                         colorbar = list(title = 'Molecular<br>Weight'),
                         symbol = 'circle', 
                         opacity = 0.7,
                         sizemode='diameter', 
                         line = list( color = '#444'))) %>%
  layout(xaxis = list(tickcolor ='rgb(250, 250, 250)',title = 'pkA', titlefont = list(color = 'rgb(0,0,0)'),zeroline=FALSE,  showticklabels = FALSE, showline = FALSE, showgrid = FALSE ),
         yaxis = list(tickcolor ='rgb(250, 250, 250)',title = 'LogP', titlefont = list(color = 'rgb(0,0,0)'),zeroline=FALSE,  showticklabels=FALSE, showline = FALSE, showgrid=FALSE),
         plot_bgcolor = 'rgb(0, 0, 0)',
         showlegend = FALSE
  )

scatter_graph <- plot_ly(df, x=~PKA, y=~LOGP, text = df$NAME,
                         type = 'scatter', mode = 'markers',
                         marker= list(color = ~MW,
                                      sizeref= 55,
                                      size = ~MW,
                                      colorscale = colorscale,
                                      colorbar = list(title = 'Molecular<br>Weight'),
                                      symbol = 'circle', 
                                      opacity = 0.7,
                                      sizemode='diameter', 
                                      line = list( color = '#444'))) %>%
  layout(xaxis = list( zeroline=FALSE, gridcolor = 'rgb(255, 255, 255)', showticklabels = FALSE,title = 'pkA' ),
         yaxis = list( zeroline=FALSE, gridcolor = 'rgb(255, 255, 255)', showticklabels = FALSE,title = 'LogP'),
         plot_bgcolor = '#d8d8d8'
  )

make_dash_table <- function(selection){
  #get subset of dataframe where the only the selected drug names are included
  df_subset <- df[is.element(df$NAME, selection),]
  table <- list()
  for(i in 1:length(selection)){
    
    table <- lapply(1:length(selection), 
                    function(i) { 
                      htmlTr( 
                        children= list( htmlTd(toString(df_subset$NAME[[i]])), 
                                        htmlTd(toString(df_subset$FORM[[i]])), 
                                        htmlTd(htmlTd(htmlImg(src=toString(df_subset$IMG_URL[[i]])))), 
                                        htmlTd(htmlA(href=toString(df_subset$PAGE[[i]]), children = "Datasheet"))) 
                      ) })
  }
  return(table)
}
#for dccDropDown options
drug_names <- unique(df$NAME)
drug_options <- lapply(unique(df$NAME), function(drug_name){ list(label = drug_name, value = drug_name) })

figure <- three_d #initial figure
starting_drug <- 'Levobupivacaine'
#get idescription from csv, initalized with starting drug
drug_description <- df$DESC[df$NAME == starting_drug]
#get img url from csv, initalized with starting drug
drug_img <- df$IMG_URL[df$NAME == starting_drug]

app$layout(
  htmlDiv(list(
    htmlDiv(list(
      #Dash icon
      htmlImg(src="https://dash.plot.ly/assets/images/logo.png")
    ),className = 'app__banner'),
    
    htmlDiv(list(
      htmlDiv(list(
        htmlDiv(list(
          #title + subtitles
          title,
          subtitle1_bold,
          subtitle1,
          htmlBr(),
          subtitle2_bold,
          subtitle2
        ))
      ), className = 'app__header'),
      htmlDiv(list(
        dccDropdown(
          id = 'chem_dropdown',
          multi = TRUE,
          value = starting_drug, 
          options = drug_options
        )
      ), className = "app__dropdown"),
      htmlDiv(list(
        htmlDiv(list(
          dccRadioItems(
            id = 'charts_radio',
            options = list(list(label = '3D Scatter', value = 'scatter3d'),
                           list(label = '2D Scatter', value = 'scatter'),
                           list(label = '2D Histogram', value = 'histogram2d')),
            labelStyle = list(display='inline'), 
            labelClassName = "radio__labels",
            inputClassName = "radio__input",
            value = 'scatter3d',
            className = 'radio__group'
          ),
          dccGraph(
            id = 'clickable-graph',
            hoverData = list(points= list(), range = NULL),
            figure = figure
          )
        ), className= 'two-thirds column'),
        htmlDiv(list(
          htmlImg(
            src = drug_img , id='chem_img'
          ),
          htmlA(
            children=starting_drug, id='chem_name', href = "https://www.drugbank.ca/drugs/DB01002",target='_blank'
          ),
          htmlP(
            children=drug_description, id='chem_desc'
          )
        ),className = "chem__desc__container")
      ), className = "container card app__content bg-white"),
      htmlDiv(list(
        htmlTable(
          id = 'table-element',
          children=make_dash_table(starting_drug),
          className = "table__container"
        )
      ), className = "container bg-white p-0")
    ),className = "app__container")
  ))
)


#####CALLBACKS#####

#change graph due to radio button
app$callback(
  output = list(id = 'clickable-graph', property = 'figure'),
  params = list(input(id = 'charts_radio', property = 'value')),
  display_graph <- function(plot_type){
    switch(plot_type, 'scatter3d' = return(three_d), 
           'scatter' = return(scatter_graph), 
           'histogram2d' = return(histo), list() )
  }
)


#callback for description
app$callback(
  output = list(id = 'chem_desc', property = 'children'),
  params = list(input(id = 'clickable-graph', property = 'hoverData')),
  function(hover){
    if(length(hover)==2){
      return(drug_description)
    }
    info <- hover$points[[1]]$text
    description <- df$DESC[df$NAME == info]
    return(description)
  }
)

#call back from img
app$callback(
  output = list(id = 'chem_img', property = 'src'),
  params = list(input(id = 'clickable-graph', property = 'hoverData')),
  function(hover){
    if(length(hover)==2){
      return(drug_img)
    }
      info <- hover$points[[1]]$text
      img <- df$IMG_URL[df$NAME == info]
      return(img)
  }
)

#callback for name URL
app$callback(
  output = list(id = 'chem_name', property = 'href'),
  params = list(input(id = 'clickable-graph', property = 'hoverData')),
  function(hover){
    if(length(hover)==2){
      return("https://www.drugbank.ca/drugs/DB01002")
    }
      info <- hover$points[[1]]$text
      page_url <- df$PAGE[df$NAME==info]
   
    return(page_url)
  }
)

#callback for name
app$callback(
  output = list(id = 'chem_name', property = 'children'),
  params = list(input(id = 'clickable-graph', property = 'hoverData')),
  function(hover){
    if(length(hover)==2){
      return(starting_drug)
    }
      info <- hover$points[[1]]$text
      name <- df$NAME[df$NAME==info]
    return(name)
  }
)

#callback to append to table
app$callback(
  output = list(id = 'table-element', property = 'children'),
  params = list(input(id = 'chem_dropdown', property ='value')),
  function(chem_dropdown_value){
    table1 <- make_dash_table(chem_dropdown_value)
    return(table1)
  }
)

app$run_server(showcase = TRUE)

