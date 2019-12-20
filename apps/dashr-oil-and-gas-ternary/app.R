library(dashHtmlComponents)
library(dashCoreComponents)
library(dash)
library(plotly)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

source('constants.R')
source('dataprep.R')

app = Dash$new()

mapbox_access_token = 'pk.eyJ1IjoiamFja2x1byIsImEiOiJjajNlcnh3MzEwMHZtMzNueGw3NWw5ZXF5In0.fk8k06T96Ml9CLGgKmk81w'
Sys.setenv('MAPBOX_TOKEN' = mapbox_access_token)
df = read.csv("data/test_compositionR.csv")
df_prod = read.csv("data/YearlyProduction_table_1.csv")

colormap <- setNames(colors[1:length(unique(df$fm_name))], unique(df$fm_name))


build_banner = function(){
  return(htmlDiv(
    id='banner',
    className='banner',
    children=list(
      htmlH6("Oil and gas ternary map"),
      htmlImg(
        src="assets/dash-logo.png"
      )
    )
  ))}

build_graph_title = function(title){
  return(htmlP(
    className='graph-title',
    children=title))}

generate_production_plot = function(processed_data){
  layout = list(
    xaxis=list(
      title='Year'
    ),
    yaxis=list(
      title='GAS Production(mcf)'
    )
  )
  
  data = list()
  
  welllist = lapply(1:length(processed_data[['well_id']]), 
                    function(x){list('well_id' = unlist(processed_data[['well_id']][x]), 
                                     'formation' = unlist(processed_data[['formation']][x]))})
  
  for (well_id in welllist) {
    well_prod = df_prod[df_prod$RecordNumber %in% well_id$well_id,]
    new_trace = list(
      x=well_prod$Year,
      y=well_prod$VolumeMCF,
      name=c(well_id$well_id),
      mode='lines+markers',
      hoverinfo='x+y+name',
      marker=list(symbol="hexagram-open", 
                  color=colormap[[as.character(well_id$formation)]],
                  line=list("width" = "0.5")),
      line=list(
        shape='spline'
      ),
      showlegend=TRUE
    )
    data = append(data,list(new_trace))
  }
  return(list('data' = data, 'layout' = layout))
  
}

generate_well_map = function(dff,selected_data,style){
  
  formations = as.list(unique(dff$fm_name))
  
  data = list()
  
  for (formation in unlist(formations)) {
    selected_index = NULL
    if(formation %in% labels(selected_data)){
      selected_index = selected_data[[formation]]
    }
    
    text_list = list(
      paste('Well ID', dff[dff$fm_name == formation,'RecordNumber']))
    
    op_list = as.list(dff[dff$fm_name == formation,'op'])
    
    text_list = lapply(1:length(text_list[[1]]), function(x){paste(op_list[[x]], "<br>", text_list[[1]][x])})
    
    
    new_trace = list(
      list(
        type="scattermapbox",
        lat = dff[dff$fm_name == formation,'nlat83'],
        lon = dff[dff$fm_name == formation,'wlon83'],
        mode = 'markers',
        marker = list('color' = colormap[[formation]], 'size' = 9),
        text= unlist(text_list),
        name = c(unlist(formation)),
        selectedpoints = selected_index,
        customdata = dff[dff$fm_name == formation,'RecordNumber']
      )
    )
    data = append(data, new_trace)
  }
  
  layout = list(
    clickmode="event+select",
    dragmode="lasso",
    showlegend=TRUE,
    autosize=TRUE,
    hovermode="closest",
    margin= list(l=0, r=0, t=0, b=0),
    mapbox = list(
      accesstoken=mapbox_access_token,
      bearing=0,
      center= list(lat=37.497562, lon=-82.755728),
      pitch=0,
      zoom=8,
      style=style
    ),
    legend=list(
      bgcolor="#1f2c56",
      orientation="h",
      font=list(color="white"),
      x=0,
      y=0,
      yanchor="bottom"
    )
  )
  
  return(list(data = data, layout = layout))
  
}



generate_ternary_map = function(dff, selected_data, contour_visible, marker_visible){
  contour_traces = list()
  j = 1
  for (key in labels(ternary_contour)) {
    contour_traces[[j]] = list(
      name=key,
      type='scatterternary',
      a= unlist(lapply(ternary_contour[[key]], function(k){k[1]})),
      b = unlist(lapply(ternary_contour[[key]], function(k){k[2]})),
      c = unlist(lapply(ternary_contour[[key]], function(k){k[3]})),
      mode='lines',
      line=list(color='#444', width=0.5),
      fill='toself',
      fillcolor=ternary_color[j][[1]],
      opacity=0.8,
      hoverinfo='all',
      showlegend=FALSE,
      visible=contour_visible
    )
    j = j+1
  }
  
  contour_text = generate_contour_text_layer(contour_visible)
  
  layout = list(
    'height' = 600,
    'margin' = list(l=50, r=10, b=40, t=40, pad=5),
    'ternary' = list('sum' = 100, 'aaxis' = list(
      'title'= list('text'= 'Quartz',
                    'font' = list('family'= 'Open Sans', 'size'= 15, 'color'= "white")), 
      'min'= -2,'linewidth'= 1.5,'ticks'= 'outside'),
      'baxis'= list('title'= list('text'= 'Carbonate', 'font'= list('family'= 'Open Sans', 'size'= 15, 
                                                                    'color'= "white")), 'min'= -2, 'linewidth'= 1.5,'ticks'= 'outside'),
      'caxis'= list('title'= list('text'= 'Clay', 'font'= list('family'= 'Open Sans', 'size'= 15, 'color'= "white")),
                    'min'= -2, 'linewidth'= 1.5, 'ticks'= 'outside')),
    "margin" = list(l=60, r=0, t=0, b=0),
    "paper_bgcolor" = "#192444",
    "plot_bgcolor" = "#192444",
    "showLegend" = FALSE,
    "font" = list("color" = "white"),
    "annotations" = list("visible" = FALSE),
    "autosize" = TRUE
  ) 
  hovertemplate = "<b> %{text}</b><br><br> Quartz: %{a:.0f}<br>Carbonate : 
  %{b:.0f}<br> Clay: %{c:.0f}<extra></extra>"
  
  formations = list(unique(dff['fm_name']))
  data_traces = list()
  j = 1
  for (key in formations[[1]]$fm_name) {
    if(length(selected_data) > 0){
      select_indices = selected_data[[key]]
    } else{
      select_indices = NULL
    }
    data_traces[[j]] = list(
      text=as.list(paste('Well ID', 
                         dff[dff$fm_name == formation,'RecordNumber'])),
      name=key,
      customdata=dff[dff$fm_name == key,'RecordNumber'],
      type='scatterternary',
      a=dff[dff$fm_name == key,'Quartz'],
      b=dff[dff$fm_name == key,'Carbonate'],
      c=dff[dff$fm_name == key,'Clay'],
      mode='markers',
      hovertemplate = hovertemplate,
      showlegend=TRUE,
      marker=list('color' = colormap[[key]], 'size' = 8, 'line' = list('color' = '#000000', 'width'= 0.2)),
      selectedpoints=select_indices,
      visible=marker_visible
    )
    j = j+1
  }
  FData = append(contour_traces, contour_text)
  FData = append(FData, data_traces)
  return(list('data' = FData, 'layout' = layout))
}

generate_contour_text_layer = function(contour_visible){
  layer = list()
  j = 1
  for (key in labels(ternary_contour)) {
    a= mean(unlist(lapply(ternary_contour[[key]], function(k){k[1]})))
    b = mean(unlist(lapply(ternary_contour[[key]], function(k){k[2]})))
    c = mean(unlist(lapply(ternary_contour[[key]], function(k){k[3]})))
    key_br = gsub(" ", '<br>', key)
    new_trace = generate_contour_text(a, b, c, key, key_br, contour_visible)
    layer[[j]] = new_trace
    j = j+1
  }
  return(layer)
}

generate_contour_text = function(a, b, c, name, text, visible){
  return(list(
    type='scatterternary',
    a=list(a),
    b=list(b),
    c=list(c),
    name=name,
    text=list(text),
    mode='text',
    hoverinfo='all',
    textposition='middle center',
    textfont=list('size' = 11, 'color' = '#000000', 'family' = 'sans-serif'),
    showlegend=FALSE,
    legendgroup='Rock type',
    visible=visible
  ))
}

generate_formation_bar = function(dff, selected_data){
  
  formations = list(unique(dff['fm_name']))
  xform <- list(categoryorder = "array", 
                categoryarray = formations[[1]]$fm_name,
                tickangle=-45, title="Formations")
  data = list()
  countt = list()
  color = list()
  
  for (i in formations[[1]]$fm_name) {
    counttt = nrow(dff[dff$fm_name == i,])
    countt = append(countt, counttt)
    colorr =unname(colormap[[i]])
    color = append(color,colorr)
  }
  if(length(selected_data) > 0){
    selected_points = list()
    select = list()
    for (i in formations[[1]]$fm_name) {
      select_indices = selected_data[[i]]
      if(length(select_indices) > 0 && is.null(select_indices) == FALSE){
        select = list(which(names(selected_data) == i))
      }
      selected_points = append(selected_points, select)
    }
    new_trace = plot_ly(
      type = 'bar',
      x = formations[[1]]$fm_name,
      y = countt,
      name = i,
      hoverinfo = 'x+y',
      marker = list(color = c(unlist(color))),
      selectedpoints = selected_points
    ) %>% layout(showlegend=FALSE,
                 hovermode="closest",
                 xaxis=xform,
                 yaxis=list(title="Well Counts"),
                 clickmode="event+select")
  }
  else{
    new_trace = plot_ly(
      type = 'bar',
      x = as.vector(as.character(formations[[1]]$fm_name)),
      y = unlist(countt),
      name = i,
      marker = list(color = c(unlist(color))),
      selectedpoints = NULL
    )%>% layout(showlegend=FALSE,
                hovermode="closest",
                xaxis=xform,
                yaxis=list(title="Well Counts"),
                clickmode="event+select")
  }
  return(new_trace)
}




get_selection = function(OurData,formation,selection_data,starting_index){
  ind = list()
  current_curve = which(sapply(as.list(unique(OurData$fm_name)), 
                               FUN=function(X){formation %in% X})) 
  for (point in selection_data[['points']]) {
    if(point[['curveNumber']] - starting_index == current_curve){
      ind = append(ind, point[['pointNumber']])
    } 
  }
  return(ind)
}

get_selectionFormation = function(OurData,formation,selection_data,starting_index){
  ind = list()
  current_curve = which(sapply(as.list(unique(OurData$fm_name)), 
                               FUN=function(X){formation %in% X})) -1
  for (point in selection_data[['points']]) {
    if(point[['curveNumber']] - starting_index == current_curve){
      ind = append(ind, point[['pointNumber']])
    } 
  }
  return(ind)
}


get_selection_by_bar = function(bar_selected_data){
  dict = list()
  if(length(bar_selected_data) > 0){
    for (point in bar_selected_data[['points']]) {
      if(length(point$x) > 0){
        dict[[(point$x)]] = as.list(seq(0, point$y))
      }
    }
  }
  return(dict)
}

app$layout(htmlDiv(
  children=list(
    htmlDiv(
      id="top-row",
      children=list(
        htmlDiv(
          className="row",
          id="top-row-header",
          children=list(
            htmlDiv(
              id="header-container",
              children=list(
                build_banner(),
                htmlP(
                  id="instructions",
                  children="Select data points from the well map, ternary map or bar graph to
                  visualize cross-filtering to other plots. Selection could be done by 
                  clicking on individual data points or using the lasso tool to capture 
                  multiple data points or bars. With the box tool from modebar, multiple 
                  regions can be selected by holding the SHIFT key while clicking and 
                  dragging."
                ),
                build_graph_title("Select Operator"),
                dccDropdown(
                  id="operator-select",
                  options= lapply(unique(df['op'])$op, function(x){list('label' = x, 'value' = x)}),
                  multi=TRUE,
                  value=list(
                    list(unique(df$op))[[1]][1],
                    list(unique(df$op))[[1]][2]
                  )
                )
                )
              )
          )
        ),
        htmlDiv(
          className="row",
          id="top-row-graphs",
          children=list(
            # Well map
            htmlDiv(
              id="well-map-container",
              children=list(
                build_graph_title("Well Map"),
                dccRadioItems(
                  id="mapbox-view-selector",
                  options=list(
                    list("label"= "basic", "value"= "basic"),
                    list("label"= "satellite", "value"= "satellite"),
                    list("label"= "outdoors", "value"= "outdoors"),
                    list(
                      "label"= "satellite-street",
                      "value"= "mapbox://styles/mapbox/satellite-streets-v9"
                    )
                  ),
                  value="basic"
                ),
                dccGraph(
                  id="well-map",
                  figure=list(
                    "layout"= list(
                      "paper_bgcolor"= "#192444",
                      "plot_bgcolor"= "#192444"
                    )
                  ),
                  config=list("scrollZoom"= TRUE, "displayModeBar"= TRUE)
                )
              )
            ),
            # Ternary map
            htmlDiv(
              id="ternary-map-container",
              children=list(
                htmlDiv(
                  id="ternary-header",
                  children=list(
                    build_graph_title(
                      "Shale Mineralogy Composition"
                    ),
                    dccChecklist(
                      id="ternary-layer-select",
                      options=list(
                        list(
                          "label"= "Well Data",
                          "value"= "Well Data"
                        ),
                        list(
                          "label"= "Rock Type",
                          "value"= "Rock Type"
                        )
                      ),
                      value=list("Well Data", "Rock Type")
                    )
                  )
                ),
                dccGraph(
                  id="ternary-map",
                  figure=list(
                    "layout"= list(
                      "paper_bgcolor"= "#192444",
                      "plot_bgcolor"= "#192444"
                    )
                  ),
                  config=list(
                    "scrollZoom"= TRUE,
                    "displayModeBar"= TRUE
                  )
                )
              )
            )
          )
        )
      )
  ),
  htmlDiv(
    className="row",
    id="bottom-row",
    children=list(
      # Formation bar plots
      htmlDiv(
        id="form-bar-container",
        className="six columns",
        children=list(
          build_graph_title("Well count by formations"),
          dccGraph(id="form-by-bar")
        )
      ),
      htmlDiv(
        # Selected well productions
        id="well-production-container",
        className="six columns",
        children=list(
          build_graph_title("Individual well annual production"),
          dccGraph(id="production-fig")
        )
      )
    )
  )
)
))
app$callback(output = list(id = "form-by-bar", property = "figure"),
             params = list(
               input("well-map", "selectedData"),
               input("ternary-map", "selectedData"),
               input("operator-select", "value")),
             function(map_selected_data, tern_selected_data, op_select){
               dff = df[df$op %in% op_select,]
               formations = list(unique(dff['fm_name']))
               ctx = app$callback_context()
               prop_id = ""
               prop_type = ""
               if (!is.null(ctx$triggered) == TRUE){
                 splitted = as.list(strsplit(ctx$triggered['prop_id']$prop_id, ".",fixed=TRUE)[1][[1]])
                 prop_id = splitted[1]
                 prop_type = splitted[2]
               }
               processed_data = list()
               if(prop_id == 'well-map' & prop_type == 'selectedData'){
                 for (formation in formations[[1]]$fm_name) {
                   if(length(map_selected_data) == 0){
                     processed_data[[formation]] = list(0)
                   } else{
                     processed_data[[formation]] = get_selection(dff, formation, map_selected_data, 0)
                   }
                 }
               }
               else if (prop_id == 'ternary-map' & prop_type == 'selectedData'){
                 for (formation in formations[[1]]$fm_name) {
                   if (length(tern_selected_data) == 0){
                     processed_data[[formation]] = list(0)
                   }
                   else{
                     processed_data[[formation]] = get_selection(dff, formation, tern_selected_data, 32)
                   }
                 }
               }
               else{
                 processed_data = processed_data
               }
               return(generate_formation_bar(dff, processed_data))
             })

app$callback(output = list(id = "ternary-map", property = "figure"),
             params = list(
               input("well-map", "selectedData"),
               input("form-by-bar", "selectedData"),
               input("form-by-bar", "clickData"),
               input("operator-select", "value"),
               input("ternary-layer-select", "value"),
               state("ternary-map", "figure")),
             function(map_selected_data, bar_selected_data, bar_click_data, op_select,
                      layer_select, curr_fig){
               marker_visible = TRUE
               contour_visible = TRUE
               dff = df[df$op %in% op_select,]
               formations = list(unique(dff['fm_name']))
               ctx = app$callback_context()
               prop_id = ""
               prop_type = ""
               if(isTRUE(ctx$triggered) == FALSE){
                 splitted = as.list(strsplit(ctx$triggered['prop_id']$prop_id, ".",fixed=TRUE)[1][[1]])
                 prop_id = splitted[1]
                 prop_type = splitted[2]
               }
               processed_data = list()
               
               if(prop_id != 'ternary-layer-select'){
                 if(prop_id == "well-map" & prop_type == "selectedData"){
                   for (formation in formations[[1]]$fm_name) {
                     if(length(map_selected_data) == 0){
                       processed_data[[formation]] = NULL
                     } else{
                       processed_data[[formation]] = get_selectionFormation(dff, formation, map_selected_data, 0)
                     }
                   }
                 }
                 else if(prop_id == "form-by-bar" & prop_type == "selectedData"){
                   processed_data = get_selection_by_bar(bar_selected_data)
                   for (formation in formations[[1]]$fm_name) {
                     if(length(bar_selected_data) == 0){
                       processed_data[[formation]] = NULL
                     } else if(!(formation %in% labels(processed_data))){
                       processed_data[[formation]] = list()
                     }
                   }
                   
                 } else if(prop_id == "form-by-bar" & prop_type == "clickData"){
                   processed_data = get_selection_by_bar(bar_click_data)
                   for (formation in formations[[1]]$fm_name) {
                     if(length(bar_click_data) == 0){
                       processed_data[[formation]] = NULL
                     } else if(!(formation %in% labels(processed_data))){
                       processed_data[[formation]] = list()
                     }
                   }
                 } else{
                   for (formation in formations[[1]]$fm_name) {
                     processed_data[formation] = NULL
                   }
                 }
                 return(generate_ternary_map(
                   dff, processed_data, contour_visible, marker_visible))
               }
               
               if (prop_id == "ternary-layer-select"){
                 if(!(length(curr_fig) == 0)){
                   if(!("Well Data" %in% layer_select)){
                     marker_visible = "legendonly"
                   }
                   if (!("Rock Type" %in% layer_select)){
                     contour_visible = "legendonly"
                   }
                   curr_fig = generate_ternary_map(dff, processed_data, 
                                                   contour_visible, marker_visible)
                   return(curr_fig)
                 } else{
                   return(curr_fig)
                 }
               }
             })

app$callback(output = list(id = "well-map", property = "figure"),
             params = list(input("ternary-map", "selectedData"),
                           input("form-by-bar", "selectedData"),
                           input("form-by-bar", "clickData"),
                           input("operator-select", "value"),
                           input("mapbox-view-selector", "value")),
             function(tern_selected_data, bar_selected_data, bar_click_data, op_select, mapbox_view){
               dff = df[df$op %in% op_select,]
               formations = list(unique(dff['fm_name']))
               ctx = app$callback_context()
               prop_id = ""
               prop_type = ""
               
               if(isTRUE(ctx$triggered) == FALSE){
                 splitted = as.list(strsplit(ctx$triggered['prop_id']$prop_id, ".",fixed=TRUE)[1][[1]])
                 prop_id = splitted[1]
                 prop_type = splitted[2]
               }
               
               processed_data = list()
               if (prop_id == 'ternary-map'){
                 for (formation in formations[[1]]$fm_name) {
                   if(length(tern_selected_data) == 0){
                     processed_data[[formation]] = NULL
                   } 
                   else{
                     processed_data[[formation]] = get_selection(dff, formation, tern_selected_data, 32)
                   }
                 }
               } 
               else if(prop_id == 'form-by-bar'){
                 bar_data = ""
                 if(prop_type == 'selectedData'){
                   bar_data = bar_selected_data
                 } else if(prop_type == "clickData"){
                   bar_data = bar_click_data
                 }
                 processed_data = get_selection_by_bar(bar_data)
                 
                 for (formation in formations[[1]]$fm_name) {
                   if(length(bar_data) == 0){
                     processed_data[[formation]] = NULL
                   } else if(!(formation %in% labels(processed_data))){
                     processed_data[[formation]] = list()
                   }
                 }
               } else{
                 for (formation in formations[[1]]$fm_name) {
                   processed_data[[formation]] = NULL
                 }}
               
               return(generate_well_map(dff, processed_data, mapbox_view))
             })

app$callback(output = list(id = "production-fig", property = "figure"),
             params = list(input("well-map", "selectedData"),
                           input("ternary-map", "selectedData"),
                           input("form-by-bar", "selectedData"),
                           input("operator-select", "value")),
             function(map_select, tern_select, bar_select, op_select){
               dff = df[df$op %in% op_select,]
               ctx = app$callback_context()
               prop_id = ""
               prop_type = ""
               
               if(isTRUE(ctx$triggered) == FALSE){
                 splitted = as.list(strsplit(ctx$triggered['prop_id']$prop_id, ".",fixed=TRUE)[1][[1]])
                 prop_id = splitted[1]
                 prop_type = splitted[2]
               }
               processed_data_init = list()
               processed_data_init[["well_id"]] = as.list(dff$RecordNumber)
               processed_data_init[["formation"]] = as.list(dff$fm_name)
               
               if(prop_id == "well-map" & prop_type == "selectedData"){
                 if(length(map_select) > 0){
                   processed_data = list('well_id' = list(), 'formation' = list())
                   for (point in map_select[['points']]) {
                     processed_data[["well_id"]] = append(processed_data[["well_id"]], point[['customdata']])
                     processed_data[['formation']] = append(processed_data[['formation']],
                                                            list(dff[dff$RecordNumber %in% point[['customdata']], 
                                                                     'fm_name']))
                   }
                 } else{
                   processed_data = processed_data_init
                 }
               } else if(prop_id == "ternary-map" & prop_type == "selectedData"){
                 if(length(tern_select) > 0){
                   processed_data = list('well_id' = list(), 'formation' = list())
                   for (point in tern_select[['points']]) {
                     processed_data[["well_id"]] = append(processed_data[["well_id"]], point$customdata)
                     processed_data[['formation']] = append(processed_data[['formation']],
                                                            list(dff[dff$RecordNumber == point$customdata, 
                                                                     'fm_name']))
                   }
                 } 
                 
                 else{
                   processed_data = processed_data_init
                 }
               } else if(prop_id == "form-by-bar" & prop_type == "selectedData"){
                 if(length(bar_select) > 0){
                   processed_data = list('well_id' = list(), 'formation' = list())
                   for (point in bar_select[['points']]) {
                     selected_form = point[['x']]
                     selected_well = list(dff[dff["fm_name"] == point[["x"]],"RecordNumber"])
                     for (well in selected_well[[1]]) {
                       processed_data[["well_id"]] = append(processed_data[["well_id"]],well)
                       processed_data[["formation"]] = append(processed_data[["formation"]],selected_form)
                     }
                   }
                 } else{
                   processed_data = processed_data_init
                 }
               } else{
                 processed_data = processed_data_init
               }
               return(generate_production_plot(processed_data))
             })


if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
