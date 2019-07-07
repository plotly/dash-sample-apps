library(dashHtmlComponents)
library(dashCoreComponents)
library(dash)

source('constants.R')
source('dataprep.R')

app <- Dash$new()

mapbox_access_token <- 'pk.eyJ1IjoieWNhb2tyaXMiLCJhIjoiY2p1MDR5c3JmMzJsbjQ1cGlhNHA3MHFkaCJ9.xb3lXp5JskCYFewsv5uU1w'

df <- read.csv("data/test_compositionR.csv")
df_prod <- read.csv("data/YearlyProduction_table_1.csv")

colormap <- list()
j <- 1
for (i in unique(df$fm_name)) {
  colormap[[i]] = colors[[j]]
  j <- j + 1
}

build_banner <- function(){
  return(htmlDiv(
    id = 'banner',
    className = 'banner',
    children = list(
      htmlH6("Oil and gas ternary map"),
      htmlImg(
        src = "assets/logo.png"
      )
    )
  ))}

build_graph_title <- function(title){
  return(htmlP(
    className='graph-title',
    children=title))}

generate_production_plot <- function(processed_data){
  layout = list(
    xaxis=list(
      title='Year'
    ),
    yaxis=list(
      title='GAS Production(mcf)'
    )
  )
  data = list()
  
  for (well_id in list(processed_data['well_id'],processed_data[formation])) {
    well_prod = df_prod[df_prod['RecordNumber'] == well_id]
    new_trace = list(
      x=well_prod['Year'],
      y=well_prod['VolumeMCF'],
      name=str(well_id),
      mode='lines+markers',
      hoverinfo='x+y+name',
      marker=list(symbol="hexagram-open", color=colormap[formation]),
      line=list(
        shape='spline'
      ),
      showlegend=TRUE
    )
  }
  data = append(data,new_trace)
  return(list('data' = data, 'layout' = layout))
  
  
  }

generate_well_map <- function(dff, selected_data, style){
  Layout= layout(annotations=list(
    list(
      x=1.1,
      y=0.85,
      align="right",
      valign="top",
      text='<b>Formations</b>',
      showarrow=False,
      xref="paper",
      yref="paper",
      font=list(
        family='sans-serif',
        size=13,
        color='#000'
      )
    )
    ),
    legend=list(
      x=1,
      y=0.8,
      traceorder='normal',
      font=list(
        family='sans-serif',
        size=12,
        color='#000'
      ),
      bordercolor='#FFFFFF',
      borderwidth=2),
    
    margin=list('l'= 5, 'r' = 5, 't' = 20, 'b' = 5, 'pad' = 0),
    height=550,
    clickmode='event+select',
    showlegend=TRUE,
    autosize=TRUE,
    hovermode='closest',
    mapbox=list(
      accesstoken=mapbox_access_token,
      bearing=0,
      center=list(
        lat=37.497562,
        lon=-82.755728
      ),
      pitch=0,
      zoom=8,
      style=style
    ))
  
  formations <- list(unique(dff['fm_name']))
  
  data <- list()
  
  for (formation in formations) {
    selected_index = NULL
    if(formation %in% selected_data){
      selected_index = selected_dta[formation]
    }
  }
  new_trace <- plot_mapbox(
    lat <- dff[dff['fm_name'] == formation,'nlat83'],
    lon <- dff[dff['fm_name'] == formation,'wlon83'],
    mode <- 'markers',
    marker <- list('color' = colormap[formation], 'size' = 9),
    #text = list(map(lambda item= 'Well ID=' + str(int(item)), dff[dff['fm_name'] == formation]['RecordNumber'])),
    name <- formation,
    selectedpoints = selected_index,
    customdata = dff[dff['fm_name'] == formation,'RecordNumber']
  ) %>% Layout
  data <- append(data,new_trace)
  
  return(list('data' = data))
  
}


generate_ternary_map <- function(dff, selected_data, contour_visible, marker_visible){
  k <- 1
  contour_traces <- list()
  for (key in labels(ternary_contour)) {
    for (k in ternary_contour[[key]]){
    trace = list(
      name = key,
      type = 'scatterternary',
      a = list(k['Quartz']),
      b = list(k['Carbonate']),
      c = list(k['Clay']),
      mode = 'lines',
      line = list(color='#444', width=0.5),
      fill = 'toself',
      fillcolor = ternary_color[[k]],
      opacity = 0.8,
      hoverinfo = 'none',
      showlegend = FALSE,
      visible = contour_visible
    )}
    k = k+1
    contour_traces = append(contour_traces,trace)
  }
  
  contour_text <- generate_contour_text_layer(contour_visible)

  layout <- list(
    'height' = 600,
    'margin' = dict(l=50, r=10, b=40, t=40, pad=5),
    'ternary' =
    list('sum' = 100,
      'aaxis' = list('title'= list('text'= 'Quartz',
        'font' = list('family'= 'Open Sans', 'size'= 15, 'color'= '#000')),
        'min'= -2,
        'linewidth'= 1.5,
        'ticks'= 'outside'),
      'baxis'= list('title'= list('text'= 'Carbonate',
        'font'= list('family'= 'Open Sans', 'size'= 15, 'color'= '#000')),
        'min'= -2,
        'linewidth'= 1.5,
        'ticks'= 'outside'),
      'caxis'= list('title'= list('text'= 'Clay',
        'font'= list('family'= 'Open Sans', 'size'= 15, 'color'= '#000')),
        'min'= -2,
        'linewidth'= 1.5,
        'ticks'= 'outside')),
    'legend'=
      list(
        x=1,
        y=0.8,
        traceorder='normal',
        font=list(
          family='sans-serif',
          size=12,
          color='#000'
        ),
        bordercolor='#FFFFFF',
        borderwidth=2)
  )
  #hovertemplate = "<b> %{text}</b><br><br> Quartz: %{a:.0f}<br>Carbonate : %{b:.0f}<br> Clay: %{c:.0f}<extra></extra>
  formations <- list(unique(dff['fm_name']))
  data_traces <- list()
  for (key in formations) {
    if(selected_data){
      select_indices = selected_data[key]
    } else{
      select_indices = NULL
    }
    new_data_trace <- list(
      #text=list(map(lambda item: 'Well ID:' + str(int(item)), dff[dff['fm_name'] == key]['RecordNumber'])),
      name = key,
      customdata = dff[dff['fm_name'] == key,'RecordNumber'],
      type = 'scatterternary',
      a = dff[dff['fm_name'] == key,'Quartz'],
      b = dff[dff['fm_name'] == key,'Carbonate'],
      c = dff[dff['fm_name'] == key,'Clay'],
      mode = 'markers',
      hovertemplate = hovertemplate,
      showlegend = TRUE,
      marker = list('color' = colormap[key], 'size' = 8, 'line' = list('color' = '#000000', 'width': 0.2)),
      selectedpoints = select_indices,
      visible=marker_visible
    )
    data_traces = append(data_traces,new_data_trace)
  }
  return(list('data' = contour_traces + contour_text + data_traces, 'layout' = layout))
}

generate_contour_text_layer <- function(contur_visible){
  layer <- list()
  for (key in labels(ternary_contour)) {
    for (k in ternary_contour[[key]]) {
      a <- mean(unlist(k['Quartz']))
      b <- mean(unlist(k['Carbonate']))
      c <- mean(unlist(k['Clay']))
    }
    key_br <- gsub(key,' ', '<br>')
    new_trace <- generate_contour_text(a, b, c, key, key_br, contour_visible)
    layer <- append(layer,new_trace)
  }
  return(layer)
}

generate_contour_text <- function(a, b, c, name, text, visible){
  return(list(
    type='scatterternary',
    a = list(a),
    b = list(b),
    c = list(c),
    name = name,
    text = text,
    mode = 'text',
    hoverinfo = 'none',
    textposition = 'middle center',
    textfont = list('size' = 11, 'color' = '#000000', 'family' = 'sans-serif'),
    showlegend = FALSE,
    legendgroup = 'Rock type',
    visible = visible
  ))
}

generate_formation_bar = function(dff, selected_Data){
  
  layout <- layout(showlegend=FALSE,
                   hovermode="closest",
                   xaxis=list(tickangle=-45, title="Formations"),
                   yaxis=list(title="Well Counts"),
                   clickmode="event+select")
  
  formations <- list(unique(dff['fm_name']))
  
  if(length(selected_Data > 0)){
    data <- list()
    for (i in formations) {
      selected_points <- list()
      selected_indices <- selected_data[i]
      if(length(select_indices) > 0){
        selected_points <- list(0)
      }
      new_trace <- plot_ly(
        type = 'Bar',
        x = list(i),
        y = length(dff[dff['fm_name'] == i,]),
        name = i,
        hoverinfo = 'x+y',
        marker = list('color' = colormap[i]),
        selectedpoints = selected_points
      )
      data <- append(data, new_trace)
    }
  }
  else{
    data <- list()
    for (i in formations) {
      new_trace <- plot_ly(
        type = 'Bar',
        x = list(i),
        y = length(dff[dff['fm_name'] == i,]),
        name = i,
        marker = list('color' = colormap[i]),
        selectedpoints = NULL
      )
      data <- append(data, new_trace)
    }
  }
  return(list('data' = data, 'layout' = layout))
}


get_selection <- function(data,formation,selection_data,starting_index){
  ind <- list()
  current_curve <- list(unique(data['fm_name', formation]))
  for (point in selection_data['points']) {
    if(point['curveNumber'] - starting_index == current_curve){
      ind <- append(ind, point['pointNumber'])
    }
  }
  return(ind)
}


get_selection_by_bar <- function(bar_selected_data) {
 dict <- list()
 if(length(bar_selected_data) > 0){
   for (point in bar_selected_data['points']) {
     if(length(point['x']) > 0){
       dict[(point["x"])] <- list(seq(0, point["y"]))
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
              className="six columns",
              id="header-container",
              children=list(
                build_banner(),
                htmlP(
                  id="instructions",
                  children="Select data points from the well map, ternary map or bar graph to visualize cross-filtering to
                  other plots. Selection could be done by clicking on individual data points or using the lasso tool 
                  to capture multiple data points or bars. With the box tool from modebar, multiple regions can be selected by holding SHIFT
                  button while click-and-drag."
                ),
                build_graph_title("Select Operator"),
                dccDropdown(
                  id="operator-select",
                  # options=list(
                  #  list("label"= i, "value"= i)
                  # for i in df["op"].unique().tolist()
                  #),
                  multi=TRUE,
                  value=list(
                    list(unique(df))["op"][0],
                    list(unique(df["op"]))[1]
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
            className="six columns",
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
                    "value"= "mapbox=//styles/mapbox/satellite-streets-v9"
                  )
                ),
                value="basic"
              ),
              dccGraph(
                id="well-map",
                config=list("scrollZoom"= TRUE, "displayModeBar"= TRUE)
              )
            )
          ),
          # Ternary map
          htmlDiv(
            id="ternary-map-container",
            className="six columns",
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
                config=list("scrollZoom"= TRUE, "displayModeBar"= TRUE)
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
               dff <- df[df['op'] %in% op_select,]
               
               formation <- list(unique(dff['fm_name']))
               ctx <- app$callback_context()
               prop_id <- ""
               prop_type <- ""
               if (isTRUE(ctx$triggered)){
                 splitted <- strsplit(ctx$triggered['prop_id'][0], ".",fixed=TRUE)
                 prop_id <- splitted[0]
                 prop_type <- splitted[1]
               }
               processed_data <- list()
               if(prop_id == 'well-map' & prop_type == 'selectedData'){
                 for (formation in formations) {
                   if(length(map_selected_data) == 0){
                     processed_data[formation] <- list(0)
                   } else{
                     processed_data[formation] <- get_selection(dff, formation, map_selected_data, 0)
                   }
                 }
               }
               else if (prop_id == 'ternary-map' & prop_type == 'selectedData'){
                 for (formation in formations) {
                   if (length(tern_selected_data) == 0){
                     processed_data[formation] = list(0)
                   }
                   else{
                     processed_data[formation] = get_selection(dff, formation, tern_selected_data, 32)
                   }
                 }
               }
              else{
                for (formation in formations) {
                  processed_data[formation] = list(0)
                }
              }
               return(generate_formation_bar(dff, processed_data))
             })


app$callback(output = list(id = "ternary-map", property = "figure"),
             params = list(
               input("well-map", "selectedData"),
               input("ternary-map", "selectedData"),
               input("form-by-bar", "clickData"),
               input("operator-select", "value"),
               input("ternary-layer-select", "values"),
               state("ternary-map", "figure")),
             function(map_selected_data,bar_selected_data,bar_click_data,op_select,layer_select, curr_fig){
               marker_visible = countour_visible
               marker_visible = TRUE
               dff = df[df['op'] %in% op_select,]
               formaions = list(unique(dff['fm_name']))
               
               ctx = app$callback()
               if(isTRUE(ctx$triggered)){
                 splitted = strsplit(ctx$triggered['prop_id'][0], ".",fixed=TRUE)
                 prop_id = splitted[0]
                 prop_type = splitted[1]
               }
               processed_data = list()
               
               if(prop_id != 'ternary-layer-select'){
                 if(prop_id == "well-map" & prop_type == "selectedData"){
                   for (formation in formations) {
                     if(length(map_selected_data) == 0){
                       processed_data[formation] = NULL
                     } else{
                       processed_data[formation] = get_selection(dff, formation, map_selected_data, 0)
                     }
                   }
                 }
               else if(prop_id == "form-by-bar" & prop_type == "selectedData"){
                 processed_data = get_selection_by_bar(bar_selected_data)
                 for (formation in formations) {
                   if(length(bar_selected_data) == 0){
                     processed_data[formation] = NULL
                   } else if(!(formation %in% processed_data)){
                     processed_data[formation] = list()
                   }
                 }
                 
               } else if(rop_id == "form-by-bar" & prop_type == "clickData"){
                 processed_data = get_selection_by_bar(bar_click_data)
                 for (formation in formations) {
                   if(length(bar_click_data) == 0){
                     processed_data[formation] = NULL
                   } else if(!(formtion %in% processed_data)){
                     processed_data[formation] = list()
                   }
                 }
               } else{
                 for (formation in formations) {
                   processed_data[formation] = NULL
                 }
               }
               return(generate_ternary_map(
                 dff, processed_data, contour_visible, marker_visible))}
               
              if (prop_id == "ternary-layer-select"){
                if(length(curr_fig) == 0){
                  if(!("Well Data" %in% layer_select)){
                    marker_visible = "legendonly"
                  }
                  if (!("Rock Type" %in% layer_select)){
                    contour_visible = "legendonly"
                  }
                  for (contour_dict in curr_fig["data"][1:32,]) {
                    contour_dict["visible"] = contour_visible
                  }
                  for (marker_dict in curr_fig["data"][1:32,]) {
                    contour_dict["visible"] = marker_visible
                  }
                  return(curr_fig)
                }else{
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
               dff = df[df['op'] %in% op_select,]
               formaions = list(unique(dff['fm_name']))
               ctx = app$callback_context()
               prop_id = ""
               prop_type = ""
               
               if(isTRUE(ctx$triggered)){
                 splitted = strsplit(ctx$triggered['prop_id'][0], ".",fixed=TRUE)
                 prop_id = splitted[0]
                 prop_type = splitted[1]
               }
               
               processed_data = list()
               if (prop_id == 'ternary-map'){
                 for (formation in formations) {
                   if(length(tern_selected_data) == 0){
                     processed_data[formation] = NULL
                   }else{
                     processed_data[formation] = get_selection(dff, formation, tern_selected_data, 32)
                   }
                 }
               } else if(prop_id == 'form-by-bar'){
                 bar_data = ""
                 if(prop_type == 'selectedData'){
                   bar_data = bar_selected_data
                 } else if(prop_type == "clickData"){
                   bar_data = bar_click_data
                 }
                 processed_data = get_selection_by_bar(bar_data)
                 for (formation in formations) {
                   if(length(bar_data) == 0){
                     processed_data[formation] = NULL
                   } else if(!(formation %in% processed_data)){
                   processed_data[formation] = list()
                 }
               }
               } else{
                 for (formation in formations) {
                   processed_data[formation] = NULL
                 }}
               
               return(generate_well_map(dff, processed_data, mapbox_view))
})

app$callback(output = list(id = "production-fig", property = "figure"),
             params = list(input("well-map", "selectedData"),
                           input("ternary-map", "selectedData"),
                           input("form-by-bar", "selectedData"),
                           input("operator-select", "value")),
             function(map_select, tern_select, bar_select, op_select){
               dff = df[df$op == op_select,]
               ctx = app$callback_context()
               prop_id = ""
               prop_type = ""
               
               if(isTRUE(ctx$triggered)){
                 splitted = strsplit(ctx$triggered['prop_id'][0], ".",fixed=TRUE)
                 prop_id = splitted[0]
                 prop_type = splitted[1]
               }
               processed_data_init = list()
               processed_data_init["well_id"] = list(dff["RecordNumber"])
               processed_data_init["formation"] = list(dff["fm_name"])
             
               if(prop_id == "well-map" & prop_type == "selectedData"){
                 if(length(map_select) > 0){
                   processed_data = list('well_id' = list(), 'formation' = list())
                   for (point in map_select['points']) {
                     processed_data['well_id'] = append(processed_data['well_id'], point[['customdata']])
                     processed_data['formation'] = append(processed_data['formation'],
                                    dff[dff['RecordNumber'] == point['customdata'], list('fm_name')][1])
                   }
                 } else{
                   processed_data = processed_data_init
                 }
               } else if(prop_id == "ternary-map" & prop_type == "selectedData"){
                 if(length(tern_select) > 0){
                   processed_data = list('well_id' = list(), 'formation' = list())
                   for (point in tern_select['points']) {
                     if('customdata' %in% point){
                       processed_data["well_id"] = append(processed_data["well_id"], point["customdata"])
                       processed_data['formation'] = append(processed_data['formation'],
                                          dff[dff['RecordNumber'] == point['customdata'], list('fm_name')][1])
                     }
                   }
                 } else{
                   processed_data = processed_data_init
                 }
               } else if(prop_id == "form-by-bar" & prop_type == "selectedData"){
                 if(length(bar_select) > 0){
                   processed_data = list('well_id' = list(), 'formation' = list())
                   for (point in bar_select['points']) {
                     selected_form = point['x']
                     selected_well = list(dff[dff["fm_name"] == point["x"],"RecordNumber"])
                     for (well in selected_well) {
                       processed_data["well_id"] = append(processed_data["well_id"],integer(well))
                       processed_data["formation"] = append(processed_data["formation"],integer(selected_form))
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

app$run_server()





























