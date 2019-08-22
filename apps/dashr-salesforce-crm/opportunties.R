#app = Dash$new()

opp = get_opportunities()

converted_opportunities = function(period, source, df){
  
  if(source == 'all_s'){
    df = df[df["IsWon"] == 1,]
  } else{
    df = df[(df["LeadSource"] == source) & (df["IsWon"] == 1),]
  }
  if (period == "W-MON"){
    df['CreatedDate'] = as.Date(df$CreatedDate) - 7
  } else if (period == "M"){
    df$CreatedDate = format(as.Date(df$CreatedDate), "%B %Y")
  }
  df = count_(df %>% group_by(CreatedDate, add=TRUE))
  
  if(length(df == 0)){
    layout = list(annotation = list(
      list(text = 'No results found',
           showarrow = FALSE)
    ))
    #return(list('data' = list(), 'layout'=layout))
  }
  
  #df['CreatedDate'] = lapply(df$CreatedDate, function(x)as.POSIXct(x, origin="1970-01-01 00:00:00 UTC"))
  trace = plot_ly(
    x=df$CreatedDate,
    y=df$n,
    type = 'scatter',
    name="converted opportunities",
    fill="tozeroy",
    fillcolor="#e6f2ff",
    mode = 'markers'
  )%>%
    layout(
      xaxis=list(showgrid=FALSE),
      margin=list(l=35, r=25, b=23, t=5, pad=4),
      paper_bgcolor="white",
      plot_bgcolor="white")
  
  return(list(data = list(trace)))
}

#heat_map_fig(opp,opp$StageName,opp$Type)

heat_map_fig= function(df,x,y){
  z = list()
  for (i in y) {
    z_row = list()
    for (j in x) {
      probability = mean(df[(df["StageName"] == j & df["Type"] == i),
                            'Probability'])
      z_row = append(z_row, probability)
    }
    z = append(z,z_row)
  }
  trace = list(
    type="heatmap", 
    z=z, x=x, y=y, 
    name="mean probability", 
    colorscale="Blues"
  )
  
  layout = list(
    margin=list(t=25, l=210, b=85, pad=4),
    paper_bgcolor="white",
    plot_bgcolor="white"
  )
  
  return(list(data= list(trace), layout = layout))
}

top_open_opportunities = function(df){
  df = data.frame(df)
  df = df[order(df['Amount']),] 
  cols = c("CreatedDate", "Name", "Amount", "StageName")
  df = df[c(1:5),cols]
  df['Name'] = lapply(df['Name'], function(x){x[1:30]})
  return(datatotable(df))
}

top_lost_opportunities = function(df){
  df = data.frame(df)
  df = df[df['StageName'] == 'Closed Lost',]
  cols = c("CreatedDate", "Name", "Amount", "StageName")
  df = df[c(1:5),cols]
  df['Name'] = lapply(df['Name'], function(x){x[1:30]})
  return(datatotable(df))
}

modal = function(){
  return(htmlDiv(
    htmlDiv(
      list(
        htmlDiv(
          list(
            
            # modal header
            htmlDiv(
              list(
                htmlSpan(
                  "New Opportunity",
                  style=list(
                    "color"= "#506784",
                    "fontWeight"= "bold",
                    "fontSize"= "20"
                  )
                ),
                htmlSpan(
                  "×",
                  id="opportunities_modal_close",
                  n_clicks=0,
                  style=list(
                    "float"= "right",
                    "cursor"= "pointer",
                    "marginTop"= "0",
                    "marginBottom"= "17"
                  )
                )
              ),
              className="row",
              style=list("borderBottom"= "1px solid #C8D4E3")
            ),
            
            
            # modal form 
            htmlDiv(
              list(
                
                # left div
                htmlDiv(
                  list(
                    htmlP(
                      list(
                        "Name"
                      ),
                      style=list(
                        "float"= "left",
                        "marginTop"= "4",
                        "marginBottom"= "2"
                      ),
                      className="row"
                    ),
                    dccInput(
                      id="new_opportunity_name",
                      placeholder="Name of the opportunity",
                      type="text",
                      value="",
                      style=list("width"= "100%")
                    ),
                    
                    htmlP(
                      list(
                        "StageName"
                      ),
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccDropdown(
                      id="new_opportunity_stage",
                      options=list(
                        list(
                          "label"= "Prospecting",
                          "value"= "Prospecting"
                        ),
                        list(
                          "label"= "Qualification",
                          "value"= "Qualification"
                        ),
                        list(
                          "label"= "Needs Analysis",
                          "value"= "Needs Analysis"
                        ),
                        list(
                          "label"= "Value Proposition",
                          "value"= "Value Proposition"
                        ),
                        list(
                          "label"= "Id. Decision Makers",
                          "value"= "Closed"
                        ),
                        list(
                          "label"= "Perception Analysis",
                          "value"= "Perception Analysis"
                        ),
                        list(
                          "label"= "Proposal/Price Quote",
                          "value"= "Proposal/Price Quote"
                        ),
                        list(
                          "label"= "Negotiation/Review",
                          "value"= "Negotiation/Review"
                        ),
                        list(
                          "label"= "Closed/Won",
                          "value"= "Closed Won"
                        ),
                        list(
                          "label"= "Closed/Lost",
                          "value"= "Closed Lost"
                        )
                      ),
                      clearable=FALSE,
                      value="Prospecting"
                    ),
                    
                    htmlP(
                      "Source",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccDropdown(
                      id="new_opportunity_source",
                      options=list(
                        list("label"= "Web", "value"= "Web"),
                        list(
                          "label"= "Phone Inquiry",
                          "value"= "Phone Inquiry"
                        ),
                        list(
                          "label"= "Partner Referral",
                          "value"= "Partner Referral"
                        ),
                        list(
                          "label"= "Purchased List",
                          "value"= "Purchased List"
                        ),
                        list("label"= "Other", "value"= "Other")
                      ),
                      value="Web"
                    ),
                    
                    htmlP(
                      list(
                        "Close Date"
                      ),
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    htmlDiv(
                      dccDatePickerSingle(
                        id="new_opportunity_date",
                        min_date_allowed=as.Date(Sys.Date()),
                        # max_date_allowed=dt(2017, 9, 19),
                        initial_visible_month=as.Date(Sys.Date()),
                        date=as.Date(Sys.Date())
                      ),
                      style=list("textAlign"= "left")
                    )
                    
                  ),
                  className="six columns",
                  style=list("paddingRight"= "15")
                ),
                
                
                # right div
                htmlDiv(
                  list(
                    htmlP(
                      "Type",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccDropdown(
                      id="new_opportunity_type",
                      options=list(
                        list(
                          "label"= "Existing Customer - Replacement",
                          "value"= "Existing Customer - Replacement"
                        ),
                        list(
                          "label"= "New Customer",
                          "value"= "New Customer"
                        ),
                        list(
                          "label"= "Existing Customer - Upgrade",
                          "value"= "Existing Customer - Upgrade"
                        ),
                        list(
                          "label"= "Existing Customer - Downgrade",
                          "value"= "Existing Customer - Downgrade"
                        )
                      ),
                      value="New Customer"
                    ),
                    
                    htmlP(
                      "Amount",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccInput(
                      id="new_opportunity_amount",
                      placeholder="0",
                      type="number",
                      value="",
                      style=list("width"= "100%")
                    ),
                    
                    htmlP(
                      "Probability",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccInput(
                      id="new_opportunity_probability",
                      placeholder="0",
                      type="number",
                      max=100,
                      step=1,
                      value="",
                      style=list("width"= "100%")
                    )
                    
                  ),
                  className="six columns",
                  style=list("paddingLeft"= "15")
                )
              ),
              className="row",
              style=list("paddingTop"= "2%")
            ),
            
            
            # submit button
            htmlSpan(
              "Submit",
              id="submit_new_opportunity",
              n_clicks=0,
              className="button button--primary add"
            )
          ),
          className="modal-content",
          style=list("textAlign"= "center")
        )
      ),
      className="modal"
    ),
    id="opportunities_modal",
    style=list("display"= "none")
  ))}

Opplayout = app$layout(
  htmlDiv(
    toJSON(get_opportunities()),  # opportunities df
    id="opportunities_df",
    style=list("display"= "none")
  ),
  # top controls
  htmlDiv(
    list(
      htmlDiv(
        dccDropdown(
          id="converted_opportunities_dropdown",
          options=list(
            list("label"= "By day", "value"= "D"),
            list("label"= "By week", "value"= "W-MON"),
            list("label"= "By month", "value"= "M")
          ),
          value="D",
          clearable=FALSE
        ),
        className="two columns"
      ),
      
      htmlDiv(
        dccDropdown(
          id="heatmap_dropdown",
          options=list(
            list("label"= "All stages", "value"= "all_s"),
            list("label"= "Cold stages", "value"= "cold"),
            list("label"= "Warm stages", "value"= "warm"),
            list("label"= "Hot stages", "value"= "hot")
          ),
          value="all_s",
          clearable=FALSE
        ),
        className="two columns"
      ),
      htmlDiv(
        dccDropdown(
          id="source_dropdown",
          options=list(
            list("label"= "All sources", "value"= "all_s"),
            list("label"= "Web", "value"= "Web"),
            list("label"= "Word of Mouth", "value"= "Word of mouth"),
            list("label"= "Phone Inquiry", "value"= "Phone Inquiry"),
            list("label"= "Partner Referral", "value"= "Partner Referral"),
            list("label"= "Purchased List", "value"= "Purchased List"),
            list("label"= "Other", "value"= "Other")
          ),
          value="all_s",
          clearable=FALSE
        ),
        className="two columns"
      ),
      
      #add button
      htmlDiv(
        htmlSpan(
          "Add new",
          id="new_opportunity",
          n_clicks=0,
          className="button button--primary add"
        ),
        className="two columns",
        style=list("float"= "right")
      )
    ),
    className="row",
    style=list("marginBottom"= "10")
  ),
  
  #indicators row
  htmlDiv(
    list(
      indicator(
        "#00cc96",
        "Won opportunities",
        "left_opportunities_indicator"
      ),
      indicator(
        "#119DFF",
        "Open opportunities",
        "middle_opportunities_indicator"
      ),
      indicator(
        "#EF553B",
        "Lost opportunities",
        "right_opportunities_indicator"
      )
    ),
    className="row"
  ),
  
  #charts row div 
  htmlDiv(
    list(
      htmlDiv(
        list(
          htmlP("Converted Opportunities count"),
          dccGraph(
            id="converted_count",
            #figure = converted_opportunities('W-MON','all_s',opp)[['data']][[1]],
            style=list("height"= "90%", "width"= "98%"),
            config=list(displayModeBar=FALSE)
          )
        ),
        className="four columns chart_div"
      ),
      
      htmlDiv(
        list(
          htmlP("Probabilty heatmap per Stage and Type"),
          dccGraph(
            id="heatmap",
            #figure = heat_map_fig(opp,opp$StageName,opp$Type),
            style=list("height"= "90%", "width"= "98%"),
            config=list(displayModeBar=FALSE)
          )
        ),
        className="eight columns chart_div"
      )
    ),
    className="row",
    style=list("marginTop"= "5px")
  ),
  
  
  # tables row div
  htmlDiv(
    list(
      htmlDiv(
        list(
          htmlP(
            "Top Open opportunities",
            style=list(
              "color"= "#2a3f5f",
              "fontSize"= "13px",
              "textAlign"= "center",
              "marginBottom"= "0"
            )
          ),
          htmlDiv(
            id="top_open_opportunities",
            style=list("padding"= "0px 13px 5px 13px", "marginBottom"= "5")
          )
          
        ),
        className="six columns",
        style=list(
          "backgroundColor"= "white",
          "border"= "1px solid #C8D4E3",
          "borderRadius"= "3px",
          "height"= "100%",
          "overflowY"= "scroll"
        )
      ),
      htmlDiv(
        list(
          htmlP(
            "Top Lost opportunities",
            style=list(
              "color"= "#2a3f5f",
              "fontSize"= "13px",
              "textAlign"= "center",
              "marginBottom"= "0"
            )
          ),
          htmlDiv(
            id="top_lost_opportunities",
            style=list("padding"= "0px 13px 5px 13px", "marginBottom"= "5")
          )
        ),
        className="six columns",
        style=list(
          "backgroundColor"= "white",
          "border"= "1px solid #C8D4E3",
          "borderRadius"= "3px",
          "height"= "100%",
          "overflowY"= "scroll"
        )
      ),
      
      
      
      modal()
    ),
    className="row",
    style=list("marginTop"= "5px", "max height"= "200px")
  )
)

app$callback(output = list(id ="left_opportunities_indicator", property = "children"),
             params = list(
               input(id = "opportunities_df", property = "children")),
             function(df){
               df = as.data.frame(fromJSON(df))
               won = sum(df[df['IsWon'] == 1, c('Amount')])
               sprintf('%.1f M', won/1000000)
             })

app$callback(output = list(id ="middle_opportunities_indicator", property = "children"),
             params = list(
               input(id = "opportunities_df", property = "children")),
             function(df){
               df = as.data.frame(fromJSON(df))
               active = sum(df[(df['IsClosed'] == 0), c('Amount')])
               sprintf('%.1f M', floor(active/1000000))
             })

app$callback(output = list(id ="right_opportunities_indicator", property = "children"),
             params = list(
               input(id = "opportunities_df", property = "children")),
             function(df){
               df = as.data.frame(fromJSON(df))
               lost = sum(df[(df["IsWon"] == 0) & (df["IsClosed"] == 1), c('Amount')])
               return(lost)
             })

app$callback(output = list(id = 'converted_count', property = 'figure'),
             params = list(
               input(id = "converted_opportunities_dropdown", property = "value"),
               input(id = "source_dropdown", property = "value"),
               input(id = "opportunities_df", property = "children")),
             function(period, source, df){
               df = data.frame(fromJSON(df))
               return(converted_opportunities(period, source,df)[['data']][[1]])
             })


app$callback(output = list(id = 'heatmap', property = 'figure'),
             params = list(
               input(id = 'heatmap_dropdown', property = 'value'),
               input(id = 'opportunities_df', property = 'children')),
             function(stage,df){
               df = data.frame(fromJSON(df))
               df = df[complete.cases(df['Type']),]
               x = list()
               y = unique(df['Type'])[['Type']]
               if (stage == 'all_s'){
                 x = unique(df['StageName'])[['StageName']]
               } else if (stage == 'cold'){
                 x = list("Needs Analysis", "Prospecting", "Qualification")
               } else if (stage == 'warm'){
                 x = list("Value Proposition", "Id. Decision Makers", "Perception Analysis")
               } else {
                 x = list("Proposal/Price Quote", "Negotiation/Review", "Closed Won")
               }
               return(heat_map_fig(df,x,y))
             })


app$callback(output = list(id ="opportunities_modal", property = "style"),
             params = list(
               input(id = "new_opportunity", property = "n_clicks")),
             function(n){
               if(n > 0){
                 return(list('display' = 'block'))
               }
               return(list('display' = 'none'))
             })

app$callback(output = list(id ="new_opportunity", property = "n_clicks"),
             params = list(
               input(id = "opportunities_modal_close", property = "n_clicks"),
               input(id = "submit_new_opportunity", property = "n_clicks")),
             function(n,n2){
               return(0)
             })

app$callback(output = list(id ="opportunities_df", property = "children"),
             params = list(
               input(id = "submit_new_opportunity", property = "n_clicks"),
               state("new_opportunity_name", "value"),
               state("new_opportunity_stage", "value"),
               state("new_opportunity_amount", "value"),
               state("new_opportunity_probability", "value"),
               state("new_opportunity_date", "date"),
               state("new_opportunity_type", "value"),
               state("new_opportunity_source", "value"),
               state("opportunities_df", "children")),
             function(n_clicks, name, stage, amount, probability, date, o_type, source, current_df){
               if(n_clicks > 0){
                 if (name == ""){
                   name = 'Not named yet'
                 }
                 add_opportunities(name, stage, amount, probability, date, o_type, source)
                 df = get_opportunities()
                 return(toJSON(df))
               }
               return(current_df)
             })

app$callback(output = list(id ="top_open_opportunities", property = "children"),
             params = list(
               input(id = "opportunities_df", property = "children")),
             function(df){
               df = data.frame(fromJSON(df))
               return(top_open_opportunities(df))
             })

app$callback(output = list(id ="top_lost_opportunities", property = "children"),
             params = list(
               input(id = "opportunities_df", property = "children")),
             function(df){
               df = data.frame(fromJSON(df))
               return(top_lost_opportunities(df))
             })

#app$run_server()