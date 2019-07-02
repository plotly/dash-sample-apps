#app = Dash$new()
leads <- getAllFields('Lead')

states = list("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DC", "DE", "FL", "GA", 
              "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", 
              "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", 
              "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", 
              "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY")

#FUNCTION FOR THE LEADS TABLE
datatotable <- function(df){
  df['CreatedDate'] = lapply(df$CreatedDate, 
  function(x)as.POSIXct(x, origin="1970-01-01 00:00:00 UTC"))
  
  df['CreatedDate'] = as.Date(format(df$CreatedDate, '%Y-%m-%d'))
  plot_ly(
    type = 'table',
    header = list(
      values = c("<b>Company</b>", names(df)),
      align = c('left', rep('center', ncol(df))),
      line = list(width = 1, color = 'white'),
      fill = list(color = 'white'),
      font = list(family = "Arial", size = 14, color = "#264e86")
    ),
    cells = list(
      values = rbind(
        rownames(df), 
        t(as.matrix(unname(df)))
      ),
      align = c('left', rep('center', ncol(df))),
      line = list(color = "white", width = 1),
      fill = list(color = c('white', 'white')),
      font = list(family = "Arial", size = 14, color = c("black"))
    ))}



choropleth_map <- function(status, df){
  if(status == 'open'){
    df = df[(df['Status'] == "Open - Not Contacted") | (df["Status"] == "Working - Contacted"),]}
  else if(status == 'converted'){
    df[(df["Status"] == "Closed - Converted"),]
  } else if(status == 'lost'){
    df = df[(df["Status"] == "Closed - Not Converted"),]
  }
  df = sqldf('SELECT COUNT(*) AS Id, State FROM df GROUP BY State')
  scl = list(list(0.0, "rgb(38, 78, 134)"), list(1.0, "#0091D5"))
  
  data = list(
    list(
      type="choropleth",
      colorscale=scl,
      locations=df$State[-c(1)],
      z=unlist(df["Id"]),
      locationmode="USA-states",
      marker=list(line=list(color="rgb(255,255,255)", width=2))
    )
  )
  layout = list(
    geo=list(
      scope="usa",
      projection=list(type="albers usa"),
      lakecolor="rgb(255, 255, 255)"
    ),
    margin=list(l=10, r=10, t=0, b=0)
  )
  
  return(list(data = data, layout = layout))
}




lead_source = function(status, df){
  if (status == 'open'){
    df = df[
      (df["Status"] == "Open - Not Contacted") | (df["Status"] == "Working - Contacted"),]
  } else if(status == 'converted'){
    df = df[df["Status"] == "Closed - Converted",]
  } else if(status == 'lost'){
    df = df[df["Status"] == "Closed - Not Converted",]
  }
  
  nb_leads = length(rownames(df))
  types = as.list(unique(df['LeadSource']))
  values = list()
  
  for (case_type in types[[1]]) {
    nb_types = dim(df[df["LeadSource"] == case_type,])[1]
    percentage = ((nb_types / nb_leads) * 100)
    values = append(values, percentage)
  }
  
  trace = plot_ly(
    labels=types[[1]],
    values=values,
    type = 'pie',
    marker=list("colors" = list("#264e86", "#0074e4", "#74dbef", "#eff0f4"))
  ) %>%
    layout(margin = list(l=15, r=10, t=0, b=65), legend=list(orientation="h"))
  
  #return(list(data=list(trace)))
}

converted_leads_count <- function(period, df){
  df = df[df["Status"] == "Closed - Converted",]
  df$CreatedDate = substr(df$CreatedDate,1,10)
  if(period == 'M'){
    df$CreatedDate = substr(df$CreatedDate,1,7)
  }
  df = count_(df %>% group_by(CreatedDate, Id, add=TRUE))
  trace = df %>% group_by(CreatedDate) %>%
    plot_ly( x = ~CreatedDate, y = ~table(CreatedDate), type = 'scatter', mode="marker",fill="tozeroy",fillcolor="#e6f2ff"
    )%>%
    layout(xaxis=list(showgrid=FALSE),
           margin=list(l=33, r=25, b=37, t=5, pad=4),
           paper_bgcolor="white",
           plot_bgcolor="white")
  
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
                  "New Lead",
                  style=list(
                    "color"= "#506784",
                    "fontWeight"= "bold",
                    "fontSize"= "20"
                  )
                ),
                htmlSpan(
                  "×",
                  id="leads_modal_close",
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
                htmlP(
                  list(
                    "Company Name"
                    
                  ),
                  style=list(
                    "float"= "left",
                    "marginTop"= "4",
                    "marginBottom"= "2"
                  ),
                  className="row"
                ),
                dccInput(
                  id="new_lead_company",
                  # placeholder="Enter company name",
                  type="text",
                  value="",
                  style=list("width"= "100%")
                ),
                htmlP(
                  "Company State",
                  style=list(
                    "textAlign"= "left",
                    "marginBottom"= "2",
                    "marginTop"= "4"
                  )
                ),
                dccDropdown(
                  id="new_lead_state",
                  options=lapply(states,function(state){list(label = state, value = state)}),
                  value="NY"
                ),
                htmlP(
                  "Status",
                  style=list(
                    "textAlign"= "left",
                    "marginBottom"= "2",
                    "marginTop"= "4"
                  )
                ),
                dccDropdown(
                  id="new_lead_status",
                  options=list(
                    list(
                      "label"= "Open - Not Contacted",
                      "value"= "Open - Not Contacted"
                    ),
                    list(
                      "label"= "Working - Contacted",
                      "value"= "Working - Contacted"
                    ),
                    list(
                      "label"= "Closed - Converted",
                      "value"= "Closed - Converted"
                    ),
                    list(
                      "label"= "Closed - Not Converted",
                      "value"= "Closed - Not Converted"
                    )
                  ),
                  value="Open - Not Contacted"
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
                  id="new_lead_source",
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
                )
              ),
              className="row",
              style=list("padding"= "2% 8%")
            ),
            
            # submit button
            htmlSpan(
              "Submit",
              id="submit_new_lead",
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
    id="leads_modal",
    style=list("display"= "none")
  ))}

leadslayout <- app$layout(
  
  htmlDiv(toJSON(getAllFields('Lead')), id="leads_df", style=list("display"= "none")),
  
  # top controls
  htmlDiv(
    list(
      htmlDiv(
        dccDropdown(
          id="converted_leads_dropdown",
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
          id="lead_source_dropdown",
          options=list(
            list("label"= "All status", "value"= "all"),
            list("label"= "Open leads", "value"= "open"),
            list("label"= "Converted leads", "value"= "converted"),
            list("label"= "Lost leads", "value"= "lost")
          ),
          value="all",
          clearable=FALSE
        ),
        className="two columns"
      ),
      
      # add button
      htmlDiv(
        htmlSpan(
          "Add new",
          id="new_lead",
          n_clicks=0,
          className="button button--primary",
          style=list(
            "height"= "34",
            "background"= "#119DFF",
            "border"= "1px solid #119DFF",
            "color"= "white"
          )
        ),
        className="two columns",
        style=list("float"= "right")
      )
    ),
    className="row",
    style=list("marginBottom"= "10")
  ),
  
  # indicators row div
  htmlDiv(
    list(
      indicator(
        "#00cc96", "Converted Leads", "left_leads_indicator"
      ),
      indicator(
        "#119DFF", "Open Leads", "middle_leads_indicator"
      ),
      indicator(
        "#EF553B",
        "Conversion Rates",
        "right_leads_indicator"
      )
    ),
    className="row"
  ),
  
  # charts row div
  htmlDiv(
    list(
      htmlDiv(
        list(
          htmlP("Leads count per state" ),
          dccGraph(
            id="map",
            figure = choropleth_map('new',leads),
            style=list("height"= "90%", "width"= "98%"),
            config=list(displayModeBar=FALSE)
          )
        ),
        className="four columns chart_div"
      ),
      
      htmlDiv(
        list(
          htmlP("Leads by source"),
          dccGraph(
            id="lead_source",
            figure = lead_source('new',leads),
            style=list("height"= "90%", "width"= "98%"),
            config=list(displayModeBar=FALSE)
          )
        ),
        className="four columns chart_div"
      ),
      
      htmlDiv(
        list(
          htmlP("Converted Leads count"),
          dccGraph(
            id="converted_leads",
            style=list("height"= "90%", "width"= "98%"),
            config=list(displayModeBar=FALSE)
          )
        ),
        className="four columns chart_div"
      )
    ),
    className="row",
    style=list("marginTop"= "5")
  ),
  
  # table div
  dccGraph(id="leads_table",
           className="row",
           style=list(
             "maxHeight"= "350px",
             "overflowY"= "scroll",
             "padding"= "8",
             "marginTop"= "5",
             "backgroundColor"="white",
             "border"= "1px solid #C8D4E3",
             "borderRadius"= "3px"
           )),
  modal()
)


app$callback(output = list(id = "left_leads_indicator", property = "children"),
             params = list(
               input(id = "leads_df", property = "children")),
             function(df){
               df = data.frame(fromJSON(df))
               converted_leads = nrow(unique(df[df["Status"] == "Closed - Converted",]))
               return(converted_leads)
             })

app$callback(output = list(id = "middle_leads_indicator", property = "children"),
             params = list(
               input(id = "leads_df", property = "children")),
             function(df){
               df = data.frame(fromJSON(df))
               open_leads = nrow(df[(df["Status"] == "Open - Not Contacted")| (df["Status"] == "Working - Contacted"),])
               return(open_leads)
             })

app$callback(output = list(id = "right_leads_indicator", property = "children"),
             params = list(
               input(id = "leads_df", property = "children")),
             function(df){
               df = data.frame(fromJSON(df))
               converted_leads = nrow(df[df["Status"] == "Closed - Converted",])
               lost_leads = nrow(df[df["Status"] == "Closed - Not Converted",])
               conversion_Rates = converted_leads/(converted_leads + lost_leads) * 100
               return(sprintf("%.2f %%", conversion_Rates))
             })



app$callback(output = list(id = "lead_source", property = "figure"),
             params = list(
               input(id = "lead_source_dropdown", property = "value"),
               input(id = "leads_df", property = "children")),
             function(status,df){
               df = data.frame(fromJSON(df))
               return(lead_source(status,df))
             })


app$callback(output = list(id = "map", property = "figure"),
             params = list(
               input(id = "lead_source_dropdown", property = "value"),
               input(id = "leads_df", property = "children")),
             function(status,df){
               df = data.frame(fromJSON(df))
               return(choropleth_map(status,df))
             })

app$callback(output = list(id = "leads_table", property = "figure"),
             params = list(
               input(id = "lead_source_dropdown", property = "value"),
               input(id = "leads_df", property = "children")),
             function(status,df){
               df = data.frame(fromJSON(df))
               if(status == 'open'){
                 df = df[
                   (df["Status"] == "Open - Not Contacted") | (df["Status"] == "Working - Contacted"),]
               } else if (status == 'converted'){
                 df = df[df["Status"] == "Closed - Converted",]
               } else if (status == 'lost'){
                 df = df[df["Status"] == "Closed - Not Converted",]
               }
               df = df[,c("CreatedDate", "Status", "Company", "State", "LeadSource")]
               return(datatotable(df))
             })


app$callback(output = list(id = "converted_leads", property = "figure"),
             params = list(
               input(id = "converted_leads_dropdown", property = "value"),
               input(id = "leads_df", property = "children")),
             function(period,df){
               df = data.frame(fromJSON(df))
               return(converted_leads_count(period, df)) #Fixxxx
             })



app$callback(output = list(id = "leads_modal", property = "style"),
             params = list(
               input(id = "new_lead", property = "n_clicks")),
             function(n){
               if(n > 0){
                 return(list('display' = 'block'))
               }
               return(list('display' = 'none'))
             })

app$callback(output = list(id = "new_lead", property = "n_clicks"),
             params = list(
               input(id = "leads_modal_close", property = "n_clicks"),
               input(id = "submit_new_lead", property = "n_clicks")),
             function(n,n2){
               return(0)
             })

app$callback(output = list(id = "leads_df", property = "children"),
             params = list(
               input(id = "submit_new_lead", property = "n_clicks"),
               state("new_lead_status", "value"),
               state("new_lead_state", "value"),
               state("new_lead_company", "value"),
               state("new_lead_source", "value"),
               state("leads_df", "children")),
             function(n_clicks, status, state, company, source, current_df){
               if(n_clicks > 0){
                 if(company == ""){
                   company = 'Not named yet'
                 }
                 add_leads(company,status,state,source)
                 df = getAllFields('Lead')
                 return(toJSON(df))
               }
               return(current_df)
             })

#app$run_server()   