
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(jsonlite)

#app = Dash$new()
colors = list('background' = "#F3F6FA",
              "background_div" = "white")

pie_chart = function(df, column, priority, origin){
  df = df[complete.cases(df$Type,df$Reason,df$Origin),]
  nb_cases = length(rownames(df))
  types = list()
  values = list()
  
  if(priority == 'all_p'){
    if (origin == 'all'){
      types = unique(df[column])
      types = lapply(as.matrix(types), function(x)x)
    } else{
      types = unique(df[df['Origin'] == origin, column])
    }
  } 
  else{
    
    if(origin == 'all'){
      types = unique(df[df['Priority'] == priority, column])
    } else{
      
      types = unique(df[(df['Priority'] == priority) & (df['Origin'] == origin),
                        column])
    }}
  
  if(is.null(types)){
    layout = list(annotations = list('text' = 'No results found',
                                     'showarrow' = FALSE))
    return(list('data' = list(),
                'layout' = layout))
  }
  
  for (case_type in types) {
    nb_type = dim(df[df[column] == case_type,])[1]
    values = append(values,(nb_type/nb_cases)*100)
  }
  
  traces = plot_ly(labels = types, values = values, type = 'pie', orientation = 'h',
                   marker = list('colors' =list("#264e86", "#0074e4", "#74dbef", "#eff0f4")))%>%
    layout(
      margin = list(l = 0, r = 0, b = 0, t=4, pad=8),
      legend = list(orientation = 'h'),
      paper_bgcolor = 'white',
      plot_bgcolor = 'white'
    )
  
  return(list(data = list(traces)))
  
}

cases_by_period = function(df, period, priority, origin){
  df = df[complete.cases(df$Type,df$Reason,df$Origin),]
  df$CreatedDate = substr(df$CreatedDate,1,10)
  stages = unique(df['Type'])
  stages = stages$Type[c(1:5)]
  
  if(priority != 'all_p'){
    df = df[df['Priority'] == priority,]
  }
  
  #Period filtering
  #df['CreatedDate'] = lapply(df$CreatedDate, function(x)as.POSIXct(x, origin="1970-01-01 00:00:00 UTC"))
  #df['CreatedDate'] = format(df$CreatedDate, '%Y-%m-%d')
  if (period == 'W-MON'){
    df$CreatedDate = as.Date(df$CreatedDate)- 7
  } else if (period == 'M'){
    df$CreatedDate = format(as.Date(df$CreatedDate), "%B %Y")
  }
  
  df = count_(df %>% group_by(CreatedDate, Type, add=TRUE))
  dates = list(unique(df$CreatedDate))
  
  co = list("Electrical"= "#264e86",
            "Other" = "#0074e4",
            "Structural" = "#74dbef",
            "Mechanical" = "#eff0f4",
            "Electronic" = "#ff7f0e"
  )
  colnames(df)[3] = 'IsDeleted'
  data = list()
  stage_rows = list()
  data_trace = df %>% group_by(CreatedDate) %>% arrange(Type) %>%
    plot_ly( x = ~CreatedDate, y = ~IsDeleted, name= ~Type ,marker = list(color = c(co)), type = 'bar'
    )%>%
    layout(barmode="stack",
           margin=list(l=40, r=25, b=40, t=0, pad=4),
           paper_bgcolor="white",
           plot_bgcolor="white")
  
  return(list('data' = list(data_trace)))
}

cases_by_account = function(cases){
  cases = cases[complete.cases(cases$AccountId),]
  colnames(accounts)[1] = 'AccountId'
  cases <- merge(cases, accounts, by="AccountId", add = TRUE)
  cases = sqldf('SELECT AccountId, Name, COUNT(*) FROM cases GROUP BY AccountID, Name')  #****
  colnames(cases)[3] = 'IsDeleted'
  cases = cases[order(cases$IsDeleted),]
  yform <- list(categoryorder = "array",
                categoryarray = c(cases$Name))
  data = plot_ly(
    x = cases$IsDeleted,
    y = cases$Name,
    type = "bar") %>%
    layout(
      yaxis = yform,
      barmode="stack",
      margin=list(l=210, r=25, b=20, t=0, pad=4),
      paper_bgcolor="white",
      plot_bgcolor="white"
    )
  return(list('data' = data))
}

modal = function(){
  contacts["Name"] = (paste(contacts$Salutation, contacts$FirstName, contacts$LastName))
  contactoptions = list()
  k = 1
  for (i in 1:nrow(contacts)) {
    contactoptions[[k]] = list('label' = contacts$Name[i], 'value' = contacts$Id[i])
    k = k+1
  }
  return(htmlDiv(
    htmlDiv(
      list(
        htmlDiv(
          list(
            # modal header
            htmlDiv(
              list(
                htmlSpan(
                  "New Case",
                  style=list(
                    "color"= "#506784",
                    "fontWeight"= "bold",
                    "fontSize"= "20"
                  )
                ),
                htmlSpan(
                  "×",
                  id="cases_modal_close",
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
                
                # left Div
                htmlDiv(
                  list(
                    htmlP(
                      "Account name",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    htmlDiv(
                      dccDropdown(
                        id="new_case_account",
                        options=lapply(accounts$Name, function(x){ list('label' = x, 'value' = x)}),
                        clearable=FALSE
                        #value=accounts.iloclist(0).Id,
                      )
                    ),
                    htmlP(
                      "Priority",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccDropdown(
                      id="new_case_priority",
                      options=list(
                        list("label"= "High", "value"= "High"),
                        list("label"= "Medium", "value"= "Medium"),
                        list("label"= "Low", "value"= "Low")
                      ),
                      value="Medium",
                      clearable=FALSE
                    ),
                    htmlP(
                      "Origin",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccDropdown(
                      id="new_case_origin",
                      options=list(
                        list("label"= "Phone", "value"= "Phone"),
                        list("label"= "Web", "value"= "Web"),
                        list("label"= "Email", "value"= "Email")
                      ),
                      value="Phone",
                      clearable=FALSE
                    ),
                    htmlP(
                      "Reason",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccDropdown(
                      id="new_case_reason",
                      options=list(
                        list(
                          "label"= "Installation",
                          "value"= "Installation"
                        ),
                        list(
                          "label"= "Equipment Complexity",
                          "value"= "Equipment Complexity"
                        ),
                        list(
                          "label"= "Performance",
                          "value"= "Performance"
                        ),
                        list(
                          "label"= "Breakdown",
                          "value"= "Breakdown"
                        ),
                        list(
                          "label"= "Equipment Design",
                          "value"= "Equipment Design"
                        ),
                        list(
                          "label"= "Feedback",
                          "value"= "Feedback"
                        ),
                        list("label"= "Other", "value"= "Other")
                      ),
                      value="Installation",
                      clearable=FALSE
                    ),
                    htmlP(
                      "Subject",
                      style=list(
                        "float"= "left",
                        "marginTop"= "4",
                        "marginBottom"= "2"
                      ),
                      className="row"
                    ),
                    dccInput(
                      id="new_case_subject",
                      placeholder="The Subject of the case",
                      type="text",
                      value="",
                      style=list("width"= "100%")
                    )
                  ),
                  className="six columns",
                  style=list("paddingRight"= "15")
                ),
                
                
                # right Div
                htmlDiv(
                  list(
                    htmlP(
                      "Contact name",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    htmlDiv(
                      dccDropdown(
                        id="new_case_contact",
                        options=contactoptions,
                        clearable=FALSE
                        #value=contacts.iloclist(0).Id,
                      )
                    ),
                    htmlP(
                      "Type",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccDropdown(
                      id="new_case_type",
                      options=list(
                        list(
                          "label"= "Electrical",
                          "value"= "Electrical"
                        ),
                        list(
                          "label"= "Mechanical",
                          "value"= "Mechanical"
                        ),
                        list(
                          "label"= "Electronic",
                          "value"= "Electronic"
                        ),
                        list(
                          "label"= "Structural",
                          "value"= "Structural"
                        ),
                        list("label"= "Other", "value"= "Other")
                      ),
                      value="Electrical"
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
                      id="new_case_status",
                      options=list(
                        list("label"= "New", "value"= "New"),
                        list(
                          "label"= "Working",
                          "value"= "Working"
                        ),
                        list(
                          "label"= "Escalated",
                          "value"= "Escalated"
                        ),
                        list("label"= "Closed", "value"= "Closed")
                      ),
                      value="New"
                    ),
                    htmlP(
                      "Supplied Email",
                      style=list(
                        "textAlign"= "left",
                        "marginBottom"= "2",
                        "marginTop"= "4"
                      )
                    ),
                    dccInput(
                      id="new_case_email",
                      placeholder="email",
                      type="email",
                      value="",
                      style=list("width"= "100%")
                    ),
                    htmlP(
                      "Description",
                      style=list(
                        "float"= "left",
                        "marginTop"= "4",
                        "marginBottom"= "2"
                      ),
                      className="row"
                    ),
                    dccTextarea(
                      id="new_case_description",
                      placeholder="Description of the case",
                      value="",
                      style=list("width"= "100%")
                    )
                  ),
                  className="six columns",
                  style=list("paddingLeft"= "15")
                )
              ),
              style=list("marginTop"= "10", "textAlign"= "center"),
              className="row"
            ),
            
            # submit button
            htmlSpan(
              "Submit",
              id="submit_new_case",
              n_clicks=0,
              className="button button--primary add"
            )
          ),
          className="modal-content",
          style=list("textAlign"= "center", "border"= "1px solid #C8D4E3")
        )
      ),
      className="modal"
    ),
    id="cases_modal",
    style=list("display"= "none")
  ))}

caseslayout = app$layout(
  
  htmlDiv(toJSON(get_cases()), id="cases_df", style=list("display"= "none")),
  
  modal(),
  # top controls
  htmlDiv(
    list(
      htmlDiv(
        dccDropdown(
          id="cases_period_dropdown",
          options=list(
            list("label"= "By day", "value"= "D"),
            list("label"= "By week", "value"= "W-MON"),
            list("label"= "By month", "value"= "M")
          ),
          value="D",
          clearable=FALSE
        ),
        className="two columns",
        style=list("marginBottom"= "10")
      ),
      htmlDiv(
        dccDropdown(
          id="priority_dropdown",
          options=list(
            list("label"= "All priority", "value"= "all_p"),
            list("label"= "High priority", "value"= "High"),
            list("label"= "Medium priority", "value"= "Medium"),
            list("label"= "Low priority", "value"= "Low")
          ),
          value="all_p",
          clearable=FALSE
        ),
        className="two columns"
      ),
      htmlDiv(
        dccDropdown(
          id="origin_dropdown",
          options=list(
            list("label"= "All origins", "value"= "all"),
            list("label"= "Phone", "value"= "Phone"),
            list("label"= "Web", "value"= "Web"),
            list("label"= "Email", "value"= "Email")
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
          id="new_case",
          n_clicks=0,
          className="button button--primary add"
          
        ),
        className="two columns",
        style=list("float"= "right")
      )
    ),
    className="row",
    style=list()
  ),
  
  # indicators div 
  htmlDiv(
    list(
      indicator(
        "#00cc96",
        "Low priority cases",
        "left_cases_indicator"
      ),
      indicator(
        "#119DFF",
        "Medium priority cases",
        "middle_cases_indicator"
      ),
      indicator(
        "#EF553B",
        "High priority cases",
        "right_cases_indicator"
      )
    ),
    className="row"
  ),
  
  
  htmlDiv(
    list(
      htmlDiv(
        list(
          htmlP("Cases Type"),
          
          dccGraph(
            id="cases_types",
            config=list(displayModeBar=FALSE),
            #figure=pie_chart(cases,'Type','all_p', 'all')[[1]][[1]],
            style=list("height"= "89%", "width"= "98%")
          )
          
        ),
        className="six columns chart_div"
      ),
      
      htmlDiv(
        list(
          htmlP("Cases Reasons"),
          
          dccGraph(
            id="cases_reasons",
            #figure=pie_chart(cases,'Reason','all_p', 'all')[[1]][[1]],
            config=list(displayModeBar=FALSE),
            style=list("height"= "89%", "width"= "98%")
          )
        ),
        className="six columns chart_div"
      )
    ),
    className="row",
    style=list("marginTop"= "5px")
  ),
  
  
  htmlDiv(
    list(
      htmlDiv(
        list(
          htmlP("Cases over Time"),
          dccGraph(
            id="cases_by_period",
            #figure = cases_by_period(cases,'W-MON','High','all')[[1]][[1]],
            config=list(displayModeBar=FALSE),
            style=list("height"= "89%", "width"= "98%")
          )
        ),
        className="six columns chart_div"
      ),
      
      htmlDiv(
        list(
          htmlP("Cases by Company"),
          dccGraph(
            id="cases_by_account",
            #cases_by_account(cases)[[1]],
            config=list(displayModeBar=FALSE),
            style=list("height"= "87%", "width"= "98%")
          )
        ),
        className="six columns chart_div"
      )
    ),
    className="row",
    style=list("marginTop"= "5px")
  )
  
)

app$callback(output=list(id="left_cases_indicator", property='children'),
             params=list(
               input(id='cases_df', property='children')),
             function(df){
               df = data.frame(fromJSON(df))
               low = length(df[(df['Priority'] == 'Low') & (df['Status'] == 'New'), 'Priority'])
               return(low)
             })

app$callback(output=list(id="middle_cases_indicator", property='children'),
             params=list(
               input(id="cases_df", property='children')),
             function(df){
               df = data.frame(fromJSON(df))
               medium = length(df[(df['Priority'] == 'Medium') & (df['Status'] == 'New'), 'Priority'])
               return(medium)
             })
app$callback(output=list(id="right_cases_indicator", property='children'),
             params=list(
               input(id="cases_df", property='children')),
             function(df){
               df = data.frame(fromJSON(df))
               high = length(df[(df['Priority'] == 'High') & (df['Status'] == 'New'), 'Priority'])
               return(high)
             })

app$callback(output = list(id = "cases_reasons", property = 'figure'),
             params = list(
               input(id = "priority_dropdown", property = 'value'),
               input(id = "origin_dropdown", property = 'value'),
               input(id = "cases_df", property = 'children')),
             function(priority,origin, df){
               df = data.frame(fromJSON(df))
               return(pie_chart(df,'Reason', priority, origin)[[1]][[1]])
             })

app$callback(output = list(id = "cases_types", property = 'figure'),
             params = list(
               input(id = "priority_dropdown", property = 'value'),
               input(id = "origin_dropdown", property = 'value'),
               input(id = "cases_df", property = 'children')),
             function(priority,origin, df){
               df = data.frame(fromJSON(df))
               return(pie_chart(df,'Type', priority, origin)[[1]][[1]])
             })

app$callback(output = list(id = "cases_by_period", property = 'figure'),
             params = list(
               input(id = "cases_period_dropdown", property = 'value'),
               input(id = "priority_dropdown", property = 'value'),
               input(id = "origin_dropdown", property = 'value'),
               input(id = "cases_df", property = 'children')),
             function(period,priority,origin,df){
               df = data.frame(fromJSON(df))
               return(cases_by_period(df,period,priority,origin)[[1]][[1]])
             })

app$callback(output = list(id = 'cases_by_account', property = 'figure'),
             params = list(
               input(id = 'cases_df', property = 'children')),
             function(df){
               df = data.frame(fromJSON(df))
               return(cases_by_account(df)[[1]])
             })

app$callback(output = list(id = 'cases_modal', property = 'style'),
             params = list(
               input(id = 'new_case', property = 'n_clicks')),
             function(n){
               if (n > 0){
                 return(list('display' = 'block'))
               }
               return(list('display' = 'none'))
             })

app$callback(output = list(id = 'new_case', property = 'n_clicks'),
             params = list(
               input(id = 'submit_new_case', property = "n_clicks")),
             function(n,n2){
               return(0)
             })

app$callback(output = list(id = 'cases_df', property = 'children'),
             params = list(
               input(id = 'submit_new_case', property = 'n_clicks'),
               state("new_case_account", "value"),
               state("new_case_origin", "value"),
               state("new_case_reason", "value"),
               state("new_case_subject", "value"),
               state("new_case_contact", "value"),
               state("new_case_type", "value"),
               state("new_case_status", "value"),
               state("new_case_description", "value"),
               state("new_case_priority", "value"),
               state("cases_df", "children")),
             function(n_clicks, account_id, origin, reason, subject, contact_id, 
                      case_type, status, description, priority, current_df){
               if(n_clicks > 0){
                 add_cases(account_id, origin, reason, subject, contact_id, case_type, status, description, priority)
                 df = get_cases()
                 return(toJSON(df))
               }
               return(current_df)
             })

#app$run_server()