library(dash)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dplyr)
library(sqldf)
library(plotly)
library(lubridate)
library(salesforcer)
source('SFManager.R')
source('functions.R')

accounts = get_accounts()
contacts = get_contacts()
users = get_users()
cases = get_cases()

source('cases.R')
source('leads.R')
source('opportunties.R')

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

app = Dash$new()

app$layout(htmlDiv(
  list(
    # header
    htmlDiv(list(
      
      htmlSpan("CRM App using Salesforce API", className='app-title'),
      
      htmlDiv(
        htmlImg(src="/assets/dash-logo.png" ,height="100%")
        ,style=list("float"="right","height"="100%"))
    ),
    className="row header"
    ),
    
    # tabs
    htmlDiv(list(
      dccTabs(
        id="tabs",
        style=list("height"="20","verticalAlign"="middle"),
        children = list(
          dccTab(label="Opportunities", value="opportunities_tab"),
          dccTab(label="Leads", value="leads_tab"),
          dccTab(label="Cases", value="cases_tab")
        ),
        value="leads_tab"
      )
      
    ),
    className="row tabs_div"
    ),
    htmlDiv(id="tab_content", className="row", style=list("margin"= "4.5% 4%"))),
  className="row",
  style=list("margin"= "0%")
))

app$callback(output=list(id="tab_content", property="children"),
             params=list(
               input(id="tabs", property="value")),
             function(tabs){
               if (tabs == "opportunities_tab"){
                 return(Opplayout)}
               else if (tabs == "cases_tab"){
                 return(caseslayout)}
               else if (tabs == "leads_tab"){
                 return(leadslayout)}
               else{
                 return(caseslayout)}
             })

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server()}

