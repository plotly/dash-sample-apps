
appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  setwd(sprintf("/app/apps/%s", appName))
}

library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dplyr)


##################### GLOBAL OBJECTS ###############################################

addButton <- htmlButton('Add', id = 'add-button')

delButton <- htmlButton('Delete', id = 'del-button')

titleInput <- dccInput(
  placeholder = 'Enter New Title Here',
  type = 'text',
  value = ''
)

trace <- dccInput(
  placeholder = 'Trace #',
  type = 'text',
  value = ''
)

trace1234 <- dccInput(
  placeholder = '1,2,3,4',
  type = 'text',
  value = ''
)

trace3333 <- dccInput(
  placeholder = '3,3,3,3',
  type = 'text',
  value = ''
)

styleOptions <-list(
  list(label = "0", value = "0"),
  list(label = "circle", value = "circle"),
  list(label = "100", value = "100"),
  list(label = "circle-open", value = "circle-open"),
  list(label = "200", value = "200"),
  
  
  
  #*****************
  list(label = "")
  list(label = "220", value = "220"),
  list(label = "star-triangle-down-dot" , value = "star-triangle-down-dot"),
  list(label = "320", value = "320"),
  list(label = "star-triangle-down-open-dot", value = "star-triangle-down-open-dot"),
  list(label = "21", value = "21"),
  list(label = "star-square", value = "star-square"),
  list(label = "121", value = "121"),
  list(label = "star-square-open", value = "star-square-open"),
  list(label = "221", value = "221"),
  list(label = "star-square-dot", value = "star-square-dot"),
  list(label = "321", value = "321"),
  list(label = "star-square-open-dot", value = "star-square-open-dot"),
  list(label = "22", value = "22"),
  list(label = "star-diamond", value = "star-diamond"),
  list(label = "122", value = "122"),
  list(label = "star-diamond-open", value = "star-diamond-open"),
  list(label = "222", value = "322"),
  list(label = "star-diamond-dot", value = "star-diamond-dot"),
  list(label = "322", value = "322"),
  list(label = "star-diamond-open-dot", value = "star-diamond-open-dot"),
  list(label = "23", value = "23"),
  list(label = "diamond-tall", value = "diamond-tall"),
  list(label = "123", value = "123"),
  list(label = "diamond-tall-open", value = "diamond-tall-open"),
  list(label = "223", value = "223"),
  list(label = "diamond-tall-dot", value = "diamond-tall-dot"),
  list(label = "323", value = "323"),
  list(label = "diamond-tall-open-dot", value = "diamond-tall-open-dot",
  list(label = "24", value = "24"),
  list(label = "diamond-wide", value = "diamond-wide"),
  list(label = "124", value = "124"),
  list(label = "diamond-wide-open", value = "diamond-wide-open"),
  list(label = "224", value = "224"),
  list(label = "diamond-wide-dot", value = "diamond-wide-dot"),
  list(label = "324", value = "324"),
  list(label = "diamond-wide-open-dot", value = "diamond-wide-open-dot"),
  list(label = "25", value = "25"),
  list(label = "hourglass", value = "hourglass"),
  list(label = "125", value = "125"),
  list(label = "hourglass-open", value = "hourglass-open"),
  list(label = "26", value = "26"),
  list(label = "bowtie", value = "bowtie"),
  list(label = "126", value = "126"),
  list(label = "bowtie-open", value = "bowtie-open"),
  list(label = "27", value = "27"),
  list(label = "circle-cross", value = "circle-cross"),
  list(label = "127", value = "127"),
  list(label = "circle-cross-open", value = "circle-cross-open"),
  list(label = "28", value = "28"),
  list(label = "circle-x", value = "circle-x"),
  list(label = "128", value = "128"),
  list(label = "circle-x-open", value = "circle-x-open"),
  list(label = "29", value = "29"),
  list(label = "square-cross", value = "square-cross"),
  list(label = "129", value = "129"),
  list(label = "square-cross-open", value = "square-cross-open"),
  list(label = "30", value = "30"),
  list(label = "square-x", value = "square-x"),
  list(label = "130", value = "130"),
  list(label = "square-x-open", value = "square-x-open"),
  list(label = "31", value = "31"),
  list(label = "diamond-cross", value = "diamond-cross"),
  list(label = "131", value = "131"),
  list(label = "diamond-cross-open", value = "diamond-cross-open"),
  list(label = "32", value = "32"),
  list(label = "diamond-x", value = "diamond-x"),
  list(label = "132", value = "132"),
  list(label = "diamond-x-open", value = "diamond-x-open"),
  list(label = "33", value = "33"),
  list(label = "cross-thin", value = "cross-thin"),
  list(label = "133", value = "133"),
  list(label = "cross-thin-open", value = "cross-thin-open"),
  list(label = "34", value = "34"),
  list(label = "x-thin", value = "x-thin"),
  list(label = "134", value = "134"),
  list(label = "x-thin-open", value = "x-thin-open"),
  list(label = "35", value = "35"),
  list(label = "asterisk", value = "asterisk"),
  list(label = "135", value = "135"),
  list(label = "asterisk-open", value = "asterisk-open"),
  list(label = "36", value = "36"),
  list(label = "hash", value = "hash"), 
  list(label = "136", value = "136"),
  list(label = "hash-open", value = "hash-open"),
  list(label = "236", value = "236"),
  list(label = "hash-dot", value = "hash-dot"),
  list(label = "336", value = "336"),
  list(label = "hash-open-dot", value = "hash-open-dot"),
  list(label = "137", value = "137"),
  list(label = "y-up-open", value = "y-up-open"),
  list(label = "38", value = "38"),
  list(label = "y-down", value = "y-down"),
  list(label = "138", value = "138"),
  list(label = "y-down-open", value = "y-down-open"),
  list(label = "39", value = "39"),
  list(label = "y-left", value = "y-left"),
  list(label = "139", value = "139"),
  list(label = "y-left-open", value = "y-left-open"),
  list(label = "40", value = "40"),
  list(label = "y-right", value = "y-right"),
  list(label = "140", value = "140"),
  list(label = "y-right-open", value = "y-right-open"),
  list(label = "41", value = "41"),
  list(label = "line-ew", value = "line-ew"),
  list(label = "141", value = "141"),
  list(label = "line-ew-open", value = "line-ew-open"),
  list(label = "42", value = "42"),
  list(label = "line-ns", value = "line-ns"),
  list(label = "142", value = "142"),
  list(label = "line-ns-open", value = "line-ns-open"),
  list(label = "43", value = "43"),
  list(label = "line-ne", value = "line-ne"),
  list(label = "143", value = "143"),
  list(label = "line-ne-open", value = "line-ne-open"),
  list(label = "44", value = "44"),
  list(label = "line-nw", value = "line-nw"),
  list(label = "144", value = "144"),
  list(label = "line-nw-open", value = "line-nw-open")
)

#############################CREATE LAYOUT VARIABLES################################

plotlyLogo <- htmlA(
  list(
    htmlImg(src = "assets/dash-logo-new.png",
            className = "logo")),
  href = "https://dashr-docs.herokuapp.com/")

styleDropDown <- dccDropdown(id = "style-dropdown",
                             options = styleOptions,
                             multi = TRUE,
                             value = "circle")

####################################################################################################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
