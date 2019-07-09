library(GEOquery)


#  Helper function to simplify data import for Speck
importSpeck <- function(filepath, 
                        header = FALSE, 
                        skip = 2) {
  textdata <- read.table(
    text = paste0(
      readLines(filepath), collapse="\n"
    ),
    header = header,
    skip = skip,
    col.names = c("symbol", "x", "y", "z"),
    stringsAsFactors = FALSE)
  return(dashTable::df_to_list(textdata))
}



#  Helper function to simplify data import for Clustergram SOFT Files

importSOFT <- function(filepath) {
  geo_data = getGEO(filename = filepath)
  
  geo_table <- Table(geo_data)
  
  row.names(geo_table) <- geo_table$ID_REF
  
  geo_table[1] <- NULL
  
  return(geo_table)
}


# app <- Dash$new()
# 
# app$layout(htmlDiv(list(
#   dccGraph(figure = heatmaply(n[1:10,])
#   )
# )))
# 
# 
# app$run_server(showcase = TRUE)

