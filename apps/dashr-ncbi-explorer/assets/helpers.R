generate_table <- function(df, nrows=20)
  
  # function generates a dash table from a supplied data frame (df)
  
  #  and number of rows (nrows) to display
{
  
  n <- min(nrows, nrow(df))
  
  rows <- lapply(seq(1, n), function(i) {
    
    htmlTr(children = lapply(as.character(df[i,]), htmlTd))
    
  })
  
  header <- htmlTr(children = lapply(names(df), htmlTh))
  
  htmlTable(
    
    children = c(list(header), rows)
  )
}

results_dataframe <- function(search_term) {
  search_term <- gsub('"', "", search_term)
  
  search <- entrez_search(db = "nuccore", term = search_term, retmax=10)
  
  search_seqs <- entrez_fetch(db = "nuccore", id = search$ids, rettype = "fasta")
  
  red <- unlist(str_extract_all(search_seqs, ">.*\\n"))
  
  titles_vector <- c()
  
  for (i in red) {
    title <- str_extract_all(i, ">.*?[:blank:]")
    title <- substring(title, 2)
    titles_vector <- c(titles_vector, title)
  }
  
  sequence_titles <- as.list(unlist(str_extract_all(search_seqs, ">.*\\n")))
  
  sequence_dataframe <- do.call(rbind.data.frame, sequence_titles)
  
  names(sequence_dataframe)[1] <- "Dataset Accession ID and Description"
  
  sequence_dataframe$Accession_ID <- titles_vector
  
  sequence_dataframe <- sequence_dataframe[, c(2,1)]
  
  return(sequence_dataframe)
}
