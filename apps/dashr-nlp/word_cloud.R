library(tm)

# remove like xxx...
removePunctWords <- function(x) {
  gsub(pattern = "xx\\w*", "", x)
}

# stoplist
stoplist_funct <- function(bank) {
  a <- strsplit(as.character(bank), " ")
  a <- tolower(unlist(a))
  # create stop list
  stop_words <- c(a, "INC.", "Ltd.", "xxxxxxxx","xxxxxxxxxxxx", "xxxxxxxxxxxxxxxx", "XXXX", "XX", "xx", "xxxx", "n't", "Trans Union", "BOA", "Citi", "account")
  stop_words <- tolower(stop_words)
  stop_words <- gsub(pattern = "\\,", "", stop_words)
  stop_words <- gsub(pattern = "\\.", "", stop_words)
  return(stop_words)
}


word_cloud <- function(words, stoplist) {

  docs <- VCorpus(VectorSource(words))
  # Remove punctuations
  #docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  # Remove your own stop word
  docs <- tm_map(docs, removeWords, stoplist)
  # Create matrix for WCloud
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m), decreasing = TRUE)
  df <- data.frame(word = names(v), n = v, freq = v / sum(v))
  df$word <- removePunctWords(df$word)
  return(df)
}
