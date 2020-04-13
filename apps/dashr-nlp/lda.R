library(topicmodels)
library(tidytext)
library(Rtsne)


build_dtm <- function(df) {
  # create stop list
  stop_words <- c("INC.", "Ltd.", "xxxxxxxx", "xxxxxxxxxxxxxxxx", "XXXX", "XX", "xx", "xxxx", "n't", "Trans Union", "BOA", "Citi", "account")
  stop_words <- tolower(stop_words)
  stop_words <- gsub(pattern = "\\,", "", stop_words)
  stop_words <- gsub(pattern = "\\.", "", stop_words)
  
  docs <- VCorpus(VectorSource(df$complaint))
  docs <- tm_map(docs, stripWhitespace)
  docs <- tm_map(docs, content_transformer(tolower))
  docs <- tm_map(docs, removeNumbers)
  docs <- tm_map(docs, removeWords, stopwords("english"))
  docs <- tm_map(docs, removePunctuation)
  docs <- tm_map(docs, stripWhitespace)
  docs <- tm_map(docs, removeWords, stop_words)
  # Create matrix for WCloud
  dtm <- DocumentTermMatrix(docs)
  
  return(dtm)
}


build_lda_topics <- function(dtm) {
  lda = LDA(dtm, k = 5, method= "Gibbs", control = list(seed = 1234))
  topics <- tidy(lda, matrix = "beta")
  return(topics)
}

tsne_topics <- function(topics) {
  tsne <- Rtsne(topics, theta = 0.99)
  tsne_df <- as.data.frame(tsne$Y)
  return(tsne_df)
}


top_terms <- function(topics) {
  ap_top_terms <- topics %>%
    group_by(topic) %>%
    top_n(3, beta) %>%
    ungroup() %>%
    arrange(topic, -beta)
  return(ap_top_terms)
}

build_lda_df <- function(df) {
  message("Build dtm")
  dtm <- build_dtm(df)
  message("Build LDA topics")
  ap_topics <- build_lda_topics(dtm)
  message("Build tsne")
  tsne_df <- tsne_topics(ap_topics)
  message("Return cbind(ap_topics, tsne_df)")
  return(cbind(ap_topics, tsne_df))
}

