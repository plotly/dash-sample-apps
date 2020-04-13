library(dplyr)

DATASOURCE <-
  read.csv(
    "https://raw.githubusercontent.com/vildly/dash-sample-apps/dash-nlp/apps/dash-nlp/data/customer_complaints_narrative_sample.csv",
    col.names = c('date', 'complaint', 'company'),
    stringsAsFactors = F
  )
DATASOURCE$date <- as.Date(DATASOURCE$date, "%m/%d/%Y")

# Convert datetime objects to pythonish timestamps (integer seconds)
DATASOURCE$date = as.numeric(as.POSIXct(DATASOURCE$date))



top_N <- function(df, N) {
  top_df <- as.data.frame(df %>%
                            group_by(company) %>%
                            summarise(n = n()) %>%
                            arrange(desc(n)))[1:N,]
  
  top_df$company <- factor(top_df$company, levels = unique(top_df$company)[order(top_df$n, decreasing = TRUE)])
  
  return(top_df)
}


get_sample_data <- function(sample_pct) {
  sample_size <- as.integer(nrow(DATASOURCE) * sample_pct / 100)
  SAMPLE <- sample_n(DATASOURCE, size = sample_size)
  return(SAMPLE)
}

filtered_data <- function(sample_pct, company, after, before) {
  # Params:
  #  sample_pct: data sammple size, %%
  #  company: filter by this company
  #  after, before: select only records between after and before
  SAMPLE <- get_sample_data(sample_pct)
  return(SAMPLE[SAMPLE$company == company & SAMPLE$date >= after & SAMPLE$date <= before, ])
}
