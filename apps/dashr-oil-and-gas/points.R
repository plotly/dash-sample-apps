library(dplyr)
library(sqldf)
options(scipen=999)
library(data.table)

dataset <- list()
columns1 <- c('API.Well.Number', 'Gas.Produced..MCF',
            'Water.Produced..bbl', 'Reporting.Year')

columns2 <- c('API.Well.Number', 'Gas.Produced..MCF',
            'Water.Produced..bbl', 'Oil.Produced..bbl', 'Reporting.Year')

df1 = read.csv('data/Oil_and_Gas_Annual_Production__1985_-_2000.csv.gz')
df1 <- df1[,columns1]

df2 = read.csv('data/Oil_and_Gas_Annual_Production__Beginning_2001.csv.gz')
df2 <- df2[,columns2]

dfIn <- merge(df1,df2, all = TRUE)
dfIn[is.na(dfIn)] <- 0
dfIn <- dfIn[!duplicated(dfIn),]
colnames(dfIn)[1] <- 'API_WellNo'
dfIn <- as.data.table(dfIn)

columns = c('Gas.Produced..MCF', 'Water.Produced..bbl',
           'Oil.Produced..bbl', 'Reporting.Year')

dataset <- setNames(split(dfIn, seq(nrow(dfIn))), dfIn[,"API_WellNo"])
#dataset[[1]]['API.Well.Number']
#saveRDS(dataset, file="points.RData")






