library(tidyr)
library(data.table)
library(pracma)

GenerateParentList <- function(time_points) {
# Generates flexible # of columns for the data-table

  parentList <- list(list())
  parentList[[1]] <- list(
    "name" = "Time (hr)",
    "id" = "time",
    "type" = "numeric")

  for (subject in seq(0, as.numeric(time_points) - 1)) {

    parentList[[length(parentList) + 1]] <- parentList[[length(parentList)]]
    # Extend List length by copying last element

    parentList[[length(parentList)]]$name <-
      sprintf("Subj%s Conc (uM)", subject + 1)
    parentList[[length(parentList)]]$id <- subject
    parentList[[length(parentList)]]$type <- "numeric"
    # Assign the values
  }
  return(parentList)
}

PkData2DashT <- function(df) {
  # Converts raw df to nested list format

  dfReshaped <- reshape(df[, -1], direction = "wide",
                        idvar = "time", timevar = "subject_index")
  # Reshape df to time subj0, sub1, subj2 ... subjn format

  names(dfReshaped) <- gsub("conc.", "", names(dfReshaped))

  return(df_to_list(dfReshaped))
  # df_to_list is defined in dashTable
}

GenerateResultColumn <- function(results) {
  # Generates flexible # of columns for the results-table

  resultsRaw <- rbindlist(results)
  # Convert results to df
  resultsDf <- as.data.frame(sapply(
    resultsRaw, function(x) as.numeric(as.character(x))))
  # Convert non-numerics to NA
  resultsDf <- resultsDf[, colSums(is.na(resultsDf)) < nrow(resultsDf)]
  # Remove Column if all NA

  resultList <- list(
    list(
      name = "Parameter",
      id = "parameter",
      type = "character"
    )
  )

  for (subject in as.numeric(names(resultsDf)[1:(ncol(resultsDf) - 1)])) {

    resultList[[length(resultList) + 1]] <- resultList[[length(resultList)]]
    # Extend List length by copying last element

    resultList[[length(resultList)]]$name <- sprintf("Subj%s", subject + 1)
    resultList[[length(resultList)]]$id <- subject
    resultList[[length(resultList)]]$type <- "numeric"
    # Assign the values
  }
  resultList[[length(resultList) + 1]] <- resultList[[length(resultList)]]
  resultList[[length(resultList)]]$name <- "Mean"
  resultList[[length(resultList)]]$id <- "mean"
  resultList[[length(resultList)]]$type <- "numeric"
  # Append Mean column

  resultList[[length(resultList) + 1]] <- resultList[[length(resultList)]]
  resultList[[length(resultList)]]$name <- "StDev"
  resultList[[length(resultList)]]$id <- "stdev"
  resultList[[length(resultList)]]$type <- "numeric"
  # Append StDev column

  return(resultList)
}

InsertAt <- function(a, pos, ...) {
# Function for inserting first position of vector
# Helper FUN #1 for GenerateResultData
  dots <- list(...)
  stopifnot(length(dots) == length(pos))
  result <- vector("list", 2 * length(pos) + 1)
  result[c(TRUE, FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos + 1)))
  result[c(FALSE, TRUE)] <- dots
  return(unlist(result))
}


CalcPk <- function(time, conc, ivCalc = FALSE, termPoints = 3) {
# Returns 6 Pk variables
# Helper FUN #2 for GenerateResultData
  if (min(time, na.rm = TRUE) > 0 && !ivCalc){
  # When doesn't exist insert 0's to first position

    time <- InsertAt(time, 0, 0)
    time <- time[ceiling(length(time) / 2):length(time)]

    conc <- InsertAt(conc, 0, 0)
    conc <- conc[ceiling(length(conc) / 2):length(conc)]
  }

  binded <- as.data.frame(cbind(time, conc))
  binded <- arrange(binded, time)
  # Combine time and conc and sort as df

  time <- binded[complete.cases(binded), "time"]
  conc <- binded[complete.cases(binded), "conc"]
  # Remove NA rows

  maxFilter <- binded$conc == max(binded$conc, na.rm = TRUE)
  maxFilter[which(is.na(maxFilter))] <- FALSE
  # Get the position of max element
  cMax <- binded[maxFilter, "conc"][1]
  # Get max conc if multiple get 1st
  tMax <- binded[maxFilter, "time"][1]
  # Get max time if multiple get 1st

  checks <- length(time) == length(conc) && !(NA %in% time) && !(NA %in% conc)
  # Checks for avoiding errors

  if (checks){
    auc0t <- trapz(time, conc)
  } else {
    auc0t <- NA
  }

  if (length(conc) > termPoints && checks) {
  # Assign variables if checks pass

    polyResult <- polyfit(time[(length(time) - termPoints + 1):length(time)],
      log(conc[(length(conc) - termPoints + 1):length(conc)]), 1)

    slope <- polyResult[1]

    rateConst <- -slope

    timeHalf <- log(2) / rateConst

    aucZeroInf <- auc0t + conc[length(conc)] / rateConst

    percentExtrap <- 100 * (conc[length(conc)] / rateConst) / aucZeroInf
  }else{
    # Handling the case where # of observations are too small or checks fail
    polyResult <- c(NA, NA)
    slope <- NA
    rateConst <- NA
    timeHalf <- NA
    aucZeroInf <- NA
    percentExtrap <- NA
  }
  return(c(timeHalf, auc0t, aucZeroInf, percentExtrap, cMax, tMax))
}

GenerateResultData <- function(records) {
# Generates results-table data

  dfRaw <- rbindlist(records)
  # Create df from records

  df <- as.data.frame(sapply(dfRaw, function(x) as.numeric(as.character(x))))
  # Convert anything non-numeric to NA

  df <- df[, colSums(is.na(df)) < nrow(df)]
  # Remove column if all NA

  parameters <- c("T½ (hr)",  "AUC_0-t (uM*hr)",
                  "AUC_0-inf (uM*hr)", "%Extrap",
                  "Cmax (uM)", "Tmax (hr)")
  # Set values for parameters column

  dfResults <- data.frame(parameter = parameters)
  # Initiate df with parameters

  for (col in select(df, -time)){
    calculatedParameters <- round(CalcPk(df$time, col), 1)
    dfResults <- cbind(dfResults, calculatedParameters)
  }
  # Assign parameter values for each subject

  names(dfResults)[2:ncol(dfResults)] <-
    names(select(df, -time))[1:(ncol(dfResults) - 1)]
  # Set names for subjects

  dfMean <- mutate(dfResults, mean = round(
    rowMeans(select(dfResults, -parameter)), 1))
  dfMeanSDev <- mutate(dfMean, stdev = round(
    apply(select(dfResults, -parameter), 1, sd), 2))
  # Mutate mean & sd columns

  return(df_to_list(dfMeanSDev))
}