library(tidyr)
library(data.table)
library(pracma)

# Generates flexible # of columns for the data-table
GenerateParentList <- function(time_points) {

  parentList <- list(list())
  parentList[[1]] <- list(
    "name" = "Time (hr)",
    "id" = "time",
    "type" = "numeric"
  )

  for (subject in seq(0, as.numeric(time_points) - 1)) {

    # Extend list length by copying last element
    parentList[[length(parentList) + 1]] <- parentList[[length(parentList)]]

    # Assign the values
    parentList[[length(parentList)]]$name <-
      sprintf("Subj%s Conc (uM)", subject + 1)
    parentList[[length(parentList)]]$id <- subject
    parentList[[length(parentList)]]$type <- "numeric"

  }

  return(parentList)
}

# Converts raw df to nested list format
PkData2DashT <- function(df) {

  # Reshape df to time subj0, sub1, subj2 ... subjn format
  dfReshaped <- reshape(df[, -1], direction = "wide",
                        idvar = "time", timevar = "subject_index")

  names(dfReshaped) <- gsub("conc.", "", names(dfReshaped))

  # df_to_list is defined in dashTable
  return(df_to_list(dfReshaped))
}

# Generates flexible # of columns for the results-table
GenerateResultColumn <- function(results) {

  # Convert results to df
  resultsRaw <- rbindlist(results)

  # Convert non-numerics to NA
  resultsDf <- as.data.frame(sapply(
    resultsRaw, function(x) as.numeric(as.character(x))))

  # Remove Column if all NA
  resultsDf <- resultsDf[, colSums(is.na(resultsDf)) < nrow(resultsDf)]

  resultList <- list(
    list(
      name = "Parameter",
      id = "parameter",
      type = "character"
    )
  )

  for (subject in as.numeric(names(resultsDf)[1:(ncol(resultsDf) - 1)])) {
    # Extend List length by copying last element
    resultList[[length(resultList) + 1]] <- resultList[[length(resultList)]]

    # Assign the values
    resultList[[length(resultList)]]$name <- sprintf("Subj%s", subject + 1)
    resultList[[length(resultList)]]$id <- subject
    resultList[[length(resultList)]]$type <- "numeric"
  }

  # Append Mean column
  resultList[[length(resultList) + 1]] <- resultList[[length(resultList)]]
  resultList[[length(resultList)]]$name <- "Mean"
  resultList[[length(resultList)]]$id <- "mean"
  resultList[[length(resultList)]]$type <- "numeric"

  # Append StDev column
  resultList[[length(resultList) + 1]] <- resultList[[length(resultList)]]
  resultList[[length(resultList)]]$name <- "StDev"
  resultList[[length(resultList)]]$id <- "stdev"
  resultList[[length(resultList)]]$type <- "numeric"

  return(resultList)
}

# Function for inserting first position of vector
# Helper FUN #1 for GenerateResultData
InsertAt <- function(a, pos, ...) {

  dots <- list(...)
  stopifnot(length(dots) == length(pos))
  result <- vector("list", 2 * length(pos) + 1)
  result[c(TRUE, FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos + 1)))
  result[c(FALSE, TRUE)] <- dots

  return(unlist(result))
}

# Returns 6 Pk variables
# Helper FUN #2 for GenerateResultData
CalcPk <- function(time, conc, ivCalc = FALSE, termPoints = 3) {

  # When doesn't exist insert 0's to first position
  if (min(time, na.rm = TRUE) > 0 && !ivCalc){
    time <- InsertAt(time, 0, 0)
    time <- time[ceiling(length(time) / 2):length(time)]
    conc <- InsertAt(conc, 0, 0)
    conc <- conc[ceiling(length(conc) / 2):length(conc)]
  }

  # Combine time and conc and sort as df
  binded <- as.data.frame(cbind(time, conc))
  binded <- arrange(binded, time)

  # Remove NA rows
  time <- binded[complete.cases(binded), "time"]
  conc <- binded[complete.cases(binded), "conc"]

  # Get the position of max element
  maxFilter <- binded$conc == max(binded$conc, na.rm = TRUE)
  maxFilter[which(is.na(maxFilter))] <- FALSE
  # Get max conc if multiple get 1st
  cMax <- binded[maxFilter, "conc"][1]
  # Get max time if multiple get 1st
  tMax <- binded[maxFilter, "time"][1]

  # Checks for avoiding errors
  checks <- length(time) == length(conc) && !(NA %in% time) && !(NA %in% conc)

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
  } else {
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

# Generates results-table data
GenerateResultData <- function(records) {

  # Create df from records
  dfRaw <- rbindlist(records)

  # Convert anything non-numeric to NA
  df <- as.data.frame(sapply(dfRaw, function(x) as.numeric(as.character(x))))

  # Remove column if all NA
  df <- df[, colSums(is.na(df)) < nrow(df)]

  # Set values for parameters column
  parameters <- c("T½ (hr)",  "AUC_0-t (uM*hr)",
                  "AUC_0-inf (uM*hr)", "%Extrap",
                  "Cmax (uM)", "Tmax (hr)")

  # Initiate df with parameters
  dfResults <- data.frame(parameter = parameters)

  # Assign parameter values for each subject
  for (col in select(df, -time)){
    calculatedParameters <- round(CalcPk(df$time, col), 1)
    dfResults <- cbind(dfResults, calculatedParameters)
  }

  # Set names for subjects
  names(dfResults)[2:ncol(dfResults)] <-
    names(select(df, -time))[1:(ncol(dfResults) - 1)]

  # Mutate mean & sd columns
  dfMean <- mutate(dfResults, mean = round(
    rowMeans(select(dfResults, -parameter)), 1))
  dfMeanSDev <- mutate(dfMean, stdev = round(
    apply(select(dfResults, -parameter), 1, sd), 2))

  return(df_to_list(dfMeanSDev))
}
