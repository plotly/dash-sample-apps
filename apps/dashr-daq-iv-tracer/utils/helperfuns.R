# ''' Labels for source/measure elements '''
GetSourceLabels <- function(source = "Voltage") {

  if (source == "Voltage") {
    # We source voltage and measure current
    sourceLabel <- "Voltage"
    measureLabel <- "Current"
  } else if (source == "Current") {
    # We source current and measure voltage
    sourceLabel <- "Current"
    measureLabel <- "Voltage"
  }
  return(c(sourceLabel, measureLabel))
}

# ''' Units for source/measure elements '''
GetSourceUnits <- function(source = "Voltage") {

  if (source == "Voltage") {
    # We source voltage and measure current
    sourceUnit <- "V"
    measureUnit <- "A"
  } else if (source == "Current") {
    # We source current and measure voltage
    sourceUnit <- "A"
    measureUnit <- "V"
  }
  return(c(sourceUnit, measureUnit))
}

# '''
#     Helper FUN for generating demo measure values
#     Generates MeasureVals from the srcVal
#     based on srcType either "I" or "V"
# '''
SourceAndMeasure <- function(srcType, srcVal,
                             vOc = 20.5, iSc = 3.45,
                             c1 = 0.000002694, c2 = 0.077976842) {

  srcVal <- srcVal * 2.1
  answer <- rep(0, length(srcVal))

  if (srcType == "I") {
    # Values of the input smaller than the short circuit current
    idxOk <- which(srcVal < iSc)
    answer[idxOk] <- c2 * vOc * log(1 + (1 - srcVal[idxOk] / iSc) / c1)
    return(round(answer, 4))
  } else if (srcType == "V") {
    idxOk <- which(srcVal < vOc)
    answer[idxOk] <- iSc * (1 - c1 * (exp(srcVal[idxOk] / (c2 * vOc)) - 1))
    return(round(answer, 4))
  }
}
