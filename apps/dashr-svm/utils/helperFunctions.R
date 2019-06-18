library(caret)
library(kernlab)
library(data.table)
library(plotly)
library(ROCR)

# split the data into train and test sets
trainTestSplit <- function(data){
  data$Y <- as.factor(data$Y)
  training_idx <- createDataPartition(y = data$Y, p = 0.6)
  train_data <- data[training_idx[[1]]]
  test_data <- data[!training_idx[[1]]]
  return(list(train = train_data, test = test_data))
}

# Fit the model
runSVM <- function(train_data, C = 1, sigma = 0.75, degree = 3,
                   kernel = c("rbf", "linear", "poly", "sigmoid"),
                   shrinking = TRUE){
  if (kernel == "rbf"){
    svm <- ksvm(
      as.matrix(train_data[, c("X1", "X2")]),
      as.matrix(train_data$Y),
      type = "C-svc",
      kernel = "rbfdot",
      C = C,
      kpar = list(sigma = sigma),
      prob.model = TRUE,
      shrinking = shrinking
    )
  } else if (kernel == "linear"){
    svm <- ksvm(
      as.matrix(train_data[, c("X1", "X2")]),
      as.matrix(train_data$Y),
      type = "C-svc",
      kernel = "vanilladot",
      C = C,
      prob.model = TRUE,
      shrinking = shrinking
    )
  } else if (kernel == "poly"){
    svm <- ksvm(
      as.matrix(train_data[, c("X1", "X2")]),
      as.matrix(train_data$Y),
      type = "C-svc",
      kernel = "polydot",
      C = C,
      kpar = list(degree = degree, scale = sigma, offset = 0),
      prob.model = TRUE,
      shrinking = shrinking
    )
  } else if (kernel == "sigmoid"){
    # tanh kernel equivalent to sklearn's sigmoid kernel (with correct settings)
    svm <- ksvm(
      as.matrix(train_data[, c("X1", "X2")]),
      as.matrix(train_data$Y),
      type = "C-svc",
      kernel = "tanh",
      C = C,
      kpar = list(scale = sigma, offset = 0),
      prob.model = TRUE,
      shrinking = shrinking
    )
  }
  return(svm)
}

# Plot the datasets and predictions
generate_main_plot <- function(data, model, threshold){
  train_data <- data[["train"]]
  test_data <- data[["test"]]
  full_data <- rbindlist(data)
  {predict(
    model, 
    as.matrix(train_data[, c("X1", "X2")]), 
    type = "decision" ) > threshold} %>% 
      as.integer() %>% 
      as.factor() -> 
        y_pred_train
  trainConfMatrix <- confusionMatrix(
    data = y_pred_train,
    reference = train_data$Y
  )
  train_accuracy <- trainConfMatrix$overall["Accuracy"]
  {predict(
    model, 
    as.matrix(test_data[, c("X1", "X2")]), 
    type = "decision" ) > threshold} %>% 
      as.integer() %>% 
      as.factor() -> 
        y_pred_test
  testConfMatrix <- confusionMatrix(
    data = y_pred_test,
    reference = test_data$Y
  )
  test_accuracy <- testConfMatrix$overall["Accuracy"]
  rngX1 <- range(full_data$X1) * 1.05
  rngX2 <- range(full_data$X2) * 1.05
  grid <- expand.grid(
    X1 = seq(rngX1[1], rngX1[2], length = 50),
    X2 = seq(rngX2[1], rngX2[2], length = 50)
  )
  grid$decision <- predict(model, grid, type = "decision")[, 1]
  scaled_threshold <- threshold *
                     (max(grid$decision) - min(grid$decision)) +
                     min(grid$decision)
  r <- max(
    abs(scaled_threshold - min(grid$decision)),
    abs(scaled_threshold - max(grid$decision))
  )
  csSeq <- seq(0, 1, length.out = 8)
  colorscale <- list(
    list(csSeq[1], "#ff744c"),
    list(csSeq[2], "#ff906d"),
    list(csSeq[3], "#ffc0a9"),
    list(csSeq[4], "#ffe7dc"),
    list(csSeq[5], "#e4fcff"),
    list(csSeq[6], "#c7feff"),
    list(csSeq[7], "#97f8ff"),
    list(csSeq[8], "#00e6ff")
  )
  contours <- plot_ly(
    data = grid,
    x = ~X1,
    y = ~X2,
    z = ~decision,
    zmin = scaled_threshold - r,
    zmax = scaled_threshold + r,
    hoverinfo = "none",
    showscale = FALSE,
    colorscale = colorscale,
    opacity = 0.9,
    contours = list(showlines = FALSE),
    type = "contour"
    ) %>%
      # add decision boundary
      add_trace(
        x = ~X1,
        y = ~X2,
        z = ~decision,
        line = list(color = "#708090", width = "3"),
        showscale = FALSE,
        hoverinfo = "none",
        name = sprintf("Threshold (%s)", round(scaled_threshold, 3)),
        contours = list(
          type = "constraint",
          showlines = FALSE,
          operation = "=",
          value = scaled_threshold
        ),
        type = "contour"
      )
  contours %>%
    add_markers(
      inherit = FALSE,
      x = train_data[train_data$Y == 0, X1],
      y = train_data[train_data$Y == 0, X2],
      showlegend = FALSE,
      name = "Train Data",
      marker = list(
        size = 10,
        color = "#ff3600"
      )
    ) %>%
    add_markers(
      inherit = FALSE,
      x = train_data[train_data$Y == 1, X1],
      y = train_data[train_data$Y == 1, X2],
      name = sprintf(
        "Train Data (accuracy=%s)",
        round(train_accuracy, 3)
      ),
      marker = list(
        size = 10,
        color = "#008dff"
      )
    ) %>%
    add_markers(
      inherit = FALSE,
      x = test_data[test_data$Y == 0, X1],
      y = test_data[test_data$Y == 0, X2],
      name = sprintf(
        "Test Data (accuracy=%s)",
        round(test_accuracy, 3)
      ),
      marker = list(
        size = 13,
        color = "#ff3600",
        symbol = 5
      )
    ) %>%
    add_markers(
      inherit = FALSE,
      x = test_data[test_data$Y == 1, X1],
      y = test_data[test_data$Y == 1, X2],
      showlegend = FALSE,
      name = "Test Data",
      marker = list(
        size = 13,
        color = "#008dff",
        symbol = 5
      )
    ) %>%
    layout(
      xaxis = list(
        ticks = "",
        title = "",
        showticklabels = FALSE,
        showgrid = FALSE,
        zeroline = FALSE,
        range = c(min(grid$X1), max(grid$X1))
      ),
      yaxis = list(
        ticks = "",
        title = "",
        showticklabels = FALSE,
        showgrid = FALSE,
        zeroline = FALSE,
        range = c(min(grid$X2), max(grid$X2))
      ),
      paper_bgcolor = "#272b38",
      plot_bgcolor = "#272b38",
      hovermode = "closest",
      legend = list(
        x = 0, y = -0.01,
        orientation = "h",
        font = list(color = "#a5b1cd")
      ),
      margin = list(t = 0, l = 0, r = 0, b = 0)
    )
}

# Pie chart of the confusion matrix
pieConfusionMatrix <- function(model, data, threshold){
  test_data <- data[["test"]]
  full_data <- rbindlist(data)
  rngX1 <- range(full_data$X1) * 1.05
  rngX2 <- range(full_data$X2) * 1.05
  grid <- expand.grid(
    X1 = seq(rngX1[1], rngX1[2], length = 50),
    X2 = seq(rngX2[1], rngX2[2], length = 50)
  )
  grid$decision <- predict(model, grid, type = "decision")[, 1]
  scaled_threshold <- threshold *
                     (max(grid$decision) - min(grid$decision)) +
                     min(grid$decision)
  {predict(
    model, 
    as.matrix(test_data[, c("X1", "X2")]), 
    type = "decision" ) > scaled_threshold} %>% 
      as.integer() %>% 
      as.factor() -> 
        y_pred_test
  testConfMatrix <- confusionMatrix(
    data = y_pred_test, reference = test_data$Y
  )
  tabl <- testConfMatrix$table
  tn <- tabl[1, 1]
  tp <- tabl[2, 2]
  fn <- tabl[1, 2]
  fp <- tabl[2, 1]
  values <- list(tp, fn, fp, tn)
  label_text <- c(
    "True Positive",
    "False Negative",
    "False Positive",
    "True Negative"
  )
  labels <- c("TP", "FN", "FP", "TN")
  colors <- c("#00c7e8", "#c7feff", "#e9b19d", "#ff744c")
  plot_ly(
    values = values,
    labels = label_text,
    type = "pie",
    hoverinfo = "label+value+percent",
    textinfo = "text+value",
    text = labels,
    rotation = 90,
    sort = FALSE,
    marker = list(colors = colors)
  ) %>%
    layout(
      title = "<b>Confusion Matrix</b>",
      margin = list(l = 10, r = 10, t = 60, b = 10),
      paper_bgcolor = "#272b38",
      plot_bgcolor = "#272b38",
      font = list(color = "#a5b1cd"),
      legend = list(
        bgcolor = "rgba(255,255,255,0)",
        orientation = "h"
      )
    )
}

# plot the ROC curve
rocCurve <- function(data, model){
  testData <- data$test
  ypredscore <- predict(model, testData[, 1:2], type = "decision")
  pred <- prediction(as.numeric(ypredscore), as.numeric(testData$Y))
  perf <- performance(pred, "tpr", "fpr")
  auc <- unlist(performance(pred, "auc")@y.values)
  # using mean fpr and tpr for final ROC plot
  fpr <- rowMeans(as.data.frame(perf@x.values))
  tpr <- rowMeans(as.data.frame(perf@y.values))
  plot_ly(
    x = fpr,
    y = tpr,
    type = "scatter",
    mode = "lines",
    name = "test Data",
    line = list(color = "#00d0e8")
  ) %>%
    layout(
      title = paste0(
        "<b>", sprintf("ROC Curve (AUC = %s)", round(auc, 3)), "</b>"
      ),
      paper_bgcolor = "#272b38",
      plot_bgcolor = "#272b38",
      xaxis = list(title = "False Positive Rate"),
      yaxis = list(title = "True Positive Rate"),
      legend = list(x = 0, y = 1.05, orientation = "h"),
      font = list(color = "#a5b1cd"),
      margin = list(l = 50, r = 10, t = 55, b = 40)
    )
}
