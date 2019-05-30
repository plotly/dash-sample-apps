library(dashCoreComponents)
library(dashHtmlComponents)
library(dashR)

appName <- Sys.getenv("DASH_APP_NAME")
pathPrefix <- sprintf("/%s/", appName)

Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
           DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

setwd(sprintf("/app/apps/%s", appName))

source("utils/helperFunctions.R")
source("utils/reusableComponents.R")

dataList <- readRDS("data/data.rds")

app <- Dash$new(name = "DashR SVM Explorer")

app$layout(
  htmlDiv(
    list(
      htmlDiv(
        className = "banner",
        children = list(
          htmlDiv(
            className = "container scalable",
            children = list(
              htmlH2(
                htmlA(
                  "Support Vector Machine (SVM) Explorer",
                  href = "https://github.com/plotly/dash-svm",
                  style = list(
                    textDecoration = "none",
                    color = "inherit"
                  )
                )
              ),
              htmlA(
                htmlImg(
                  src = "https://user-images.githubusercontent.com/37411533/56820587-55d56e80-681a-11e9-91b6-026e1551b338.png"
                  ),
                href = "https://plot.ly/products/dash/"
              )
            ))
        )),
        htmlDiv(
          id = "body",
          className = "container scalable",
          children = list(
            htmlDiv(
              id = "app-container",
              className = "row",
              children = list(
                htmlDiv(
                id = "left-column",
                className = "three columns",
                style = list(
                  minWidth = "24.5%",
                  maxHeight = "calc(100vh - 85px)",
                  overflowY = "auto",
                  overflowX = "hidden"
                ),
                children = list(
                  card(
                    id = "first-card",
                    children = list(
                    namedDropdown(
                      name = "Select Dataset",
                      id = "dropdown-select-dataset",
                      options = list(
                        list(label = "Moons", value = "moons"),
                        list(label = "Linearly Separable", value = "linear"),
                        list(label = "Circles", value = "circles")
                      ),
                      clearable = FALSE,
                      searchable = FALSE,
                      value = "moons"
                    ),
                    namedSlider(
                      name = "Sample Size",
                      id = "slider-dataset-sample-size",
                      min = 100,
                      max = 500,
                      step = 100,
                      marks = as.list(
                        setNames(
                          seq(100, 500, 100),
                          seq(100, 500, 100)
                        )
                      ),
                      value = 300
                    ),
                    namedSlider(
                      name = "Noise Level",
                      id = "slider-dataset-noise-level",
                      min = 0,
                      max = 1,
                      step = 0.1,
                      marks = as.list(
                        setNames(
                          seq(0, 1, 0.2),
                          seq(0, 1, 0.2)
                        )
                      ),
                      value = 0.2
                    )
                  )),
                  card(
                    id = "button-card",
                    children = list(
                      namedSlider(
                        name = "Threshold",
                        id = "slider-threshold",
                        min = 0,
                        max = 1,
                        value = 0.5,
                        step = 0.1
                      ),
                      htmlButton(
                        "Reset Threshold",
                        id = "button-zero-threshold",
                        style = list(display = "table", margin = "0 auto")
                      )
                    )
                  ),
                  card(
                    id = "last-card",
                    children = list(
                    namedDropdown(
                      name = "Kernel",
                      id = "dropdown-svm-parameter-kernel",
                      options = list(
                        list(
                          label = "Radial basis function (RBF)",
                          value = "rbf"
                        ),
                        list(label = "Linear", value = "linear"),
                        list(label = "Polynomial", value = "poly"),
                        list(label = "Sigmoid", value = "sigmoid")
                      ),
                      value = "rbf",
                      clearable = FALSE,
                      searchable = FALSE
                    ),
                    namedSlider(
                      name = "Cost (C)",
                      id = "slider-svm-parameter-C-power",
                      min = -2,
                      max = 4,
                      marks = as.list(setNames(10 ** (-2:4), -2:4)),
                      value = 0
                    ),
                    formattedSlider(
                      id = "slider-svm-parameter-C-coef",
                      min = 1,
                      max = 9,
                      value = 1
                    ),
                    namedSlider(
                      name = "Degree",
                      id = "slider-svm-parameter-degree",
                      min = 2,
                      max = 10,
                      value = 3,
                      step = 1,
                      marks = as.list(setNames(seq(2, 10, 2),
                                               seq(2, 10, 2)))
                    ),
                    namedSlider(
                      name = "Gamma",
                      id = "slider-svm-parameter-gamma-power",
                      min = -5,
                      max = 0,
                      value = -1,
                      marks = as.list(setNames(10 ** (-5:0), -5:0))
                    ),
                    formattedSlider(
                      id = "slider-svm-parameter-gamma-coef",
                      min = 1,
                      max = 9,
                      value = 5
                    ),
                    namedRadioItems(
                      name = "Shrinking",
                      id = "radio-svm-parameter-shrinking",
                      labelStyle = list(
                        marginRight = "7px",
                        display = "inline-block"
                      ),
                      options = list(
                        list(label = " Enabled", value = TRUE),
                        list(label = " Disabled", value = FALSE)
                      ),
                      value = TRUE
                    )
                  )),
                  htmlDiv(
                    dccMarkdown(
                      paste0(
                         "[Click here](https://github.com/plotly/dash-svm)",
                         " to visit the project repo, ",
                         "and learn about how to use the app"
                      )
                    ),
                    style = list(
                      margin = "20px 20px",
                      textAlign = "center",
                      color = "#a5b1cd"
                    )
                  )
                )
              ),
              htmlDiv(
                id = "div-graphs",
                children = dccGraph(
                  id = "main_figure",
                  style = list(display = "none")
                )
              )
            )
          )
        )
      )
    )
  )
)

app$callback(
  output("slider-svm-parameter-gamma-coef", "marks"),
  list(input("slider-svm-parameter-gamma-power", "value")),
  function(power){
    s <- 10 ** power
    as.list(setNames(lapply(seq(1, 10, 2) * s, round, 8), seq(1, 10, 2)))
  }
)

app$callback(
  output("slider-svm-parameter-C-coef", "marks"),
  list(input("slider-svm-parameter-C-power", "value")),
  function(power){
    s <- 10 ** power
    as.list(setNames(lapply(seq(1, 10, 2) * s, round, 8), seq(1, 10, 2)))
  }
)

app$callback(
  output("slider-threshold", "value"),
  list(input("button-zero-threshold", "n_clicks"),
       state("main_figure", "figure")),
  function(n_clicks, figure){
    if (!identical(n_clicks, list(NULL))){
      Z <- unlist(figure$data[[1]]$z)
      value <- - (min(Z) / (max(Z) - min(Z)))
    } else {
        value <- 0.4959986285375595
    }
    return(value)
  }
)

app$callback(
  output("slider-svm-parameter-degree", "disabled"),
  list(input("dropdown-svm-parameter-kernel", "value")),
  function(kernel){
    kernel != "poly"
  }
)

app$callback(
  output("slider-svm-parameter-gamma-coef", "disabled"),
  list(input("dropdown-svm-parameter-kernel", "value")),
  function(kernel){
    !kernel %in% c("rbf", "poly", "sigmoid")
  }
)

app$callback(
  output("slider-svm-parameter-gamma-power", "disabled"),
  list(input("dropdown-svm-parameter-kernel", "value")),
  function(kernel){
    !kernel %in% c("rbf", "poly", "sigmoid")
  }
)

app$callback(
  output("div-graphs", "children"),
  list(
    input("slider-svm-parameter-degree", "value"),
    input("slider-svm-parameter-C-coef", "value"),
    input("slider-svm-parameter-C-power", "value"),
    input("slider-svm-parameter-gamma-coef", "value"),
    input("slider-svm-parameter-gamma-power", "value"),
    input("dropdown-select-dataset", "value"),
    input("slider-dataset-noise-level", "value"),
    input("dropdown-svm-parameter-kernel", "value"),
    input("slider-dataset-sample-size", "value"),
    input("slider-threshold", "value"),
    input("radio-svm-parameter-shrinking", "value")
  ),
  function(degree,
           C_coef,
           C_power,
           gamma_coef,
           gamma_power,
           dataset,
           noise,
           kernel,
           n_points,
           threshold,
           shrinking){
    set.seed(42)
    if (dataset == "moons"){
      noise <- matrix(rnorm( (n_points * 2), sd = noise), ncol = 2)
      dat <- dataList$moons[sample(500, n_points)]
      dat[, c("X1", "X2") := dat[, c("X1", "X2")] + noise]
    } else if (dataset == "circles"){
      noise <- matrix(rnorm( (n_points * 2), sd = noise), ncol = 2)
      dat <- dataList$circles[sample(500, n_points)]
      dat[, c("X1", "X2") := dat[, c("X1", "X2")] + noise]
    } else if (dataset == "linear"){
      noise <- matrix(rnorm( (n_points * 2), sd = noise + 0.05), ncol = 2)
      dat <- dataList$linear[sample(500, n_points)]
      dat[, c("X1", "X2") := dat[, c("X1", "X2")] + noise]
    }
    C <- C_coef * 10 ** C_power
    gamma <- gamma_coef * 10 ** gamma_power
    splitData <- trainTestSplit(dat)
    svm <- runSVM(
      train_data = splitData[["train"]],
      kernel = kernel,
      degree = degree,
      sigma = gamma,
      C = C,
      shrinking = shrinking
    )
    f <- generate_main_plot(splitData, svm, threshold = threshold)
    pieConf <- pieConfusionMatrix(svm, splitData, threshold = threshold)
    ROC <- rocCurve(splitData, svm)
    list(
      htmlDiv(
        className = "six columns",
        id = "svm-graph-container",
        style = list(marginTop = "1rem"),
        children = list(
          dccGraph(
            id = "main_figure",
            figure = f,
            style = list(
              height = "calc(100vh - 90px)",
              margin = "0 1.66rem"
            ),
            config = list(displayModeBar = FALSE)
          )
        )
      ),
      htmlDiv(
        className = "three columns",
        id = "graphs-container",
        style = list(
          minWidth = "22%",
          height = "calc(100vh - 90px)",
          marginTop = "1rem",
          "-moz-user-select" = "none",
          "-webkit-user-select" = "none",
          "-ms-user-select" = "none"
        ),
        children = list(
          dccGraph(
            id = "graph-line-roc-curve",
            style = list(height = "40%"),
            figure = ROC,
            config = list(displayModeBar = FALSE)
          ),
          dccGraph(
            id = "graph-pie-confusion-matrix",
            style = list(height = "55%", marginTop = "5%"),
            figure = pieConf,
            config = list(displayModeBar = FALSE))
        )
      )
    )
  }
)

app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
