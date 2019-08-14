library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(rjson)

appName <- Sys.getenv("DASH_APP_NAME")

if (!appName == "") {
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("/app/apps/%s", appName))
}

########################################## Reading Data ############################################

# 1st page
dfFundData <- read.csv("data/17534.csv")
dfPerfSummary <- read.csv("data/17530.csv",
                          na.strings = c(NA, ""),
                          stringsAsFactors = FALSE)
dfCalYear <- read.csv("data/17528.csv")
dfPerfPc <- read.csv("data/17532.csv")

# 2nd page
dfFundInfo <- read.csv("data/17544.csv")
dfFundCharacteristics <- read.csv("data/17542.csv")
dfFundFacts <- read.csv("data/17540.csv")
dfBondAllocation <- read.csv("data/17538.csv")

# Plots
absReturnPlot <- fromJSON(file = "data/17553.json")
sectorAllocationPlot <- fromJSON(file = "data/17560.json")
currencyWeightPlot <- fromJSON(file = "data/17555.json")
creditAllocationPlot <- fromJSON(file = "data/17557.json")
####################################################################################################


########################################## Plot Edits ##############################################

absReturnPlot$data[[1]]$line$color <- "#129EFF"
absReturnPlot$layout$plot_bgcolor <- "white"
absReturnPlot$layout$height <- 175
absReturnPlot$layout$font$family <- "HelveticaNeue"
# Arranging annotation positions
for (i in 1:length(currencyWeightPlot$layout$annotations)) {
  currencyWeightPlot$layout$annotations[[i]]$font$family <- "HelveticaNeue"
}
absReturnPlot$layout$legend$bgcolor <- "#ecf7fd"
absReturnPlot$data[[1]]$line$width <- 2
absReturnPlot$data[[2]]$line$width <- 2

sectorAllocationPlot$data[[1]]$marker$color <- "#119dff"
sectorAllocationPlot$layout$height <- 195
sectorAllocationPlot$layout$margin$b <- 58
sectorAllocationPlot$data[[1]]$x[2] <- "Asset-Back Securities"
sectorAllocationPlot$data[[1]]$x[3] <- "Residential Mortgages"
sectorAllocationPlot$data[[1]]$x[7] <- "Commercial Mortgages"
sectorAllocationPlot$layout$font$family <- "HelveticaNeue"


currencyWeightPlot$data[[1]]$marker$color <- "#119dff"
currencyWeightPlot$layout$height <- 250
currencyWeightPlot$layout$margin$t <- 0
currencyWeightPlot$layout$margin$b <- 20
currencyWeightPlot$layout$font$family <- "HelveticaNeue"
# Arranging annotation positions
for (i in 1:length(currencyWeightPlot$layout$annotations)) {
  a <- 20
  if (i == 1) {
    a <- 30
  } else if (i == 3 | i == 9) {
    a <- 25
  }
  currencyWeightPlot$layout$annotations[[i]]$x <-
    currencyWeightPlot$layout$annotations[[i]]$x - a
  currencyWeightPlot$layout$annotations[[i]]$font$family <- "HelveticaNeue"
}

creditAllocationPlot$data[[1]]$marker$color <- "#c2ebff"
creditAllocationPlot$data[[2]]$marker$color <- "#119dff"
creditAllocationPlot$layout$height <- 195
creditAllocationPlot$layout$margin$b <- 20
creditAllocationPlot$layout$margin$t <- 0
creditAllocationPlot$layout$font$family <- "HelveticaNeue"
creditAllocationPlot$layout$legend$bgcolor <- "#ecf7fd"
creditAllocationPlot$layout$legend$x <- 0.6
creditAllocationPlot$layout$legend$y <- 0.75
# Arranging annotation positions
for (i in 1:length(creditAllocationPlot$layout$annotations)) {
  a <- 5
  creditAllocationPlot$layout$annotations[[i]]$x <-
    creditAllocationPlot$layout$annotations[[i]]$x - a
  creditAllocationPlot$layout$annotations[[i]]$font$family <- "HelveticaNeue"
}
####################################################################################################


########################################## Helper Fun/Object #######################################

# Return a dash definition of an HTML table for dataframe
MakeDashTable <- function(df) {
  table <- lapply(1:nrow(df), function(row) {
    htmlRow <- lapply(1:length(df[row, ]), function(i) {
      htmlTd(df[row, i])
      })
    htmlTr(htmlRow)
  })
  return(table)
}

modifiedPerfTable <- MakeDashTable(dfPerfSummary)

rowToInsert <- list(
  "",
  htmlTr(
    list(
      htmlTd(list()),
      htmlTd(list("Cumulative"), colSpan = 4,
             style = list("text-align" = "center")),
      htmlTd(list("Annualised"), colSpan = 4,
             style = list("text-align" = "center"))
      ),
    style = list("background" = "white", "font-weight" = "1000"),
  )
)
modifiedPerfTable <- c(rowToInsert, modifiedPerfTable)
####################################################################################################


########################################## App Start ###############################################

app <- Dash$new(name = "Dash Monthly Fund Report")

app$layout(htmlDiv(list(

  htmlDiv(list( # Page 1

    htmlDiv(list( # Subpage 1

      # Header
      htmlDiv(
        className = "row",
        children = htmlImg(
          className = "logo", src = "assets/logo.png"
        )
      ),
      htmlDiv(
        className = "row",
        children = list(
          htmlDiv(
            className = "nine columns header-left",
            children = list(
              htmlH1(
                "Dash Financial Absolute Return Bond II Portfolio"
              ),
              htmlH2("A sub-fund of Dash Monthly Fund, SICAV")
            )
          ),
          htmlDiv(
            className = "three columns header-right",
            children = list(
              htmlH1(
                children = list(
                  htmlSpan("03", className = "light-blue"),
                  htmlSpan("17", className = "blue")
                )
              ),
              htmlH6("Monthly Fund Update")
            )
          )
        )
      ),
      htmlBr(),
      # Investor Profile
      htmlDiv(
        className = "spec-row",
        children = list(
          htmlDiv(
            className = "five columns div-investor-profile",
            children = list(
              htmlH5("Investor Profile"),
              htmlH6("Investor objective"),
              htmlP("Capital appreciation and income."),
              htmlH6("Position in your overall investment portfolio"),
              htmlP("The fund can complement your portfolio.")
            )
          ),
          htmlDiv(
            className = "seven columns div-fund-designed-for",
            children = list(
              htmlH4("The fund is designed for:"),
              htmlBr(),
              htmlP("The fund is designed for investors who are looking
                    for a flexible global investment and sub-investment
                    grade fixed income portfolio that has the ability
                    to alter its exposure with an emphasis on interest
                    rates, currencies and credit markets and that seeks
                    to generate returns through different market conditions
                    with a riskier investment strategy than GS Strategic
                    Absolute Return Bond I Portfolio.")
            )
          )
        )
      ),
      htmlBr(),
      # Content Left
      htmlDiv(
        list(
          htmlDiv(
            list(
              htmlH6(
                "Performance (%)",
                className = "gs-header"
              ),
              htmlTable(MakeDashTable(dfPerfPc)),
              htmlH6(
                "Fund Data",
                className = "gs-header"),
              htmlTable(MakeDashTable(dfFundData), className = "table-tiny")
            ), className = "five columns"
          ),
          # Content Right
          htmlDiv(
            list(
              htmlH6(
                list("Performance (Indexed)"),
                className = "gs-header title-perf-ind"
              ),
              dccGraph(figure = absReturnPlot,
                       style = list("padding-left" = "20px",
                                    "border" = "0",
                                    "width" = "100%",
                                    "height" = "25%")
              ),
              htmlP("This is an actively managed fund that is not designed
                    to track its reference benchmark. Therefore the performance
                    of the fund and the performance of its reference benchmark
                    may diverge. In addition stated reference benchmark returns
                    do not reflect any management  or other charges to the fund,
                    whereas stated returns of the fund do.",
                    className = "title-perf-ind"),
              htmlP("*Past performance does not guarantee future results,
                     which may vary. The value of investments and the income
                     derived from investments will fluctuate and can go down
                     as well as up. A loss of capital may occur.",
                     className = "title-perf-ind")
            ), className = "seven columns"
          )
        ), className = "row"
      ),
      # Performance Summary & Calendar Performance
      htmlDiv(
        list(
              htmlH6(
                "Performance Summary (%)",
                className = "gs-header"),
              htmlTable(modifiedPerfTable, className = "table-tiny"),
              htmlH6(
                "Calendar Year Performance (%)",
                className = "gs-header"
              ),
              htmlTable(MakeDashTable(dfCalYear), className = "table-tiny")
        ), className = "row"
      )
    ), className = "subpage")
  ), className = "page"),

  htmlDiv(list( # Page 2

    htmlDiv(list( # Subpage 2

      # Header
      htmlDiv(
        className = "row",
        children = htmlImg(
          className = "logo", src = "assets/logo.png"
        )
      ),
      htmlDiv(
        className = "row",
        children = list(
          htmlDiv(
            className = "nine columns header-left",
            children = list(
              htmlH1(
                "Dash Financial Absolute Return Bond II Portfolio"
              ),
              htmlH2("A sub-fund of Dash Monthly Fund, SICAV")
            )
          ),
          htmlDiv(
            className = "three columns header-right",
            children = list(
              htmlH1(
                children = list(
                  htmlSpan("03", className = "light-blue"),
                  htmlSpan("17", className = "blue")
                )
              ),
              htmlH6("Monthly Fund Update")
            )
          )
        )
      ),
      htmlBr(),
        # Content pg2 left
        htmlDiv(
          list(
            htmlDiv(
              list(
                htmlH6(
                  "Financial Information",
                  className = "gs-header" ),
                htmlTable(MakeDashTable(dfFundInfo)),
                htmlH6(
                  "Fund Characteristics",
                  className = "gs-header"),
                htmlTable(MakeDashTable(dfFundCharacteristics)),
                htmlH6(
                  "Fund Facts",
                  className = "gs-header"),
                htmlTable(MakeDashTable(dfFundFacts)),
                htmlH6(
                  "Country Bond Allocation (%)",
                  className = "gs-header"),
                htmlTable(MakeDashTable(dfBondAllocation))

              ), className = "six columns"
            ),
            # Content pg2 right
            htmlDiv(
              list(
                htmlH6(
                  "Sector Allocation (%)",
                  className = "gs-header title-pg2-right"),
                dccGraph(figure = sectorAllocationPlot,
                         style = list("padding-left" = "50px",
                                      "border" = "0",
                                      "width" = "100%")),
                htmlH6(
                  "Top 10 Currency Weights (%)",
                  className = "gs-header title-pg2-right"),
                dccGraph(figure = currencyWeightPlot,
                         style = list("padding-left" = "50px",
                                      "border" = "0",
                                      "width" = "300px"
                                      )),
                htmlH6(
                  "Credit Allocation (%)",
                  className = "gs-header title-pg2-right"),
                dccGraph(figure = creditAllocationPlot,
                         style = list("padding-left" = "50px",
                                      "border" = "0",
                                      "width" = "100%"))
              ), className = "six columns"
            )
          ), className = "row"
        )), className = "subpage"
      )
    ), className = "page"
  )
)))
####################################################################################################

if (!appName == ""){
  app$run_server(host = "0.0.0.0", port = Sys.getenv("PORT", 8050))
} else {
  app$run_server(debug = TRUE)
}
