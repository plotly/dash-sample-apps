library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(rjson) #

setwd("/Users/Caner/Desktop/plotly/dashR-monthly-fund-report")

# Reading data
# 1st page
dfFundData <- read.csv("data/17534.csv")
dfPerfSummary <- read.csv("data/17530.csv", na.strings = c(NA, ""), stringsAsFactors = FALSE)
dfCalYear <- read.csv("data/17528.csv")
dfPerfPc <- read.csv("data/17532.csv")

# 2nd page
dfFundInfo <- read.csv("data/17544.csv")
dfFundCharacteristics <- read.csv("data/17542.csv")
dfFundFacts <- read.csv("data/17540.csv")
dfBondAllocation <- read.csv("data/17538.csv")

# Plot Edits
absReturnPlot <- fromJSON(file = "data/17553.json")
absReturnPlot$data[[1]]$line$color <- '#129EFF'
absReturnPlot$layout$plot_bgcolor <- "white"
absReturnPlot$layout$height <- 195
absReturnPlot$layout$font$family <- "HelveticaNeue"
# Arranging annotation positions
for (i in 1:length(currencyWeightPlot$layout$annotations)) {
  currencyWeightPlot$layout$annotations[[i]]$font$family <- "HelveticaNeue"
}
absReturnPlot$layout$legend$bgcolor <- "#ecf7fd"
absReturnPlot$data[[1]]$line$width <- 2
absReturnPlot$data[[2]]$line$width <- 2

sectorAllocationPlot <- fromJSON(file = "data/17560.json")
sectorAllocationPlot$data[[1]]$marker$color <- "#119dff"
sectorAllocationPlot$layout$height <- 195
sectorAllocationPlot$layout$margin$b <- 58
sectorAllocationPlot$data[[1]]$x[2] <- "Asset-Back Securities"
sectorAllocationPlot$data[[1]]$x[3] <- "Residential Mortgages"
sectorAllocationPlot$data[[1]]$x[7] <- "Commercial Mortgages"
sectorAllocationPlot$layout$font$family <- "HelveticaNeue"


currencyWeightPlot <- fromJSON(file = "data/17555.json")
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

creditAllocationPlot <- fromJSON(file = "data/17557.json")
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


# """ Return a dash definition of an HTML table for dataframe """
MakeDashTable <- function(df) {
  table <- lapply(1:nrow(df), function(row) {
    htmlRow <- lapply(1:length(df[row, ]), function(i) {htmlTd(df[row, i])})
    htmlTr(htmlRow)
  })
  return(table)
}

modifiedPerfTable <- MakeDashTable(dfPerfSummary)

rowToInsert <- list(
  "", # -> WHY THIS IS NECESSARY?
  htmlTr(
    list(
      htmlTd(list()),
      htmlTd(list("Cumulative"), colSpan = 4, style = list("text-align" = "center")),
      htmlTd(list("Annualised"), colSpan = 4, style = list("text-align" = "center"))
      ),
    style = list("background" = "white", "font-weight" = "600"),
  )
)
modifiedPerfTable <- c(rowToInsert, modifiedPerfTable)


app <- Dash$new(name = "Dash Monthly Fund Report")

app$layout(htmlDiv(list(

  htmlDiv(list( # Page 1

    htmlA(list("Print PDF"),
          className = "button no-print"
    ),

    htmlDiv(list( # Subpage 1

      # Header
      htmlDiv( # -> header1
        className = "row",
        children = htmlImg(
          className = "logo", src = "assets/logo.png"
        )
      ), # -> //header1
      htmlDiv(
        className = "row",
        children = list(
          htmlDiv(
            className = "nine columns header-left",
            children = list(
              htmlH1(
                "Dash Monthly Fund Absolute Return Bond II Portfolio"
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
      # Row 2
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
      # Row 2.5
      htmlDiv(
        list(
          htmlDiv(
            list(
              htmlH6(
                "Performance (%)",
                className = "gs-header"#gs-text-header padded
              ),
              htmlTable(MakeDashTable(dfPerfPc)),
              htmlH6(
                "Fund Data",
                className = "gs-header"),#gs-text-header padded
              htmlTable(MakeDashTable(dfFundData), className = "tiny-table")
            ), className = "five columns"
          ),
          htmlDiv(
            list(
              htmlH6(
                list("Performance (Indexed)"),
                className = "gs-header title-perf-ind"#gs-table-header padded
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
      # Row 3
      htmlDiv(
        list(
          htmlDiv(
            list(
              htmlH6(
                "Performance Summary (%)",
                className = "gs-header"),#gs-table-header padded
              htmlTable(modifiedPerfTable, className = "table-tiny")
            )
          ),
          htmlDiv(
            list(
              htmlH6(
                "Calendar Year Performance (%)",
                className = "gs-header"#gs-table-header padded
              ),
              htmlTable(MakeDashTable(dfCalYear))
            )
          )
        ), className = "row"
      )



    ), className = "subpage") # // Subpage 1


  ), className = "page"), # // Page 1

  htmlDiv(list( # Page 2
    htmlA(list("Print PDF"),
          className = "button no-print"
    ),

    htmlDiv(list( # Subpage 2

      # Header
      htmlDiv( # -> header2
        className = "row",
        children = htmlImg(
          className = "logo", src = "assets/logo.png"
        )
      ), # -> //header2
      htmlDiv(
        className = "row",
        children = list(
          htmlDiv(
            className = "nine columns header-left",
            children = list(
              htmlH1(
                "Dash Monthly Fund Absolute Return Bond II Portfolio"
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
        # Row 2
        htmlDiv(
          list(
            htmlDiv(
              list(
                htmlH6(
                  "Financial Information",
                  className = "gs-header" ),#gs-text-header padded"
                htmlTable(MakeDashTable(dfFundInfo)),
                htmlH6(
                  "Fund Characteristics",
                  className = "gs-header"), #gs-text-header padded
                htmlTable(MakeDashTable(dfFundCharacteristics)),
                htmlH6(
                  "Fund Facts",
                  className = "gs-header"), #gs-text-header padded
                htmlTable(MakeDashTable(dfFundFacts)),
                htmlH6(
                  "Country Bond Allocation (%)",
                  className = "gs-header"),#gs-table-header padded
                htmlTable(MakeDashTable(dfBondAllocation))

              ), className="six columns"
            ),
            # Column 2
            htmlDiv(
              list(
                htmlH6(
                  "Sector Allocation (%)",
                  className = "gs-header title-pg2-right"),#gs-table-header padded
                dccGraph(figure = sectorAllocationPlot,
                         style = list("padding-left" = "50px",
                                      "border" = "0",
                                      "width" = "100%",
                                      "height" = "30%")),
                htmlH6(
                  "Top 10 Currency Weights (%)",
                  className="gs-header title-pg2-right"), #gs-table-header padded
                dccGraph(figure = currencyWeightPlot,
                         style = list("padding-left" = "50px",
                                      "border" = "0",
                                      "width" = "300px"
                                      )),
                #htmlBr(),
                htmlH6(
                  "Credit Allocation (%)",
                  className = "gs-header title-pg2-right"),#gs-table-header padded
                dccGraph(figure = creditAllocationPlot,
                         style = list("padding-left" = "50px",
                                      "border" = "0",
                                      "width" = "100%",
                                      "height" = "30%"))
              ), className="six columns"
            )

          ), className = "row"
        )



        ), className = "subpage"
      ) # -> subpage 2 end


    ), className = "page"
  ) # -> # Page 2 end


)))

app$run_server(port = 8896)




