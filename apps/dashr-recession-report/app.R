appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)

  setwd(sprintf("/app/apps/%s", appName))
}

library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dplyr)

################################ LOAD DATA & CREATE GLOBAL OBJECTS #################################

dfJobs <- read.csv("data/nyt_255_ces.csv", stringsAsFactors = FALSE,
                   na.strings = c("NA", ""))
dfWages <- read.csv("data/nyt_255_wages.csv", stringsAsFactors = FALSE,
                   na.strings = c("NA", ""))

####################################################################################################

######################################## DATA PREP #################################################

elm <- list()
series <- rep(list(elm), 6)
names(series) <- c("absolute_wages", "absolute_jobs", "relative_jobs",
                   "relative_wages", "meta", "job_growth")
# Initiated nested parent list of length 6

for (i in seq(length(dfJobs$cescode))) {
  jobsMeta <- dfJobs[i, 1:11]
  jobsRow <- dfJobs[i, 38:ncol(dfJobs)]
  # Selecting jobs starting from March 2006 for
  # making  jobs data consistent with wages data.
  name <- jobsMeta$nytlabel
  cescode <- as.character(jobsMeta$cescode)
  # Extracted category names and their codes

  wagesRow <- select(filter(dfWages, seriesid == cescode), -seriesid)
  # Filtered out cescode for all dates

  series$absolute_wages[[cescode]] <- as.list(wagesRow)
  series$absolute_jobs[[cescode]] <- as.list(jobsRow)
  # Lists of wages & jobs for each industry

  series$relative_job[[cescode]] <- 100 * (jobsRow[length(jobsRow)] -
    jobsRow[1]) / jobsRow[length(jobsRow)]
  # % Change between first & last # of jobs
  series$relative_wages[[cescode]] <- 100 * (wagesRow[length(wagesRow)] -
    wagesRow[1]) / wagesRow[length(wagesRow)]
  # % Change between first & last # of jobs

  series$meta[[cescode]] <- jobsMeta
  # All other info in jobs data except # of jobs

  series$job_growth[[cescode]] <- jobsRow[length(jobsRow)] - jobsRow[1]
  # Difference # of jobs btw last and first observed dates
}

####################################################################################################

#################################### MAIN FUNCTION START ###########################################

CreateFigure <- function(highlight_cescode = character(),
                         skip_labels = character(),
                         show_only = character()) {
  # Creates wages vs. # of jobs figure for 255 industries.
  #
  # Args:
  #   highlight_cescode: Vector of job codes to be highlighted
  #   skip_labels: Vector of job labels to hide for highlighted series
  #   show_only: Vector of job labels to show for highlighted series
  #
  # Returns:
  #   List of data & layout for the figure property of dccGraph

  maxWage <- max(dfWages["X2006.03.01"])
  minWage <- min(dfWages["X2006.03.01"])
  # Getting min & max wages for earliest date

  scale <- c("rgb(202,0,32)", "rgb(244,165,130)", "rgb(247,247,247)",
             "rgb(146,197,222)", "rgb(5,113,176)")
  # Colors for lines

  traces <- list(list())
  annotations <- list(list())
  # Initate empty traces & annotations for appending later

  for (cescode in names(series[["absolute_wages"]])) {
  # Main loop for populating figure's data for jobs
    growth <- as.numeric(series$relative_job[[cescode]])
    absJobs <- series$absolute_jobs[[cescode]]
    initialWage <- series$absolute_wages[[cescode]][[1]]
    relativeWage <- (initialWage - minWage) / maxWage
    # Assign necessary list elements to variables

    x <- seq(from = relativeWage, to = relativeWage + 0.2,
             length.out = length(series$absolute_wages[[cescode]]))
    # Wages data

    relativeGrowthAcrossTime <- (unlist(absJobs) -
        absJobs[[1]]) / absJobs[[length(absJobs)]]

    # Relative Growth Time component for y series

    maxJobsCreated <- max(unlist(series$job_growth))
    minJobsCreated <- min(unlist(series$job_growth))
    # min & max to derive Relative Growth  Industry component

    jobs <- series$absolute_jobs[[cescode]]
    # Number of jobs for y series

    relativeGrowthAcrossIndustry <-
    ( (jobs[[length(jobs)]] - jobs[[1]]) - minJobsCreated ) / maxJobsCreated
    # Relative Growth  Industry component for y series

    y <- as.numeric(relativeGrowthAcrossTime +
          (relativeGrowthAcrossIndustry * 2))
    # y series ranging between [-1,4] indicating job growth

    if (growth > 20) {
      color <- scale[5]
      legendgroup <- ">20%"
      name <- "Greater than 20%"
    } else if (growth > 10) {
      color <- scale[4]
      legendgroup <- ">10%"
      name <- "Between 10% and 20%"
    } else if (growth > -10) {
      color <- "lightgrey"
      legendgroup <- ">-10%"
      name <- "Between -10% and 10%"
    } else if (growth > -20) {
      color <- scale[2]
      legendgroup <- ">-20%"
      name <- "Between -20% and -10%"
    } else {
      color <- scale[1]
      legendgroup <- "<-20%"
      name <- "Less than -20%"
    }
    # Set line colors and legend entries based on growth

    if (length(highlight_cescode) > 0 &&
        !(cescode %in% highlight_cescode)) {
    # When highlight code entered but cescode not in highlight
      color <- "lightgrey"
    }

    hoverinfo <- "text"
    # Standard hoverinfo
    width <- 1
    # Standard line width

    if (length(highlight_cescode) > 0 &&
        cescode %in% highlight_cescode) {
      # When highlight code entered and cescode in highlight
      width <- 2.5
      # Thicker line
      if (color == "lightgrey") {
        color <- "grey"
      }
      hoverinfo <- "text"
    } else if (length(highlight_cescode) > 0) {
      width <- 0.5
      # Cescodes not selected thinner
      hoverinfo <- "none"
      # Cescodes not selected no hoverinfo
    }

    label <- series$meta[[cescode]]$nytlabel
    # Assign label to variable

    text <- list()
    for (j in seq(ncol(wagesRow))) {
      text[[j]] <- paste("<b>", label, "</b><br>",
        "Overall job growth: ", round(growth, 1), "%<br>",
        "Wages: $", series$absolute_wages[[cescode]][[j]], " / hour<br>",
        "Jobs: ", jobs[[j]], "k<br>", "Date: ", gsub("X", "", names(jobs)[[j]]),
        sep = "")
    }
    # Prepare the text to be shown when hovering over lines

    ### Legend Start
    if (length(traces) <= 1) {
      legendsShown <- list("dummy")
    # Set legendsShown to a dummy value for the first iteration
    }
    if (length(traces) >= 2) {
    # When some trace is already appended
      for (fig in traces[-1]) {
      # Exluding the first empty trace
        if (!(fig$legendgroup %in% legendsShown)) {
        # When legend is seen for the first time, add it to legendsShown
         legendsShown <- append(legendsShown, fig$legendgroup)
        }
      }
    }
    if (legendgroup %in% legendsShown) {
      showlegend <- FALSE
    } else {
      showlegend <- TRUE
    }
    # Showlegend if not shown before
    ### Legend End

    tracestoAppend <- list(
      "x" = x,
      "y" = y,
      "mode" = "lines",
      "line" = list("color" = color, "width" = width),
      "text" = text,
      "legendgroup" = legendgroup,
      "name" = name,
      "hoverinfo" = hoverinfo,
      "showlegend" = showlegend
    )
    # Collect trace data in list

    traces <- append(traces, list(tracestoAppend))
    # Append current trace data to main traces list

    if (length(highlight_cescode) > 0 && (cescode %in% highlight_cescode) &&
        (length(skip_labels) > 0 && !(label %in% skip_labels)) ||
        (length(show_only) > 0 && (label %in% show_only)) &&
        length(traces) >= 2
        ) {
    # 1- highlight entered, cescode in highlight
    # 2- skip label entered and label is not in skip OR
    # show_only entered and label is in show only
    # 3- there should be at least 1 trace added traces >=2
      annotationstoAppend <- list(
        "x" = traces[[length(traces)]]$x[length(x)],
        "xref" = "x",
        "xanchor" = "left",
        "y" = traces[[length(traces)]]$y[length(y)],
        "yref" = "y",
        "yanchor" = "middle",
        "showarrow" = FALSE,
        "text" = label,
        "font" = list("size" = 12),
        "bgcolor" = "rgba(255, 255, 255, 0.5)"
      )
      # Annotate job industry for highlight series

      annotations <- append(annotations, list(annotationstoAppend))
      # Append annotation to annotations list
    }
  }# End of main for loop

  if (length(traces) > 1) {
  # Removing first element list(list()) (empty nested list)
  # from populated traces
    traces <- traces[-1]
  }
  if (length(highlight_cescode) == 0) {
  # Moving trace to highlight to first place in traces list
    for (i in seq(length(traces))) {
      t <- traces[[i]]
      if (t$showlegend) {
        if (t$legendgroup == ">20%") {
          traceToMove <- traces[[i]]
          traces[[i]] <- NULL
          traces <- append(traces, list(traceToMove), 0)
        } else if (t$legendgroup == ">10%") {
          traceToMove <- traces[[i]]
          traces[[i]] <- NULL
          traces <- append(traces, list(traceToMove), 1)
        } else if (t$legendgroup == ">-10%") {
          traceToMove <- traces[[i]]
          traces[[i]] <- NULL
          traces <- append(traces, list(traceToMove), 2)
        } else if (t$legendgroup == ">-20%") {
          traceToMove <- traces[[i]]
          traces[[i]] <- NULL
          traces <- append(traces, list(traceToMove), 3)
        } else if (t$legendgroup == "<-20%") {
          traceToMove <- traces[[i]]
          traces[[i]] <- NULL
          traces <- append(traces, list(traceToMove), 4)
        }
      }
    }
  } else {
    # Below is done to display thick line at the most outer level
    for (i in seq(length(traces))) {
      t <- traces[[i]]
      if (t$line[["width"]] == 2.5) {
        traces <- append(traces, list(t), length(traces))
      }
    }
  }

  if (length(highlight_cescode) == 0) {
    annotations <- list(
      list(
        "x" = 0.8, "xref" = "paper", "xanchor" = "left",
        "y" = 0.95, "yref" = "paper", "yanchor" = "bottom",
        "text" = "<b>Job Growth</b>",
        "showarrow" = FALSE
      )
    )
    # Set legend title
  }

  layout <- list(
    "xaxis" = list(
      "showgrid" = FALSE,
      "showline" = FALSE,
      "zeroline" = FALSE,
      "showticklabels" = FALSE,
      "ticks" = "",
      "title" = "← Lower Wages        Industries        Higher Wages →"
    ),
    "yaxis" = list(
      "showgrid" = FALSE,
      "showline" = FALSE,
      "zeroline" = FALSE,
      "showticklabels" = FALSE,
      "ticks" = "",
      "title" = "Jobs"
    ),
    "showlegend" = ifelse(length(highlight_cescode) > 0, FALSE, TRUE),
    "hovermode" = "closest",
    "legend" = list(
      "x" = 0.8,
      "y" = 0.95,
      "xanchor" = "left"
    ),
    "annotations" = annotations,
    "margin" = list("t" = 20, "b" = 20, "r" = 0, "l" = 20),
    "font" = list("size" = 12)
    )
  # Layout Created

  return (list("data" = traces, "layout" = layout ))
}

####################################################################################################

################################### INITIATE APPLICATION ###########################################

app <- Dash$new()

####################################################################################################

####################################### CREATE LAYOUT ##############################################

options <- list()
for (c in sort(unique(dfJobs$nytlabel))){
  optionItem <- list(label = c, value = c)
  options <- append(options, list(optionItem))
}
# Dropdown option for user custom selection

app$layout(htmlDiv(
  list(
    dccMarkdown(
      " # How the Recession Reshaped the Economy, in 255 Charts

Five years since the end of the Great Recession,
the economy has finally regained the nine million jobs it lost.
But not all industries recovered equally.
Each line below shows how the number of jobs has changed for
a particular industry over the past 10 years.
Scroll down to see how the recession reshaped the nation’s job market,
industry by industry.

> This interactive report is a rendition of a
> [New York Times original](https://www.nytimes.com/interactive/2014/06/05/upshot/how-the-recession-reshaped-the-economy-in-255-charts.html).
> This app demonstrates how to build high-quality, interactive
> reports using the Dash framework in Python.

***",
      className = "container",
      style = list("maxWidth" = "650px")),
    dccMarkdown(
      " ## A Mixed Recovery

Industries in the health care and energy sectors grew substantially
over the last five years, while jobs in real estate and
construction continued to shrink.
Industries that paid in the middle of the wage spectrum
generally lost jobs. And while the economy overall
is back to its pre-recession level, it hasn't added the
roughly 10 million jobs needed to keep up with growth
in the working-age population.",
      className = "container",
      style = list("maxWidth" = "650px")),
    dccGraph(figure = CreateFigure(),
             id = "overview", style = list("height" = "90vh")),
    dccMarkdown(
      " ***
## More Bad — and Good — Jobs
Americans often lament the quality of jobs today, and some
low-paying industries — like **fast food**, where annual average pay
is less than $22,000 — are growing.
But so are some high-paying sectors, such as **consulting**,
**computing** and **biotech**.",
      className = "container",
      style = list("maxWidth" = "650px")),
    dccGraph(
      figure = CreateFigure(
        as.character(unlist(
          select(filter(dfJobs, alternacategory == "wages"), "cescode"))),
        skip_labels = c(
          "Full-service restaurants",
          "Engineering and drafting services",
          "Gasoline stations",
          "Supermarkets and other grocery stores",
          "Used merchandise stores"
        )
      ), id = "mixed", style = list("height" = "90vh")),
    dccMarkdown(
    "***
## The Medical Economy
The middle-wage industries that have added jobs are
overwhelmingly in health care.
Labs, home-care providers and dentist offices all pay
between $18 and $29 an hour on average —
and all have grown.
But these gains have not offset losses in other middle-wage
industries, such as airlines and construction.",
      className = "container",
      style = list("maxWidth" = "650px")),
    dccGraph(
      figure = CreateFigure(
        as.character(unlist(
          select(filter(dfJobs, alternacategory == "health"), "cescode"))),
        show_only = c(
          "Home health care services",
          "Offices of physicians",
          "General and surgical hospitals",
          "Diagnostic imaging centers",
          "Offices of dentists",
          "Outpatient care centers, except mental health"
        )
      ), id = "health", style = list("height" = "90vh")),
    dccMarkdown(
      "***
## A Long Housing Bust
Home prices have rebounded from their crisis lows,
but home building remains at historically low levels.
Overall, industries connected with construction and
real estate have lost 19 percent of their jobs since
the recession began — hundreds of thousands more than
health care has added.",
      className = "container",
      style = list("maxWidth" = "650px")),
    dccGraph(
      figure = CreateFigure(
        as.character(unlist(
          select(filter(dfJobs, alternacategory == "housing"), "cescode"))),
        show_only = c(
          "General contractors for new homes",
          "Land subdivision",
          "For-sale home builders",
          "Building finishing contractors (drywall, painting)",
          "Interior design services",
          "Architectural services",
          "Wood product manufacturing",
          "For"
        )
      ), id = "real-estate", style = list("height" = "90vh")),
    dccMarkdown(
      "***
## Black Gold Rush
While it took a hit from the recession,
oil and gas extraction — and its associated jobs —
have been booming,
transforming economies in resource-rich places like West Texas
and North Dakota.
Many of these industries have average salaries above $70,000.",
      className = "container",
      style = list("maxWidth" = "650px")),
    dccGraph(
      figure = CreateFigure(
        as.character(unlist(
          select(filter(dfJobs, alternacategory == "oil"), "cescode"))),
        skip_labels = c(
          "Petroleum and coal product manufacturing"
        )
      ), id = "oil", style = list("height" = "90vh")),
    dccMarkdown(
      "***
## Digital Revolution
Bookstores, printers and publishers of newspapers and magazines
have lost a combined 400,000 jobs since the recession began.
Internet publishers — including web-search firms — offset
only a fraction of the losses, adding 76,000 jobs.
Electronic shopping and auctions made up the
fastest-growing industry, tripling in employment in 10 years.",
      className = "container",
      style = list("maxWidth" = "650px")),
    dccGraph(
      figure = CreateFigure(
        as.character(unlist(
          select(filter(dfJobs, alternacategory == "media"), "cescode"))),
        skip_labels = c(
          "Video and photography services",
          "Book publishers (and those that also publish e-books)",
          "Data processing, web hosting and related services"
        )
      ), id = "media", style = list("height" = "90vh")),
    dccMarkdown(
      "***
## Grooming Boom
In the midst of recession, Americans held on to simple luxuries —
for themselves and their pets.
Nail salons, which made up one of the most resilient industries,
were closely rivaled by pet boarding, grooming and training.",
      className = "container",
      style = list("maxWidth" = "650px")),
    dccGraph(
      figure = CreateFigure(
        as.character(unlist(
          select(filter(dfJobs, alternacategory == "booming"), "cescode"))),
        skip_labels = c(
          "Barber shops and beauty salons"
        )
      ), id = "booming", style = list("height" = "90vh")),
    dccMarkdown(
      "***
## And More
Discover patterns yourself by filtering through industries with
the dropdown below.",
      className = "container",
      style = list("maxWidth" = "650px")),
    htmlDiv(
      dccDropdown(
        options = options,
        value = list("Florists", "Full-service restaurants"),
        multi = TRUE,
        id = "category-filter"
      ), className = "container", style = list("maxWidth" = "650px")),
    htmlDiv(children = list(dccGraph()), id = "filtered-content"),
    dccMarkdown(
        "***
## Made with Dash
This report was written in R with the Dash framework.
Dash abstracts away all of the server logic, API calls, CSS,
HTML, and Javascript that is usually required to produce a rich web
application. This application was written in a single R file
containing around 500 lines of code. This includes the data analysis,
markup, and interactive visualizations. See for yourself below.

Interested in what you see?
[Get in touch](https://plot.ly/products/consulting-and-oem/).",
        className = "container",
        style = list("maxWidth" = "650px")),
    htmlDiv(
      dccMarkdown(children = paste("```r", "\n",
        readChar("app.R", file.info("app.R")$size), "\n",
        "```", sep = "")
        ),
      className = "container",
      style = list("maxWidth" = "650px", "borderLeft" = "thin solid lightgrey")
    )
)))

####################################################################################################

########################################## CALLBACK ################################################
app$callback(output = list(id = "filtered-content", property = "children"),
             params = list(input(id = "category-filter", property = "value")),

  function(selected_values) {
    if (length(selected_values) > 0) {
      selectedValues <- unlist(selected_values)
      cescodes <- dfJobs[which(dfJobs$nytlabel %in% selectedValues), "cescode"]
    } else {
      cescodes <- c("-")
      skip_labels <- c("-")
    }
    graph <- dccGraph(
      figure = CreateFigure(
        highlight_cescode = cescodes,
        skip_labels = c("-")))
    return(graph)
  }
)

####################################################################################################

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
