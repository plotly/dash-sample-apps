library(dash)
library(dashDaq)
library(dashCoreComponents)
library(dashHtmlComponents)

appName <- Sys.getenv("DASH_APP_NAME")

if (!appName == "") {
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

source("utils/helperfuns.R")

# Define the app
app <- Dash$new(name = "DashR DAQ IV Tracer")

# Font and background colors associated with each theme
bannerColor <- list("dark" = "#23262e", "light" = "#ffffff")
bkgColor <- list("dark" = "#23262e", "light" = "#f6f6f7")
gridColor <- list("dark" = "#53555B", "light" = "#969696")
textColor <- list("dark" = "#95969A", "light" = "#595959")
cardColor <- list("dark" = "#2D3038", "light" = "#FFFFFF")
accentColor <- list("dark" = "#FFD15F", "light" = "#ff9827")

# ''' Generate the layout of the app '''
GenerateMainLayout <- function(theme = "light",
                               srcType = "Voltage",
                               modeVal ="Single measure",
                               fig = NULL,
                               measSrc = 0,
                               measDisplay = 0,
                               srcKnob = 0,
                               sourceToggle = FALSE,
                               modeToggle = FALSE) {

  sourceLabel <- GetSourceLabels(srcType)[1]
  measureLabel <- GetSourceLabels(srcType)[2]
  sourceUnit <- GetSourceUnits(srcType)[1]
  measureUnit <- GetSourceUnits(srcType)[2]

  if (is.null(fig)) {
    data <- list()
  } else {
    data <- fig[["data"]]
  }

  htmlLayout <- list(
    htmlDiv(
      id = "page-body-content",
      className = "row flex-display",
      children = list(
        htmlDiv(
          id = "figure-card",
          className = "six columns",
          style = list("backgroundColor" = cardColor[[theme]]),
          children = list(
            htmlP("IV Curve"),
            dccGraph(
              id = "IV_graph",
              style = list("width" = "100%"),
              figure = list(
                "data" = data,
                "layout" = list(
                  paper_bgcolor = cardColor[[theme]],
                  plot_bgcolor = cardColor[[theme]],
                  automargin = TRUE,
                  font = list(color = textColor[[theme]], size = 12),
                  xaxis = list(
                    "color" = gridColor[[theme]],
                    "gridcolor" = gridColor[[theme]]
                  ),
                  yaxis = list(
                    "color" = gridColor[[theme]],
                    "gridcolor" = gridColor[[theme]]
                  )
                )
              )
            ),
            htmlDiv(
              id = "bottom-card",
              style = list(
                "backgroundColor" = cardColor[[theme]],
                "color" = textColor[[theme]]
              ),
              children = list(# Display the sourced and measured values
                htmlDiv(
                  id = "measure-div",
                  children = list(
                    daqLEDDisplay(
                      id = "source-display",
                      label = sprintf("Applied %s (%s)",
                                      sourceLabel, sourceUnit),
                      value = measSrc,
                      color = accentColor[[theme]]
                    ),
                    daqLEDDisplay(
                      id = "measure-display",
                      label = sprintf("Measured %s (%s)",
                                      measureLabel, measureUnit),
                      value = measDisplay,
                      color = accentColor[[theme]]
                    )
                  )
                )
              )
            )
          )
        ),
        htmlDiv(
          id = "control-container",
          className = "five columns",
          children = list(
            htmlDiv(
              # Controls and options for the IV tracer
              id = "up-control-card",
              style = list(
                "backgroundColor" = cardColor[[theme]],
                "color" = textColor[[theme]]
              ),
              children = list(
                htmlDiv(
                  id = "control-sections",
                  children = list(
                    htmlDiv(
                      className = "IV-source-options",
                      children = list(
                        htmlLabel("Sourcing",
                                  title = paste0("Choose whether you want to ",
                                    "source voltage and measure current, or ",
                                    "source current and measure voltage",
                                    sep = "")
                        ),
                        daqToggleSwitch(
                          id = "source-choice-toggle",
                          label = list("Voltage", "Current"),
                          style = list("width" = "150px",
                                       "margin" = "auto"),
                          value = sourceToggle
                        )
                      )
                    ),
                    htmlDiv(
                      className = "measure-options",
                      children = list(
                        htmlLabel(
                          "Measure mode",
                          title = paste0("Choose if you want to do single ",
                          "measurement or to start a sweep mode",
                          sep = "")
                        ),
                        daqToggleSwitch(
                          id = "mode-choice-toggle",
                          label = list("Single measure", "Sweep"),
                          style = list("width" = "150px"),
                          value = modeToggle
                        )
                      )
                    )
                  )
                ),
                daqStopButton(
                  id = "clear-graph_btn",
                  buttonText = "Clear Graph",
                  className = "daq-button",
                  size = 120
                ),
                daqIndicator(
                  id = "clear-graph_ind",
                  color = accentColor[[theme]],
                  value = FALSE
                )
              )
            ),
            htmlDiv(
              id = "mid-control-card",
              style = list(
                "backgroundColor" = cardColor[[theme]],
                "color" = textColor[[theme]],
                "marginTop" = "10px"
              ),
              children = list(
                # Sourcing controls
                htmlDiv(
                  id = "source-div",
                  children = list(
                    # To perform single measures
                    # adjusting the source with a knob
                    htmlDiv(
                      id = "single_div",
                      className = "single_div_toggle_style",
                      children = list(
                        daqKnob(
                          id = "source-knob",
                          size = 100,
                          value = srcKnob,
                          min = 0,
                          max = 10,
                          color = accentColor[[theme]],
                          label = sprintf("%s (%s)", sourceLabel, sourceUnit)
                        ),
                        daqLEDDisplay(
                          id = "source-knob-display",
                          label = "Knob readout",
                          value = srcKnob,
                          color = accentColor[[theme]]
                        )
                      )
                    ),
                    # To perform automatic sweeps of the source
                    htmlDiv(
                      id = "sweep_div",
                      className = "sweep_div_toggle_style",
                      children = list(
                        htmlDiv(
                          className = "sweep-div-row",
                          children = list(
                            htmlDiv(
                              className = "sweep-div-row",
                              style = list("width" = "98%"),
                              children = list(
                                htmlDiv(
                                  id = "sweep-title",
                                  children = htmlP(
                                    sprintf("%s sweep", sourceLabel)
                                  )
                                ),
                                htmlDiv(
                                  children = list(
                                    daqIndicator(
                                      id = "sweep-status",
                                      label = "Sweep active",
                                      color = accentColor[[theme]],
                                      value = FALSE
                                    )
                                  ),
                                  title = "Indicates if the sweep is running"
                                )
                              )
                            )
                          )
                        ),
                      htmlDiv(
                        className = "sweep-div-row",
                        children = list(
                          htmlDiv(
                            className = "sweep-div-row-inner",
                            children = list(
                              "Start",
                              daqPrecisionInput(
                                id = "sweep-start",
                                precision = 4,
                                label = sprintf(" %s", sourceUnit),
                                labelPosition = "right",
                                value = 0,
                                style = list("marginLeft" = "5px")
                              )
                            ),
                            title = "The lowest value of the sweep"
                          ),
                          htmlDiv(
                            className = "sweep-div-row-inner",
                            children = list(
                              "Stop",
                              daqPrecisionInput(
                                id = "sweep-stop",
                                precision = 4,
                                label = sprintf(" %s", sourceUnit),
                                labelPosition = "right",
                                value = 10
                              )
                            ),
                            title = "The highest value of the sweep"
                          )
                        )
                      ),
                      htmlDiv(
                        className = "sweep-div-row",
                        children = list(
                          htmlDiv(
                            className = "sweep-div-row-inner",
                            children = list(
                              "Step",
                              daqPrecisionInput(
                                id = "sweep-step",
                                precision = 4,
                                label = sprintf(" %s", sourceUnit),
                                labelPosition = "right",
                                min = 0.2,
                                value = 0.2
                              )
                            ),
                            title = "The increment of the sweep",
                          ),
                          htmlDiv(
                            className = "sweep-div-row-inner",
                            children = list(
                              "Time of a step",
                              daqNumericInput(
                                id = "sweep-dt",
                                value = 0.2,
                                min = 0.01,
                                style = list("margin" = "5px")
                              ),
                              "s"
                            ),
                            title = "The time spent on each increment"
                          )
                        )
                      )
                  )
                )
              )
            ),
            # Measure button and indicator
            htmlDiv(
              id = "trigger-div",
              children = list(
                daqStopButton(
                  id = "trigger-measure_btn",
                  buttonText = "Single measure",
                  className = "daq-button",
                  size = 120
                ),
                daqIndicator(
                  id = "measure-triggered",
                  color = accentColor[[theme]],
                  value = FALSE,
                  label = "Measure active"
                )
              )
            )
              )
            )
          )
        )
      )
    )
  )

  if (theme == "dark") {
    return(daqDarkThemeProvider(children = htmlLayout))
  }
  if (theme == "light") {
    return(children = htmlLayout)
  }
}

GenerateModal <- function() {

  return(htmlDiv(
    id = "markdown",
    className = "modal",
    style = list("display" = "none"),
    children = htmlDiv(
      id = "markdown-container",
      className = "markdown-container",
      style = list("color" = textColor[["light"]],
                   "backgroundColor" = cardColor[["light"]]),
      children = list(
        htmlDiv(
          className = "close-container",
          children = htmlButton(
            "Close",
            id = "markdown_close",
            n_clicks = 0,
            className = "closeButton",
            style = list("color" = textColor[["dark"]])
          )
        ),
        htmlDiv(
          className = "markdown-text",
          children = dccMarkdown(
            '
            **What is this mock app about?**

            This is an app to show the graphic elements of Dash DAQ used to create an
            interface for an IV curve tracer using a Keithley 2400 SourceMeter.
            Please note that this application does not interface with a physical instrument;
            the values displayed are simulated using an IV curve model for demonstration purposes.

            **How to use the application**

            First, select either current or voltage as the source using the toggle
            switch below the **Sourcing** label.
            Then choose if you want to operate in a **single measurement mode** or in a **sweep mode**.

            ***Single measurement mode***

            When selecting single measurement mode, a knob will appear to specify
            current (amperes) or voltage (volts). Use this knob to adjust the value
            of the source, then click the "Single Measure" button to display the result.
            Repeating this process using a range of source values will reveal the full IV curve.

            ***Sweep mode***

            Set the sweep parameters `start`, `stop` and `step` as well as the time
            spent on each step, then click on the button `Start Sweep`, the result of the
            sweep will be displayed on the graph.

            The data is never erased unless the button `Clear Graph` is pressed, or if the
            source type is changed.

            ***Dark/light theme***

            Click on theme toggle on top of the page to view dark/light layout.

            You can check out the Dash DAQ components at [
            dashdaq.io](https://www.dashdaq.io/)
            '
          )
        )
      )
    )
  ))
}

app$layout(htmlDiv(
  id = "main-page",
  className = "container",
  style = list("backgroundColor" = bkgColor[["light"]]),
  children = list(
    dccLocation(id = "url", refresh = FALSE),
    dccInterval(id = "refresher", interval = 1000000),
    htmlDiv(
      id = "header",
      className = "banner",
      style = list("backgroundColor" = bannerColor[["light"]],
                   "color" = textColor[["light"]]),
      children = list(
        htmlImg(src = "/assets/dash-logo.png",
                className = "logo three columns"),
        htmlH6("Dash DAQ: IV Curve Tracer", className = "title six columns"),
        htmlDiv(
          className = "three columns",
          style = list("float" = "left"),
          children = list(
            daqToggleSwitch(
              id = "toggleTheme",
              label = list("Light", "Dark"),
              style = list(
                "margin" = "auto",
                "width" = "65%",
                "color" = textColor[["light"]]
              ),
              value = FALSE,
              size = 35,
            )
          )
        )
      )
    ),
    htmlDiv(
      id = "intro-banner",
      className = "intro-banner",
      style = list("color" = "#FFFFFF",
                   "backgroundColor" = accentColor[["light"]]),
      children = htmlDiv(
        className = "intro-banner-content",
        children = list(
          htmlP(
            children = "This app uses graphic elements of Dash DAQ to create
            an interface for an IV curve tracer using a Keithley 2400 SourceMeter.
            Please note that this application does not interface with a physical instrument;
            the values displayed are simulated using an IV curve model for demonstration purposes.",
            className = "intro-banner-text"
          ),
          htmlButton(
            id = "learn-more-button",
            children = "Learn More",
            n_clicks = 0,
            style = list(
              "borderColor" = bkgColor[["light"]],
              "color" = "#FFFFFF",
              "backgroundColor" = accentColor[["light"]]
            )
          )
        )
      )
    ),
    htmlDiv(
      id = "page-content",
      children = GenerateMainLayout(),
      className = "flex-display",
      style = list("backgroundColor" = bkgColor[["light"]], "padding" = "2%")
    ),
    GenerateModal(),
    dccStore(id = "store-sourced-values", data = list(NA)),
    dccStore(id = "store-measured-values", data = list(NA)),
    dccStore(id = "store-source", data = c("V", FALSE)),
    dccStore(id = "store-nclick", data = 0),
    dccStore(id = "store-mode-choice", data = FALSE)
  )
))


# ======= Dark/light themes callbacks =======
app$callback(
  output = list(id = "page-content", property = "children"),
  params = list(
    input(id = "toggleTheme", property = "value"),
    state(id = "source-choice-toggle", property = "value"),
    state(id = "mode-choice-toggle", property = "value"),
    state(id = "IV_graph", property = "figure"),
    state(id = "source-display", property = "value"),
    state(id = "measure-display", property = "value"),
    state(id = "source-knob", property = "value"),
    state(id = "source-choice-toggle", property = "value"),
    state(id = "mode-choice-toggle", property = "value")
  ),

  # ''' Update the theme of the daq components '''
  function(value, srcChoice, modeChoice, fig,
           measSrc, measDisplay, srcKnob,
           sourceToggle, modeToggle) {

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    if (modeChoice) {
      modeVal <- "Sweep"
    } else {
      modeVal <- "Single Measure"
    }

    if (value) {
      return(GenerateMainLayout("dark", srcType, modeVal, fig,
                                measSrc, measDisplay, srcKnob,
                                sourceToggle, modeToggle))
    } else {
      return(GenerateMainLayout("light", srcType, modeVal, fig,
                                measSrc, measDisplay, srcKnob,
                                sourceToggle, modeToggle))
    }
  }
)


app$callback(
  output = list(id = "page-content", property = "style"),
  params = list(
    input(id = "toggleTheme", "value"),
    state(id = "page-content", property = "style")
  ),

  # ''' Update the theme of the app '''
  function(value, styleList) {

    if (value) {
      theme <- "dark"
    } else {
      theme <- "light"
    }

    styleList[["color"]] <- textColor[[theme]]
    styleList[["backgroundColor"]] <- bkgColor[[theme]]

    return(styleList)
  }
)


app$callback(
  output = list(id = "header", property = "style"),
  params = list(
    input(id = "toggleTheme", "value"),
    state(id = "header", property = "style")
  ),

  # ''' Update the theme of the header '''
  function(value, styleList) {

    if (value) {
      theme <- "dark"
    } else {
      theme <- "light"
    }

    styleList[["color"]] <- textColor[[theme]]
    styleList[["backgroundColor"]] <- bkgColor[[theme]]

    return(styleList)
  }
)


app$callback(
  output = list(id = "intro-banner", property = "style"),
  params = list(
    input(id = "toggleTheme", "value"),
    state(id = "intro-banner", property = "style")
  ),

  # ''' Update the theme of the banner '''
  function(value, styleList) {

    if (value) {
      theme <- "dark"
    } else {
      theme <- "light"
    }

    styleList[["color"]] <- "#FFFFFF"

    return(styleList)
  }
)


app$callback(
  output = list(id = "markdown-container", property = "style"),
  params = list(
    input(id = "toggleTheme", "value"),
    state(id = "markdown-container", property = "style")
  ),

  # ''' Update the theme of markdown '''
  function(value, styleList) {

    if (value) {
      theme <- "dark"
    } else {
      theme <- "light"
    }

    styleList[["color"]] <- textColor[[theme]]
    styleList[["backgroundColor"]] <- cardColor[[theme]]

    return(styleList)
  }
)


app$callback(
  output = list(id = "main-page", property = "style"),
  params = list(
    input(id = "toggleTheme", "value"),
    state(id = "main-page", property = "style")
  ),

  # ''' Update the theme of entire page '''
  function(value, styleList) {

    if (value) {
      theme <- "dark"
    } else {
      theme <- "light"
    }

    styleList[["color"]] <- textColor[[theme]]
    styleList[["backgroundColor"]] <- bkgColor[[theme]]

    return(styleList)

  }
)


# ======= Callback for modal popup =======
app$callback(
  output = list(id = "markdown", property = "style"),
  params = list(
    input(id = "learn-more-button", "n_clicks"),
    input(id = "markdown_close", property = "n_clicks")
  ),

  function(buttonClick, closeClick) {

    if (buttonClick > closeClick) {
      return(list("display" = "block"))
    } else {
      return(list("display" = "none"))
    }

    styleList[["color"]] <- textColor[[theme]]
    styleList[["backgroundColor"]] <- bkgColor[[theme]]

    return(styleList)
  }
)


# ======= Callbacks for changing labels =======
# ======= Label for single measures, sweep mode, displays =======
app$callback(
  output = list(id = "source-knob", property = "label"),
  params = list(
    input(id = "source-choice-toggle", "value")
  ),

  function(srcChoice) {

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    sourceLabel <- GetSourceLabels(srcType)[1]

    return(sourceLabel)
  }
)

app$callback(
  output = list(id = "source-knob-display", property = "label"),
  params = list(
    input(id = "source-choice-toggle", "value")
  ),

  function(srcChoice) {

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    sourceLabel <- GetSourceLabels(srcType)[1]
    sourceUnit <- GetSourceUnits(srcType)[1]

    return(sprintf("Value : %s (%s)", sourceLabel, sourceUnit))
  }
)


app$callback(
  output = list(id = "sweep-start", property = "label"),
  params = list(
    input(id = "source-choice-toggle", "value")
  ),

  function(srcChoice) {

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    sourceUnit <- GetSourceUnits(srcType)[1]

    return(sprintf("(%s)", sourceUnit))
  }
)


app$callback(
  output = list(id = "sweep-stop", property = "label"),
  params = list(
    input(id = "source-choice-toggle", "value")
  ),

  function(srcChoice) {

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    sourceUnit <- GetSourceUnits(srcType)[1]

    return(sprintf("(%s)", sourceUnit))
  }
)


app$callback(
  output = list(id = "sweep-step", property = "label"),
  params = list(
    input(id = "source-choice-toggle", "value")
  ),

  function(srcChoice) {

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    sourceUnit <- GetSourceUnits(srcType)[1]

    return(sprintf("(%s)", sourceUnit))
  }
)


app$callback(
  output = list(id = "sweep-title", property = "children"),
  params = list(
    input(id = "source-choice-toggle", "value")
  ),

  function(srcChoice) {

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    sourceLabel <- GetSourceLabels(srcType)[1]

    return(htmlP(sprintf("%s sweep", sourceLabel)))
  }
)


app$callback(
  output = list(id = "source-display", property = "label"),
  params = list(
    input(id = "source-choice-toggle", "value")
  ),

  function(srcChoice) {

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    sourceLabel <- GetSourceLabels(srcType)[1]
    sourceUnit <- GetSourceUnits(srcType)[1]

    return(sprintf("Applied %s (%s)", sourceLabel, sourceUnit))

  }
)


app$callback(
  output = list(id = "measure-display", property = "label"),
  params = list(
    input(id = "source-choice-toggle", "value")
  ),

  function(srcChoice) {

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    measureLabel <- GetSourceLabels(srcType)[2]
    measureUnit <- GetSourceUnits(srcType)[2]

    return(sprintf("Measured %s (%s)", measureLabel, measureUnit))

  }
)


app$callback(
  output = list(id = "trigger-measure_btn", property = "buttonText"),
  params = list(
    input(id = "mode-choice-toggle", property = "value")
  ),

  # ''' Update the measure button upon choosing single or sweep '''
  function(modeChoice) {

    if (modeChoice) {
      return("Start sweep")
    } else {
      return("Single measure")
    }
  }
)


# ======= Callbacks to change elements in the layout =======
app$callback(
  output = list(id = "single_div", property = "style"),
  params = list(
    input(id = "mode-choice-toggle", property = "value")
  ),

  # '''toggle the layout for single measure'''
  function(modeChoice) {

    if (modeChoice) {
      return(list("display" = "none"))
    } else {
      return(list("display" = "flex",
                  "flexDirection" = "column",
                  "alignItems" = "center",
                  "justifyContent" = "space-around"))
    }
  }
)


app$callback(
  output = list(id = "sweep_div", property = "style"),
  params = list(
    input(id = "mode-choice-toggle", property = "value")
  ),

  # '''toggle the layout for sweep'''
  function(modeChoice) {

    if (modeChoice) {
      return(list("display" = "flex",
                  "flexDirection" = "column",
                  "alignItems" = "center",
                  "justifyContent" = "space-around"))

    } else {
      return(list("display" = "none"))
    }
  }
)


# ======= Applied/measured values display =======
app$callback(
  output = list(id = "store-source", property = "data"),
  params = list(
    input(id = "source-choice-toggle", property = "value"),
    state(id = "store-source", property = "data"),
    state(id = "source-knob", property = "value")
  ),

  # '''
  # This callback updates the store-source's
  # source value ("V"/"I") and is_sourced_changed T/F
  # status similar to python version
  # '''
  function(srcType, hiddenSource, knobVal) {

    if (srcType) {
      srcT <- "I"
    } else {
      srcT <- "V"
    }

    if (srcT != hiddenSource[[1]]) {
      hiddenSource[[1]] <- srcT
      hiddenSource[[2]] <- TRUE
      return(hiddenSource)
    } else {
      hiddenSource[[2]] <- FALSE
      return(hiddenSource)
    }
  }
)


app$callback(
  output = list(id = "source-knob", property = "value"),
  params = list(
    input(id = "source-choice-toggle", property = "value"),
    state(id = "store-source", property = "data"),
    state(id = "source-knob", property = "value")
  ),

  # '''
  # Reset the knob value when
  # source changed
  # '''
  function(srcType, hiddenSource, knobVal) {
    if (srcType) {
      srcT <- "I"
    } else {
      srcT <- "V"
    }
    if (srcT != hiddenSource[[1]]) {
      return(0)
    } else {
      return(knobVal)
    }
  }
)


# ======= Interval callbacks =======
app$callback(
  output = list(id = "refresher", property = "interval"),
  params = list(
    input(id = "sweep-status", property = "value"),
    state(id = "mode-choice-toggle", property = "value"),
    state(id = "sweep-dt", property = "value")
  ),

  # ''' Change the interval to high frequency for sweep '''
  function(swpOn, modeChoice, sweepDt) {


    if (modeChoice && swpOn) {
      return (sweepDt * 1000)
    } else {
      return(1000000)
    }
  }
)


app$callback(
  output = list(id = "refresher", property = "n_intervals"),
  params = list(
    input(id = "trigger-measure_btn", property = "n_clicks"),
    input(id = "mode-choice-toggle", property = "value"),
    state(id = "sweep-status", property = "value"),
    state(id = "refresher", property = "n_intervals")

  ),

  # ''' Reset the n_interval of the dccInterval once a sweep is done'''
  function(clicks, modeChoice, swpOn, nInterval) {


    if (modeChoice && swpOn) {
      return(nInterval)
    } else {
      return(0)
    }
  }
)


app$callback(
  output = list(id = "sweep-status", property = "value"), #daqIndicator
  params = list(
    input(id = "trigger-measure_btn", property = "n_clicks"),
    input(id = "source-display", property = "value"),
    state(id = "measure-triggered", property = "value"),
    state(id = "sweep-status", property = "value"),
    state(id = "sweep-stop", property = "value"),
    state(id = "sweep-step", property = "value"),
    state(id = "mode-choice-toggle", property = "value")
  ),

  # '''
  # Decide whether to turn on or off the sweep
  # when single mode is selected, it is off by default
  # when sweep mode is selected, it enables the sweep if it wasn't on
  # otherwise it stops the sweep once the sourced value gets higher or equal
  # than the sweep limit minus the sweep step
  #'''
  function(clicks, sourcedVal, measTriggered,
           swpOn, swpStop, swpStep, modeChoice) {

    if (!(modeChoice)) {
      return(FALSE)
    } else {
      if (swpOn) {
        # The condition of continuation is to source lower than the sweep
        # limit minus one sweep step

        answer <- sourcedVal <= swpStop - swpStep

        return(answer)
      } else {
          if (!(measTriggered)) {
            return(FALSE)
          }
          # Initiate a sweep
          return(TRUE)
      }
    }
  }
)


# ======= Measurements callbacks =======
app$callback(
  output = list(id = "source-knob-display", property = "value"),
  params = list(
    input(id = "source-knob", property = "value")
  ),

  # ''' Set the value of the knob on a LED display '''
  function(knobVal) {
    return(knobVal)
  }
)


app$callback(
  output = list(id = "measure-triggered", property = "value"),
  params = list(
    input(id = "trigger-measure_btn", property = "n_clicks"),
    input(id = "mode-choice-toggle", property = "value"),
    state(id = "store-nclick", property = "data"),
    state(id = "store-mode-choice", property = "data")
  ),

  # '''
  # Controls if a measure can be made or not
  # the indicator 'measure-triggered' can be set to TRUE only by a click
  # on the 'trigger-measure_btn' button or by the 'refresher' interval
  # '''
  function(nClick, trigger, storeNclick, storeMode) {

    if (is.null(unlist(nClick))) {
      nClick <- 0
    }

    if (is.null(unlist(storeNclick))) {
      storeNclick <- 0
    }

    if (nClick != storeNclick && trigger == storeMode) {
      # It was triggered by a click on the trigger-measure_btn button
      return(TRUE)
    } else {
      # It was triggered by a change of the mode
      return(FALSE)
    }
  }
)


app$callback(
  output = list(id = "store-mode-choice", property = "data"),
  params = list(
    input(id = "mode-choice-toggle", property = "value"),
    state(id = "store-mode-choice", property = "data")
  ),

  # ''' Store & Update modeChoice '''
  function(currentMode, prevMode) {

    if (currentMode == prevMode) {
      answer <- prevMode
    } else {
      answer <- currentMode
    }
    return(answer)
  }
)


app$callback(
  output = list(id = "source-display", property = "value"),
  params = list(
    input(id = "refresher", property = "n_intervals"),
    input(id = "measure-triggered", property = "value"),
    state(id = "source-knob", property = "value"),
    state(id = "source-display", property = "value"),
    state(id = "sweep-start", property = "value"),
    state(id = "sweep-stop", property = "value"),
    state(id = "sweep-step", property = "value"),
    state(id = "mode-choice-toggle", property = "value"),
    state(id = "sweep-status", property = "value")
  ),


  # ''' Set the source value to the instrument '''
  function(nInterval, measTriggered, knobVal,
           oldSourceDisplayVal, swpStart,
           swpStop, swpStep, modeChoice, swpOn) {

    # Default answer
    answer <- oldSourceDisplayVal

    if (!(modeChoice)) {
      answer <- knobVal
    } else {
      if (measTriggered) {
        if (swpOn) {
          answer <- swpStart + (nInterval - 1) * swpStep
          if (answer > swpStop) {
            answer <- oldSourceDisplayVal
          }
        }
      }
    }
    return(answer)
  }
)


app$callback(
  output = list(id = "measure-display", property = "value"),
  params = list(
    input(id = "source-display", property = "value"),
    state(id = "measure-triggered", property = "value"),
    state(id = "measure-display", property = "value"),
    state(id = "source-choice-toggle", property = "value"),
    state(id = "mode-choice-toggle", property = "value"),
    state(id = "sweep-status", property = "value")
  ),

  # '''
  # Read the measured value from the instrument
  # check if a measure should be made
  # initiate a measure of the KT2400
  # read the measure value and return it
  # by default it simply returns the value previously available
  # '''
  function(srcVal, measTriggered, measOldVal, srcChoice, modeChoice, swpOn) {

    if (srcChoice) {
      srcType <- "I"
    } else {
      srcType <- "V"
    }

    measuredValue <- measOldVal

    if (!(modeChoice)) { # Single measure
      if (measTriggered) {
        # Initiate a measurement
        measuredValue <- SourceAndMeasure(srcType, srcVal)
      }
    } else { # Sweep
        if (measTriggered && swpOn) {
          # Initiate a measurement
          measuredValue <- SourceAndMeasure(srcType, srcVal)
        }
    }
    return(round(measuredValue, 2))
  }
)


app$callback(
  output = list(id = "store-sourced-values", property = "data"),
  params = list(
    input(id = "source-display", property = "value"),
    input(id = "clear-graph_btn", property = "n_clicks"),
    state(id = "measure-triggered", property = "value"),
    state(id = "source-choice-toggle", property = "value"),
    state(id = "mode-choice-toggle", property = "value"),
    state(id = "sweep-status", property = "value"),
    state(id = "store-sourced-values", property = "data"),
    state(id = "store-nclick", property = "data")

  ),


  # ''' Save sourcedValue to store-sourced-value to retrieve for graph later '''
  function(srcVal, nClicks, measTriggered, srcChoice,
           modeChoice, swpOn, sourceStore, storeNclick) {

    if (srcChoice) {
      srcType <- "I"
    } else {
      srcType <- "V"
    }

    if (is.list(nClicks)){
      nClicks <- 0
    }
    if (is.list(storeNclick)){
      storeNclick <- 0
    }

    if (nClicks == storeNclick) { # clear button not clicked
      if (!(modeChoice)) { # single measure
        if (measTriggered) {
          # Save the sourced value
          if (is.null(sourceStore[[1]])) {
            sourceStore <- list(srcVal)
          } else {
            sourceStore <- c(sourceStore, srcVal)
          }
        }
      } else { # sweep mode
        if (measTriggered && swpOn) {
          # Save the sourced value
          if (is.null(sourceStore[[1]])) {
            sourceStore <- list(srcVal)
          } else {
            sourceStore <- c(sourceStore, srcVal)
          }
        }
      }
    } else { # clear button clicked
        sourceStore <- NA
    }
    return(sourceStore)
  }
)


app$callback(
  output = list(id = "store-measured-values", property = "data"),
  params = list(
    input(id = "source-display", property = "value"),
    input(id = "clear-graph_btn", property = "n_clicks"),
    state(id = "measure-triggered", property = "value"),
    state(id = "measure-display", property = "value"),
    state(id = "source-choice-toggle", property = "value"),
    state(id = "mode-choice-toggle", property = "value"),
    state(id = "sweep-status", property = "value"),
    state(id = "store-measured-values", property = "data"),
    state(id = "store-nclick", property = "data")
  ),

  # ''' Save measuredValue to dccStore to retrieve for graph later '''
  function(srcVal, nClicks, measTriggered, measuredValue,
           srcChoice, modeChoice, swpOn, measStore, storeNclick) {

    if (srcChoice) {
      srcType <- "I"
    } else {
      srcType <- "V"
    }

    if (is.list(nClicks)) {
      nClicks <- 0
    }
    if (is.list(storeNclick)) {
      storeNclick <- 0
    }
    if (nClicks == storeNclick) { # clear button not clicked
      if (!(modeChoice)) {
        if (measTriggered) {
          # Initiate a measurement
          measuredValue <- SourceAndMeasure(srcType = srcType, srcVal = srcVal)
          # Save the measured value
          if (is.null(measStore[[1]])) {
            measStore <- list(measuredValue)
          } else {
            measStore <- c(measStore, measuredValue)
          }
        }
      } else {
        if (measTriggered && swpOn) {
          # Save the measured value
          measuredValue <- SourceAndMeasure(srcType = srcType, srcVal = srcVal)
          if (is.null(measStore[[1]])) {
            measStore <- list(measuredValue)
          } else {
            measStore <- c(measStore, measuredValue)
          }
        }
      }
    } else { # clear button clicked
        measStore <- NA
    }
    return(measStore)
  }
)


# ======= Graph related callbacks =======
app$callback(
  output = list(id = "store-nclick", property = "data"),
  params = list(
    input(id = "clear-graph_btn", property = "n_clicks"),
    state(id = "store-nclick", property = "data")
  ),

  # ''' Update storedNclick value '''
  function(nClick, storeNclick) {

    if (is.list(nClick)) {
      nClick <- 0
    }
    if (nClick != storeNclick) {
      storeNclick <- nClick
    }
    return(storeNclick)
  }
)


app$callback(
  output = list(id = "clear-graph_ind", property = "value"),
  params = list(
    input(id = "clear-graph_btn", property = "n_clicks"),
    input(id = "source-knob", property = "value"),
    input(id = "measure-triggered", property = "value"),
    state(id = "store-nclick", property = "data")
  ),

  # '''
  # Turn on indicator light, if clear-graph_btn clicked
  # Turn off indicator light, if source-knob or measure-triggered values change
  # '''
  function(nClick, knobVal, measTriggered, storeNclick) {

    if (is.list(nClick)) {
      nClick <- 0
    }
    if (is.list(storeNclick)) {
      storeNclick <- 0
    }

    if (nClick != storeNclick) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
)


app$callback(
  output = list(id = "IV_graph", property = "figure"),
  params = list(
    input(id = "measure-display", property = "value"),
    input(id = "clear-graph_ind", property = "value"),
    state(id = "toggleTheme", property = "value"),
    state(id = "measure-triggered", property = "value"),
    state(id = "IV_graph", property = "figure"),
    state(id = "source-choice-toggle", property = "value"),
    state(id = "mode-choice-toggle", property = "value"),
    state(id = "sweep-status", property = "value"),
    state(id = "store-sourced-values", property = "data"),
    state(id = "store-measured-values", property = "data")
  ),

  # ''' Update the IV-Graph '''
  function(measuredVal,
           clearGraph, # to fire callback on reset
           theme,
           measTriggered,
           graphData,
           srcChoice,
           modeChoice,
           swpOn,
           sourceStore,
           measStore) {

    if (theme) {
      theme <- "dark"
    } else {
      theme <- "light"
    }

    if (srcChoice) {
      srcType <- "Current"
    } else {
      srcType <- "Voltage"
    }

    # Labels for sourced and measured quantities
    sourceLabel <- GetSourceLabels(srcType)[1]
    measureLabel <- GetSourceLabels(srcType)[2]
    sourceUnit <- GetSourceUnits(srcType)[1]
    measureUnit <- GetSourceUnits(srcType)[2]

    # Sorting data
    if (!is.null(sourceStore[[1]])) {
      xdata <- sort(unlist(sourceStore), decreasing = TRUE)
      ydata <- sort(unlist(measStore), decreasing = FALSE)
    } else {
      xdata <- NA
      ydata <- NA
    }

    figureData <- list(
      list(
        x = xdata,
        y = ydata,
        mode = "lines+markers",
        name = "IV curve",
        line = list("color" = accentColor[[theme]], "width" = 2)
      )
    )
    ivFigure <- list(
      "data" = figureData,
      "layout" = list(
        xaxis = list(
          "title" = sprintf("Applied %s (%s)", sourceLabel, sourceUnit),
          "color" = textColor[[theme]],
          "gridcolor" = gridColor[[theme]]
        ),
        yaxis = list(
          "title" = sprintf("Measured %s (%s)", measureLabel, measureUnit),
          "color" = textColor[[theme]],
          "gridcolor" = gridColor[[theme]]
        ),
        font = list(color = textColor[[theme]], size = 12),
        automargin = TRUE,
        plot_bgcolor = cardColor[[theme]],
        paper_bgcolor = cardColor[[theme]]
      )
    )

    if (!(modeChoice)) { # Single measure case
      if (measTriggered) {

        # Single data point does not draw graph,
        # duplicating data to show single data point
        if (length(xdata) == 1 && !(NA %in% xdata)) {
          xdata <- rep(xdata, 2)
          ydata <- rep(ydata, 2)
          figureData[[1]][["x"]] <- xdata
          figureData[[1]][["y"]] <- ydata
          ivFigure[["data"]] <- figureData
        }
      }
    }
    return(ivFigure)
  }
)

if (!appName == ""){
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(debug = TRUE)
}
