library(dash)
library(dashDaq)
library(dashCoreComponents)
library(dashHtmlComponents)

setwd("/Users/Caner/Desktop/plotly/dashR-daq-iv-tracer/")

# List to store information useful to callbacks
localVars <- list()

localVars['nClicks'] <- 0
localVars['nClicksClearGraph'] <- 0

# Define the app
app <- Dash$new()

# '''labels for source/measure elements'''
GetSourceLabels <- function(source = "Voltage") {

  if (source == "Voltage") {
    # we source voltage and measure current
    sourceLabel <- "Voltage"
    measureLabel <- "Current"
  } else if (source == "Current") {
    # we source current and measure voltage
    sourceLabel <- "Current"
    measureLabel <- "Voltage"
  }
  return(c(sourceLabel, measureLabel))
}

# '''units for source/measure elements'''
GetSourceUnits <- function(source = "Voltage") {

  if (source == "Voltage") {
    # we source voltage and measure current
    sourceUnit <- "V"
    measureUnit <- "A"
  } else if (source == "Current") {
    # we source current and measure voltage
    sourceUnit <- "A"
    measureUnit <- "V"
  }
  return(c(sourceUnit, measureUnit))
}

# Font and background colors associated with each theme
bannerColor = list("dark" = "#23262e", "light" = "#ffffff")
bkgColor = list("dark" = "#23262e", "light" = "#f6f6f7")
gridColor = list("dark" = "#53555B", "light" = "#969696")
textColor = list("dark" = "#95969A", "light" = "#595959")
cardColor = list("dark" = "#2D3038", "light" = "#FFFFFF")
accentColor = list("dark" = "#FFD15F", "light" = "#ff9827")

# """generate the layout of the app"""
GenerateMainLayout <- function(theme = "light",
                               srcType = "Voltage",
                               modeVal ="Single measure",
                               fig = NULL) {

  sourceLabel <- GetSourceLabels(srcType)[1]
  measureLabel <- GetSourceLabels(srcType)[2]
  sourceUnit <- GetSourceUnits(srcType)[1]
  measureUnit <- GetSourceUnits(srcType)[2]

  # As the trigger-measure btn will have its n_clicks reset by the reloading
  # of the layout we need to reset this one as well
  localVars['nClicks'] <- 0
  localVars['nClicksClearGraph'] <- 0

  # Doesn't clear the data of the graph
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
                "color" = textColor[[theme]],
                "marginTop" = "10px"
              ),
              children = list(# Display the sourced and measured values
                htmlDiv(
                  id = "measure-div",
                  children = list(
                    daqLEDDisplay(
                      id = "source-display",
                      label = sprintf("Applied %s (%s)", sourceLabel, sourceUnit),
                      value = 0.00,
                      color = accentColor[[theme]]
                    ),
                    daqLEDDisplay(
                      id = "measure-display",
                      label = sprintf("Measured %s (%s)", measureLabel, measureUnit),
                      value = 0.0000,
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
              # controls and options for the IV tracer
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
                                  title = "Choose whether you want to source voltage and measure current, or source current and measure voltage"
                        ),
                        daqToggleSwitch(
                          id = "source-choice-toggle",
                          label = list("Voltage", "Current"),
                          style = list("width" = "150px",
                                       "margin" = "auto"),
                          value = FALSE
                        )
                      )
                    ),
                    htmlDiv(
                      className = "measure-options",
                      children = list(
                        htmlLabel(
                          "Measure mode",
                          title = "Choose if you want to do single measurement or to start a sweep mode"
                        ),
                        daqToggleSwitch(
                          id = "mode-choice-toggle",
                          label = list("Single measure", "Sweep"),
                          style = list("width" = "150px"),
                          value = FALSE
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
                    # To perform single measures adjusting the source with a knob
                    htmlDiv(
                      id = "single_div",
                      className = "single_div_toggle_style",
                      children = list(
                        daqKnob(
                          id = "source-knob",
                          size = 100,
                          value = 0.00,
                          min = 0,
                          max = 10,
                          color = accentColor[[theme]],
                          label = sprintf("%s (%s)", sourceLabel, sourceUnit)
                        ),
                        daqLEDDisplay(
                          id = "source-knob-display",
                          label = "Knob readout",
                          value = 0.00,
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
            interface for an IV curve tracer using a Keithley 2400 SourceMeter. This mock
            demo does not actually connect to a physical instrument the values displayed
            are generated from an IV curve model for demonstration purposes.

            **How to use the mock app**

            First choose if you want to source (apply) current or voltage, using the toggle switch under **Sourcing** label.
            Then choose if you want to operate in a **single measurement mode** or in a **sweep mode**.

            ***Single measurement mode***

            Adjust the value of the source with the knob at the bottom of the graph area
            and click on the `SINGLE MEASURE` button, the measured value will be displayed.
            Repetition of this procedure for different source values will reveal the full
            IV curve.

            ***Sweep mode***

            Set the sweep parameters `start`, `stop` and `step` as well as the time
            spent on each step, then click on the button `START SWEEP`, the result of the
            sweep will be displayed on the graph.

            The data is never erased unless the button `CLEAR GRAPH` is pressed, or if the
            source type is changed.

            ***Dark/light theme***

            Click on theme toggle on top of the page to view dark/light layout.

            You can purchase the Dash DAQ components at [
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
  style = list("backgroundColor" = bkgColor[["light"]], "height" = "100vh"),
  children = list(
    dccLocation(id = "url", refresh = FALSE),
    dccInterval(id = "refresher", interval = 1000000),
    htmlDiv(
      id = "header",
      className = "banner",
      style = list("backgroundColor" = bannerColor[["light"]],
                   "color" = textColor[["light"]]),
      children = list(
        htmlImg(src = "/assets/dash-daq-logo.png",
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
      style = list("color" = "#FFFFFF", "backgroundColor" = accentColor[["light"]]),
      children = htmlDiv(
        className = "intro-banner-content",
        children = list(
          htmlP(
            children = "This app uses graphic elements of Dash DAQ to create
            an interface for an IV curve tracer using a Keithley 2400 SourceMeter.
            The mock demo does not actually connect to a physical instrument,
            the values displayed are generated from an IV curve model for
            demonstration purposes.",
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
    GenerateModal()
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
    state(id = "source-display", property = "value"), # Keep measure LED display while changing themes
    state(id = "measure-display", property = "value")
  ),

  # '''update the theme of the daq components'''
  function(value, srcChoice, modeChoice, fig, measSrc, measDisplay) { # CHECK LATER IF last 2 args necessary !

    # print(c("value = ", value))
    # print(c("srcChoice = ", srcChoice))
    # print(c("modeChoice = ", modeChoice))
    # print(c("fig = ", fig))
    # print(c("measSrc = ", measSrc))
    # print(c("measDisplay = ", measDisplay))

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
      return(GenerateMainLayout("dark", srcType, modeVal, fig))
    } else {
      return(GenerateMainLayout("light", srcType, modeVal, fig))
    }

  }
)


app$callback(
  output = list(id = "page-content", property = "style"),
  params = list(
    input(id = "toggleTheme", "value"),
    state(id = "page-content", property = "style")
  ),

  # '''update the theme of the app'''
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

  # '''update the theme of the header'''
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

  # '''update the theme of the banner'''
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

  # '''update the theme of markdown'''
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

  # '''update the theme of entire page'''
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
    input(id = "source-choice-toggle", "value"),
    input(id = "mode-choice-toggle", property = "value") # Only required to fire callback CHECK IF NECESSARY!
  ),

  function(srcChoice, modeChoice) {

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
    input(id = "source-choice-toggle", "value"),
    input(id = "mode-choice-toggle", property = "value")
  ),

  function(srcChoice, modeChoice) {

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
    input(id = "source-choice-toggle", "value"),
    input(id = "mode-choice-toggle", property = "value")
  ),

  function(srcChoice, modeChoice) {

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
    input(id = "source-choice-toggle", "value"),
    input(id = "mode-choice-toggle", property = "value")
  ),

  function(srcChoice, modeChoice) {

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
    input(id = "source-choice-toggle", "value"),
    input(id = "mode-choice-toggle", property = "value")
  ),

  function(srcChoice, modeChoice) {

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
    input(id = "source-choice-toggle", "value"),
    input(id = "mode-choice-toggle", property = "value")
  ),

  function(srcChoice, modeChoice) {

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
    input(id = "source-choice-toggle", "value"),
    input(id = "mode-choice-toggle", property = "value")
  ),

  function(srcChoice, modeChoice) {

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
    input(id = "source-choice-toggle", "value"),
    input(id = "mode-choice-toggle", property = "value")
  ),

  function(srcChoice, modeChoice) {

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

  # '''update the measure button upon choosing single or sweep'''
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
  output = list(id = "source-knob", property = "value"),
  params = list(
    input(id = "source-choice-toggle", property = "value")
  ),

  # '''
  # modification upon source-change change
  # the source type in local_vars reset the
  # knob to zero reset the measured values on the graph
  #'''
  function(srcType) {

    if (srcType) {
      localVars['source'] <- "A" # -> Updating localVars like in the python version
                                 # -> also code `local_vars.is_source_being_changed` later IF NECESSARY
    } else {
      localVars['source'] <- "V"
    }
    return(0) # -> Simply setting the source-knob value to 0 if fireback is called. This handles knob
              # -> RESET GRAPH AFTER CODING THE CALLBACKS FOR IT.
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

  # '''change the interval to high frequency for sweep'''
  function(swpOn, modeChoice, sweepDt) {
    # print(swpOn)
    # print(modeChoice)
    # print(sweepDt)

    if (modeChoice && swpOn) {
      return (sweepDt * 1000) # -> ! TEST THIS AFTER connecting sweep with button
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

  # '''reset the n_interval of the dccInterval once a sweep is done'''
  function(clicks, modeChoice, swpOn, nInterval) {
    # print(swpOn)
    # print(modeChoice)
    # print(sweepDt)

    if (modeChoice && swpOn) {
      return(nInterval)
    } else {
      localVars['nRefresh'] <- 0 # For Parity w python. Delete later if redundant!
      return(0)
    }
  }
)


app$callback(
  output = list(id = "sweep-status", property = "value"),
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
  # decide whether to turn on or off the sweep
  # when single mode is selected, it is off by default
  # when sweep mode is selected, it enables the sweep if is wasn't on
  # otherwise it stops the sweep once the sourced value gets higher or equal
  # than the sweep limit minus the sweep step
  #'''
  function(clicks, sourcedVal, measTriggered, swpOn, swpStop, swpStep, modeChoice) {

    # print(class(sourcedVal))
    # print(swpStop)
    # print(swpStep)

    if (!(modeChoice)) {
      return(FALSE)
    } else {
      if (swpOn) {
        # The condition of continuation is to source lower than the sweep
        # limit minus one sweep step

        #print("in answer")
        answer <- sourcedVal <= swpStop - swpStep
        return(answer)
      } else {

        #print("inside else")
        if (!(measTriggered)) { # -> measTriggered NOT implemented YET CHECK BACK after coding
          # The 'trigger-measure_btn' wasn't pressed yet
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

  # '''set the value of the knob on a LED display'''
  function(knobVal) {
    return(knobVal)

  }
)


app$callback(
  output = list(id = "measure-triggered", property = "value"),
  params = list(
    input(id = "trigger-measure_btn", property = "n_clicks"),
    input(id = "mode-choice-toggle", property = "value")
  ),

  # '''
  # controls if a measure can be made or not
  # the indicator 'measure-triggered' can be set to TRUE only by a click
  # on the 'trigger-measure_btn' button or by the 'refresher' interval
  #'''
  function(nClick, trigger) {
    #print(is.null(unlist(nClick)))
    #print(class(nClick))

    if (is.null(unlist(nClick))) {
      nClick <- 0
    }
    #print(localVars[['nClicks']])

    if (nClick != localVars[['nClicks']]) {
      # It was triggered by a click on the trigger-measure_btn button
      localVars[['nClicks']] <- nClick
      #print(nClick)
      return(TRUE)
    } else {
      # It was triggered by a change of the mode
      return(FALSE)
    }
  }
)


app$callback(
  output = list(id = "sweep-status", property = "value"),
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


  # '''set the source value to the instrument'''
  function(nInterval, measTriggered, knobVal,
           oldSourceDisplayVal, swpStart,
           swpStop, swpStep, modeChoice, swpOn) {
    print(oldSourceDisplayVal)
  }
)




app$run_server(port = 8896, debug = TRUE)




