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

# colors for layout and figure
colorsRaw <- read.table("colors.txt", sep = ' ', comment.char = '',
  stringsAsFactors = FALSE)

colors <- setNames(as.list(colorsRaw[, 2]), colorsRaw[, 1])

app <- Dash$new()


############################
# Controls
############################

# integration time, microseconds
intTime <- htmlDiv(children = list(
  htmlDiv(children = list('int. time (\U{00B5}s)'), className = 'option-name'),
  daqNumericInput(
    id = 'integration-time-input',
    value = 1000,
    size = 150,
    min = 1000,
    max = 650000000,
    disabled = FALSE
  )
),
id = 'integration-time')

# scans to average over
nscansAvg <- htmlDiv(children = list(
  htmlDiv(children = list('number of scans'), className = 'option-name'),
  daqNumericInput(
    id = 'nscans-to-average-input',
    value = 1,
    size = 150,
    min = 1,
    max = 100,
    disabled = FALSE
  )
),
id = 'nscans-to-average')

# strobe
strobeEnable <- htmlDiv(children = list(
  htmlDiv(children = list('strobe'), className = 'option-name'),
  daqBooleanSwitch(
    id = 'continuous-strobe-toggle-input',
    color = colors[['accent']],
    on = FALSE,
    disabled = FALSE
  )
),
id = 'continuous-strobe-toggle')

# strobe period
strobePeriod <- htmlDiv(children = list(
  htmlDiv(children = list('strobe pd. (\U{00B5}s)'), className = 'option-name'),
  daqNumericInput(
    id = 'continuous-strobe-period-input',
    value = 1,
    size = 150,
    min = 1,
    max = 100,
    disabled = FALSE
  )
),
id = 'continuous-strobe-period')


lightOptions <- list(list('label' = 'Lamp 1 at 127.0.0.1', 'value' = 'l1'),
                     list('label' = 'Lamp 2 at 127.0.0.1', 'value' = 'l2'))

# light sources
lightSources <- htmlDiv(children = list(
  htmlDiv(children = list('light source'), className = 'option-name'),
  dccDropdown(
    id = 'light-source-input',
    placeholder = "select light source",
    options = lightOptions,
    value = "l2",
    disabled = FALSE
  )
),
id = 'light-source')

controls <- list(intTime, nscansAvg, strobeEnable, strobePeriod, lightSources)


############################
# Layout
############################

pageLayout <- list(htmlDiv(id = 'page', children = list(

  # banner
  htmlDiv(
    id = 'logo',
    title = 'Dash by Plotly',
    style = list(
      'position' = 'absolute',
      'left' = '10px',
      'top' = '10px',
      'zIndex' = 100
    ),
    children = list(htmlA(
      htmlImg(src = "/assets/logo-white.png",
              style = list('height' = '50px')),
      href = "https://plot.ly/dash"
    ))
  ),

  # plot
  htmlDiv(id = 'graph-container',
          children = list(htmlDiv(
            children = list(
              htmlDiv(id = 'graph-title',
                      children = list("ocean optics USB2000+")),
              dccGraph(id = 'spec-readings', animate = TRUE),
              dccInterval(
                id = 'spec-reading-interval',
                interval = 1 * 1000,
                n_intervals = 0,
                max_intervals = 300 # stop after 5 mins.
                # otherwise server has to handle callbacks for idle app
              )
            )
          ))),

  # power button
  htmlDiv(id = 'power-button-container',
          title = 'Turn the power on to begin viewing the data and controlling the spectrometer.',
          children = list(
            daqPowerButton(
              id = 'power-button',
              size = 50,
              color = colors[['accent']],
              on = TRUE
            )
          )),

  # status box
  htmlDiv(
    id = 'status-box',
    children = list(

      # light intensity
      htmlDiv(className = 'status-box-title',
              children = list("light intensity")),
      htmlDiv(
        id = 'light-intensity-knob-container',
        title = 'Controls the intensity of the light source, if any.',
        children = list(
          daqKnob(
            id = 'light-intensity-knob',
            size = 110,
            color = colors[['accent']],
            value = 0
          )
        )
      ),

      # autoscale
      htmlDiv(className = 'status-box-title',
              children = list("autoscale plot")),
      htmlDiv(
        id = 'autoscale-switch-container',
        title = 'Controls whether the plot automatically resizes to fit the spectra.',
        children = list(
          daqBooleanSwitch(id = 'autoscale-switch',
                           on = TRUE,
                           color = colors[['accent']])
        )
      ),

      # submit button
      htmlDiv(
        id = 'submit-button-container',
        title = 'Sends all of the control values below the graph to the spectrometer.',
        children = list(
          htmlButton(
            'update',
            id = 'submit-button',
            n_clicks = 0,
            n_clicks_timestamp = 0
          )
        )
      ),

      # displays whether the parameters were successfully changed
      htmlDiv(
        id = 'submit-status',
        title = 'Contains information about the success or failure of your commands.'
      )
    )
  ),

  # all controls
  htmlDiv(
    id = 'controls',
    title = 'All of the spectrometer parameters that can be changed.',
    children = controls
  ),

  # hidden-div integration time
  htmlDiv(id = 'hidden-div-int-time',
          children = list("1000"),
          hidden = TRUE
  ),

  # hidden-div light source
  htmlDiv(id = 'hidden-div-light-source',
          children = list("l2"),
          hidden = TRUE
  ),

  # about the app
  htmlDiv(id = 'infobox',
          children = list(
            htmlDiv("about this app",
                    id = 'infobox-title'),
            dccMarkdown(
            '
            This is a demo app created to act as an interface for an Ocean Optics
            spectrometer. The options above are used to control various
            properties of the instrument; the integration time, the number of
            scans to average over, the strobe and strobe period, and the
            light source.

            Clicking \"Update\" after putting in the desired settings will
            result in them being sent to the device. A status message
            will appear below the button indicating which commands, if any,
            were unsuccessful; below the unsuccessful commands, a list of
            successful commands can be found.

            (Note that the box containing the status information is
            scrollable.)


            The dial labelled \"light intensity\" will affect the current
            selected light source, if any. The switch labelled \"autoscale
            plot\" will change the axis limits of the plot to fit all of the
            data. Please note that the animations and speed of the graph will
            improve if this feature is turned off, and that it will not be
            possible to zoom in on any portion of the plot if it is turned
            on.
            '
            )
          )
        )
)))

app$layout(htmlDiv(id = 'main', children = pageLayout))


############################
# Helper Variables & FUNs
############################

# input list control elements
controlValues <- lapply(1:length(controls), function(i) {
  input(id = controls[[i]]$props$children[[2]]$props$id,
        property = ifelse(
          length(controls[[i]]$props$children[[2]]$props$value) > 0,
          "value",
          "on"
        ))
})

# state list control elements
controlValuesState <- lapply(1:length(controls), function(i) {
  state(id = controls[[i]]$props$children[[2]]$props$id,
        property = ifelse(
          length(controls[[i]]$props$children[[2]]$props$value) > 0,
          "value",
          "on"
        ))
})

# generates intensities
SampleSpectrum <- function(scale, knobIntensity, x) {
  scale * (exp(-1 * ( (x - 500) / 5) ** 2) +
             0.01 * runif(length(x))) + knobIntensity * 10
}


############################
# Callbacks
############################

# disable/enable the update button depending on whether options have changed
app$callback(
  output = list(id = "submit-button", property = "style"),
  params = c(controlValues,
             list(input(id = "submit-button", "n_clicks_timestamp"))),

  function(...) {
    args <- list(...)

    argsTime <- args[[length(args)]]
    now <- as.numeric(Sys.time()) * 1000

    disabled <- list(
      'color' = colors[['accent']],
      'backgroundColor' = colors[['background']],
      'cursor' = 'not-allowed'
    )

    enabled <- list(
      'color' = colors[['background']],
      'backgroundColor' = colors[['accent']],
      'cursor' = 'pointer'
    )

    # if the button was recently clicked (less than a second ago), then
    # it's safe to say that the callback was triggered by the button; so
    # we have to "disable" it
    if (now - argsTime < 500 && argsTime > 0) {
      return(disabled)
    } else {
      return(enabled)
    }
  }
)


# disable/enable controls
app$callback(
  output = list(id = "controls", property = "children"),
  params = list(input(id = "power-button", property = "on")),

  function(pwr_on) {
    lapply(1:length(controls),
           function(i) {
             controls[[i]]$props$children[[2]]$props$disabled <- !pwr_on
             controls[[i]]
           })
  }
)


# disable/enable intensity knob
app$callback(
  output = list(id = "light-intensity-knob", property = "disabled"),
  params = list(
    input(id = "power-button", property = "on"),
    input(id = "light-source-input", property = "value")
  ),

  function(pwr, lsi) {
    return(!(pwr && !is.null(lsi[[1]]) && lsi != ""))
  }
)


# send user-selected options to spectrometer
app$callback(
  output = list(id = "submit-status", property = "children"),
  params = c(list(input(id = "submit-button", "n_clicks")),
             controlValuesState,
             list(state(id = "power-button", "on"))),

  function(...) {
    args <- list(...)

    # handle `NULL` args`
    nullCheck <- sapply(args, class)
    if ("list" %in% nullCheck) {
      args[[which(nullCheck == "list")]] <- "NULL"
    }

    # power-button off case
    if (!args[[length(args)]]) {
      return (
        list(
          "Press the power button to the top-right of the app,
  then press the \"update\" button above to apply your options to the spectrometer."
        )
      )
    }

    # if submit-button clicked
    if (args[[1]] > 0) {
      sumList <- list()
      nextLine <- list(htmlBr(), htmlBr())

      if ("l1" %in% args ) {
        sumList <- list(
          "The following parameters were not successfully updated:")
        sumList <- c(sumList, nextLine)
        sumList <- c(sumList, list("LIGHT SOURCE: Lamp not found."))
        sumList <- c(sumList, nextLine, list(htmlHr(), htmlBr()))
      }

      sumList <-  c(sumList,
                    list("The following parameters were successfully updated:"),
                    nextLine)

      # append updated controls to sumList
      for (i in 1:length(controls)) {

        if (!(args[[i + 1]] == "l1")) {
          argName <- controls[[i]]$props$children[[1]]$props$children[[1]]
          sumList <- c(sumList,
                       toupper(argName),
                       list(": "),
                       list(as.character(args[[i + 1]])),
                       list(htmlBr()))
        }
      }
      return(sumList)
    }
  }
)


# store `int. time` for plot after submit button clicked
app$callback(
  output = list(id = "hidden-div-int-time", property = "children"),
  params = list(
    input(id = "submit-button", "n_clicks_timestamp"),
    state(id = "integration-time-input", property = "value"),
    state(id = "power-button", property = "on")
  ),

  function(clicks, value, on) {
    if (on) {
      list(value)
    } else {
      list(1000)
    }
  }
)


# store `light source` for plot after submit button clicked
app$callback(
  output = list(id = "hidden-div-light-source", property = "children"),
  params = list(
    input(id = "submit-button", "n_clicks_timestamp"),
    state(id = "light-source-input", property = "value"),
    state(id = "power-button", property = "on")
  ),

  function(clicks, value, on) {
    if (on) {
      # Handle the case when value is removed
      if ( is.list(value) ) {
        list("NULL")
      } else {
      list(value)
      }
    } else {
      list("l2")
    }
  }
)


# update the plot
app$callback(
  output = list(id = "spec-readings", property = "figure"),
  params = list(
    input(id = "spec-reading-interval", property = "n_intervals"),
    state(id = "power-button", property = "on"),
    state(id = "autoscale-switch", property = "on"),
    state(id = "light-intensity-knob", property = "value"),
    state(id = "hidden-div-int-time", property = "children"),
    state(id = "hidden-div-light-source", property = "children")
  ),

  function(n, pwr, auto_range, knob_intensity, scale, lsi) {

    lsi <- unlist(lsi)
    knobIntensity <- ifelse(lsi == "l2", as.numeric(unlist(knob_intensity)), 0)

    scale <- as.numeric(unlist(scale))

    wavelengths <- seq(400, 900, length.out = 5000L)

    if (pwr) {
      intensities <- SampleSpectrum(scale = scale, knobIntensity = knobIntensity, wavelengths)
    } else {
      intensities <- rep(0, length(wavelengths))
    }

    # start creating layout variables
    xAxis <- list(
      'title' = 'Wavelength (nm)',
      'titlefont' = list(
        'family' = 'Helvetica, sans-serif',
        'color' = colors[['secondary']]
      ),
      'tickfont' = list(
        'color' = colors[['tertiary']]
      ),
      'dtick' = 100,
      'color' = colors[['secondary']],
      'gridcolor' = colors[['grid-colour']]
    )

    yAxis <- list(
      'title' = 'Intensity (AU)',
      'titlefont' = list(
        'family' = 'Helvetica, sans-serif',
        'color' = colors[['secondary']]
      ),
      'tickfont' = list(
        'color' = colors[['tertiary']]
      ),
      'color' = colors[['secondary']],
      'gridcolor' = colors[['grid-colour']]
    )

    # auto-scale feature
    if (pwr && auto_range) {
      xAxis[['range']] <- list(min(wavelengths),
                               max(wavelengths))
      yAxis[['range']] <- list(min(intensities),
                               max(intensities))

    }

    layout <- list(
      height = 600,
      font = list('family' = 'Helvetica Neue, sans-serif',
                  'size' = 12),
      margin = list('t' = 20),
      titlefont = list(
        'family' = 'Helvetica, sans-serif',
        'color' = colors[['primary']],
        'size' = 26),
      xaxis = xAxis,
      yaxis = yAxis,
      paper_bgcolor = colors[['background']],
      plot_bgcolor = colors[['background']]
      )

    return(
      figure <- list(data = list(
        list(
          x = as.list(wavelengths),
          y = as.list(intensities),
          name = 'Spectrometer readings',
          mode = 'lines',
          line = list('width' = 1, 'color' = colors[['accent']])
        )
      ),
      layout = layout)
    )
  }
)

if (!appName == ""){
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(debug = TRUE)
}
