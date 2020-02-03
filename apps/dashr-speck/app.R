library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashBio)
library(dashDaq)
library(dashTable)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

header_colors <- function() {
  list(
    bg_color = '#5673FA',
    font_color = 'white'
  )
}

DATAPATH <- './sample_data/speck_'

##### main code ####
app <- Dash$new()

default_sliders <- list(
  htmlDiv(
    className = 'app-controls-block',
    children = list(
      htmlDiv(
        "Atom radius",
        className='app-controls-name'
      ),
      dccSlider(
        id = 'speck-atom-radius',
        className = 'control-slider',
        max = 1,
        step = 0.01,
        value = 0.6,
        updatemode = 'drag'
      )
    )
  ),
  htmlDiv(
    className = 'app-controls-block',
    children = list(
      htmlDiv(
        "Relative atom radius",
        className = 'app-controls-name'
      ),
      dccSlider(
        id = 'speck-relative-atom-radius',
        className = 'control-slider',
        max = 1,
        step = 0.01,
        value = 1.0,
        updatemode = 'drag'
      )
    )
  ),
  htmlDiv(className='app-controls-block', 
          children = list(
            htmlDiv(
              "Ambient occlusion",
              className = 'app-controls-name'
            ),
            dccSlider(
              id = 'speck-ao',
              className = 'control-slider',
              max = 1,
              step = 0.01,
              value = 0.75
            )
          )
  ),
  htmlDiv(className = 'app-controls-block', 
          children = list(
            htmlDiv(
              "Brightness",
              className='app-controls-name'
            ),
            dccSlider(
              id = 'speck-brightness',
              className = 'control-slider',
              max = 1,
              step = 0.01,
              value = 0.5,
              updatemode = 'drag'
            )
          )
  ),
  htmlDiv(className = 'app-controls-block', 
          children = list(
            htmlDiv(
              "Outline",
              className = 'app-controls-name'
            ),
            dccSlider(
              id = 'speck-outline',
              className = 'control-slider',
              max = 1,
              step = 0.01,
              value = 0.0,
              updatemode = 'drag'
            )
          )
  ),
  htmlHr(),
  dccChecklist(
    id = 'speck-show-hide-bonds',
    options = list(
      list(
        label = 'Show bonds',
        value = 'True'
      )
    ),
    value = list()
  ),
  htmlDiv(className = 'app-controls-block', 
          children = list(
            htmlDiv(
              'Bond scale',
              className = 'app-controls-name'
            ),
            dccSlider(
              id = 'speck-bond-scale',
              className = 'control-slider',
              max = 1,
              step = 0.01,
              value = 0.5,
              updatemode = 'drag'
            )
          )
  )
)

tab1 <-   dccTab(
  label = 'About',
  value = 'what-is',
  children = htmlDiv(
    className = 'control-tab',
    children = list(
      htmlH4(className = 'what-is', children = 'What is Speck?'),
      htmlP('Speck is a WebGL-based molecule renderer. By
            using ambient occlusion, the resolution of
            the rendering does not suffer as you zoom in.'),
      htmlP('You can toggle between molecules using the menu under the
            \"Data\" tab, and control parameters related to
            the appearance of the molecule in the \"View\" tab.
            These parameters can be controlled at a low level
            with the sliders provided, or preset views can be
            applied for a higher-level demonstration of changing
            atom styles and rendering.')
    )
  )
)

tab2 <- dccTab(
  label = 'Data',
  value = 'datasets',
  children = htmlDiv(
    className = 'control-tab',
    children = list(
      htmlDiv(
        className = 'app-controls-block',
        children = list(
          htmlDiv(
            className = 'app-controls-name',
            children = 'Preloaded'
          ),
          dccDropdown(
            id = 'speck-molecule-dropdown',
            className = 'speck-dropdown',
            options = list(
              list(
                label = 'DNA',
                value = paste0(DATAPATH, 'dna.xyz')
              ),
              list(
                label = 'Caffeine',
                value = paste0(DATAPATH, 'caffeine.xyz')
              ),
              list(
                label = 'Methane',
                value = paste0(DATAPATH, 'methane.xyz')
              ),
              list(
                label = 'Testosterone',
                value = paste0(DATAPATH, 'testosterone.xyz')
              ),
              list(
                label = 'Gold nanoparticle',
                value = paste0(DATAPATH, 'au.xyz')
              ),
              list(
                label = 'Thiolated gold nanoparticle',
                value = paste0(DATAPATH, 'au_thiol.xyz')
              ),
              list(
                label = 'Benzene',
                value = paste0(DATAPATH, 'benzene.xyz')
              ),
              list(
                label = 'Protein (4E0O)',
                value = paste0(DATAPATH, '4E0O.xyz')
              )
            ),
            value = paste0(DATAPATH, 'dna.xyz')
          )
        )
      ),
      htmlDiv(id='speck-preloaded-uploaded-alert'),
      dccUpload(
        id = 'speck-file-upload',
        className = 'control-upload',
        children = htmlDiv(list("Drag and drop .xyz files, or click to select files."))
      ),
      htmlBr(),
      htmlA(
        htmlButton(
          'Download sample .xyz data',
          id = 'speck-file-download',
          className = 'control-download'
        ),
        href = 'assets/sample_data/speck_4QCI.xyz',
        download = '4QCI.xyz'
      )
    )
  )
)

tab3 <- dccTab(
  label = 'View',
  value = 'view-options',
  children = htmlDiv(
    className = 'control-tab', 
    children = list(
      dccChecklist(
        id = 'speck-enable-presets',
        options = list(
          list(label = 'Use presets', 
               value = 'True')
        ),
        value = list()
      ),
      htmlHr(),
      htmlDiv(id = 'speck-controls-detailed', 
              children = default_sliders),
      htmlDiv(
        id = 'speck-controls-preset',
        className = 'speck-controls',
        children = list(
          htmlDiv(className = 'app-controls-block', 
                  children = list(
                    htmlDiv(className = 'app-controls-name',
                            children = 'Rendering style'),
                    dccDropdown(
                      id = 'speck-preset-rendering-dropdown',
                      className = 'speck-dropdown',
                      options = list(
                        list(label = 'Default/reset',
                             value = 'default'),
                        list(label = 'Toon',
                             value = 'toon')
                      ),
                      value='default'
                    )
                  )
          ),
          htmlDiv(className = 'app-controls-block', 
                  children = list(
                    htmlDiv(className = 'app-controls-name', children='Atom style'),
                    dccDropdown(
                      id = 'speck-preset-atom-style-dropdown',
                      className = 'speck-dropdown',
                      options = list(
                        list(label = 'Ball-and-stick',
                             value = 'stickball'),
                        list(label = 'Licorice',
                             value = 'licorice')
                      ),
                      value = 'default'
                    )
                  )
          )
        )
      )
    )
  )
)

header <- htmlDiv(
  id = "app-page-header",
  style = list(
    width = "100%",
    background = header_colors()[["bg_color"]],
    color = header_colors()[["font_color"]]
  ),
  children = list(
    htmlA(
      id = "dashbio-logo",
      children = list(
        htmlImg(src='assets/plotly-dash-bio-logo.png', height = '36', width = '190',
                style = list('top' = '10', 'margin-left' = '10px'))
      ),
      href = "/Portal"
    ),
    htmlH2("Speck"),
    htmlA(
      id = "gh-link",
      children = list("View on GitHub"),
      href = "https://github.com/plotly/dash-sample-apps/blob/master/apps/dashr-speck/app.R",
      style = list(color = "white", border = "solid 1px white")
    ),
    htmlImg(
      src = "assets/GitHub-Mark-Light-64px.png"
    )
  )
)

#### dash app layout ####
app$layout(
  header,
  htmlBr(),
  htmlDiv(
    id = 'speck-body',
    className = 'app-body',
    children = list(
      dccLoading(className = 'dashbio-loading', 
                 children = htmlDiv(
                   id = 'speck-container',
                   children = list(
                     dashbioSpeck(
                       id = 'speck',
                       view = list('resolution' = 600, 'zoom' = 0.3),
                       scrollZoom = TRUE
                     )
                   )
                 )
      ),
      htmlDiv(
        id = 'speck-control-tabs',
        className = 'control-tabs',
        children = list(
          dccTabs(
            id = 'speck-tabs',
            value = 'what-is',
            children = list(
              tab1,
              tab2,
              tab3
            )
          )
        )
      ),
      dccStore(
        id='speck-store-preset-rendering',
        data = list()
      ),
      dccStore(
        id='speck-store-preset-atom-style',
        data = list()
      )
    )
  )
)

#### dash app callback ####
app$callback(
  output = list(id = 'speck-controls-detailed', property = 'style'),
  params = list(input(id = 'speck-enable-presets', property = 'value')),
  function(presets_enable){
    list(
      display = if(length(presets_enable) > 0) 'none' else 'inline-block'
    )
  }
)

app$callback(
  output = list(id = 'speck-controls-preset', property = 'style'),
  params = list(input(id = 'speck-enable-presets', property = 'value')),
  function(presets_enable){
    list(
      display = if(length(presets_enable) == 0) 'none' else 'inline-block'
    )
  }
)

app$callback(
  output = list(id = 'speck-molecule-dropdown', property = 'value'),
  params = list(input(id = 'speck-file-upload', property = 'contents'),
                state(id = 'speck-molecule-dropdown', property = 'value')),
  function(upload_contents, current){
    
    if(is.null(unlist(upload_contents))) {
      current
    } else {
      list()
    }
  }
)

app$callback(
  output = list(id = 'speck-preloaded-uploaded-alert', property = 'children'),
  params = list(input(id = 'speck-molecule-dropdown', property = 'value'),
                input(id = 'speck-file-upload', property = 'contents'),
                state(id = 'speck-file-upload', property = 'filename')),
  function(molecule_fname, upload_contents, upload_fname){
    if(!is.null(unlist(molecule_fname)) && !is.null(unlist(upload_contents))) {
      paste0("Warning: you have uploaded a dataset ", 
             upload_fname, 
             "To view the dataset,\nplease ensure that the \"Preloaded\" dropdown has\n
             been cleared.")
    } else ""
  }
)

read_xyz <- function(filepath, 
                     header = FALSE, 
                     skip = 2) {
  textdata <- read.table(
    text = paste0(
      readLines(filepath), collapse="\n"
    ),
    header = header,
    skip = skip,
    col.names = c("symbol", "x", "y", "z"),
    stringsAsFactors = FALSE)
  return(dashTable::df_to_list(textdata))
}

get_xyz <- function(decoded_char, skip = 2) {
  data <- strsplit(decoded_char, "\n")[[1]][-(1:skip)]
  lapply(data, 
         function(x) {
           info <- Filter(f = function(y) y != "",strsplit(x, " ")[[1]])
           if(length(info) != 4) stop("input file should have four columns: symbol, x, y and z")
           list(
             symbol = info[1],
             x = as.numeric(info[2]),
             y = as.numeric(info[3]),
             z = as.numeric(info[4])
           )
         })
}

app$callback(
  output = list(id = 'speck', property = 'data'),
  params = list(input(id = 'speck-molecule-dropdown', property = 'value'),
                input(id = 'speck-file-upload', property = 'contents')),
  function(molecule_fname, upload_contents) {
    if(!is.null(unlist(upload_contents)) && is.null(unlist(molecule_fname))) {
      
      content_string <- unlist(strsplit(upload_contents, ","))[2]
      decoded <- jsonlite::base64_dec(content_string)
      get_xyz(rawToChar(decoded))
    } else if (!is.null(unlist(molecule_fname))) {
      read_xyz(molecule_fname)
    } else list()
  }
)

app$callback(
  output = list(id = 'speck', property = 'view'),
  params = list(
    input(id = 'speck-enable-presets', property = 'value'),
    input(id = 'speck-atom-radius', property = 'value'),
    input(id = 'speck-relative-atom-radius', property = 'value'),
    input(id = 'speck-show-hide-bonds', property = 'value'),
    input(id = 'speck-bond-scale', property = 'value'),
    input(id = 'speck-ao', property = 'value'),
    input(id = 'speck-brightness', property = 'value'),
    input(id = 'speck-outline', property = 'value')
  ), function(presets_enabled,
              atom_radius,
              relative_atom_radius,
              show_bonds,
              bond_scale,
              ambient_occlusion,
              brightness,
              outline) {
    list(
      atomScale = atom_radius,
      relativeAtomScale = relative_atom_radius,
      bonds = (length(unlist(show_bonds)) > 0),
      bondScale = bond_scale,
      ao = ambient_occlusion,
      brightness = brightness,
      outline = outline
    )
  }
)

app$callback(
  output = list(id = 'speck-store-preset-rendering', property = 'data'),
  params = list(input(id = 'speck-preset-rendering-dropdown', 
                      property = 'value')),
  function(render) render
)

app$callback(
  output = list(id = 'speck-store-preset-atom-style', property = 'data'),
  params = list(input(id = 'speck-preset-atom-style-dropdown', 
                      property = 'value')),
  function(atomstyle) atomstyle
)

app$callback(
  output = list(id = 'speck', property = 'presetView'),
  params = list(
    input(id = 'speck-store-preset-atom-style', property = 'modified_timestamp'),
    input(id = 'speck-store-preset-rendering', property = 'modified_timestamp'),
    state(id = 'speck-preset-rendering-dropdown', property = 'value'),
    state(id = 'speck-preset-atom-style-dropdown', property = 'value')
  ),
  function(atomstyle_ts, 
           render_ts,
           render, 
           atomstyle) {

    if(is.null(unlist(atomstyle_ts)) && is.null(unlist(render_ts))) {
      return('default')
    }
    if(!is.null(unlist(atomstyle_ts)) && is.null(unlist(render_ts))) {
      atomstyle
    } else if(is.null(unlist(atomstyle_ts)) && !is.null(unlist(render_ts))) {
      render
    } else {
      if(render_ts > atomstyle_ts || is.null(unlist(atomstyle))) {
        if(is.null(unlist(render))) 'default' else render
      } else atomstyle
    }
  }
)

app$callback(
  output = list(id = 'speck-preset-atom-style-dropdown', property = 'value'),
  params = list(
    input(id = 'speck-preset-rendering-dropdown', property = 'value'),
    state(id = 'speck-preset-atom-style-dropdown', property = 'value')),
    function(render, current) {
      if(render == 'default') list() else current
    }
)
#### dash app run ####
app$run_server(host = "0.0.0.0", port=8050)
