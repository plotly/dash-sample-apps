appName <- Sys.getenv("DASH_APP_NAME")

if (appName != "") {
 pathPrefix <- sprintf("/%s/", appName)

 Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
            DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}


# Load Necessary Packages
library('dash')
library('dashCoreComponents')
library('dashHtmlComponents')
library('plyr')
library('plotly')
library('dashBio')
library('dashDaq')

#################################################################################

app <- Dash$new()
  
  
  
  
  listofchromosomes = (sprintf("%d",seq(1:22)))
  
  listofchromosomes = c(listofchromosomes, "X")
  
  listofchromosomes = c(listofchromosomes, "Y")
  
  
chromosomes <- c(as.character(seq(1:22)), "X", "Y")

listOfOptions <- lapply(chromosomes, function(x) {
  list(label = x, value = x)
})
  
  chromosome_div <- function(id_tag = 'chr',
                             name_tag = 'Chr',
                             startone = 0,
                             stopone = 1,
                             starttwo = 0,
                             stoptwo = 1) {
    return(htmlDiv(list(
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = sprintf('%s Start-one', name_tag)),
        dccInput(
          id=sprintf('%s-startone', id_tag),
          placeholder= sprintf('%s Startone', name_tag),
          type = 'number',
          value = startone,
          className = 'ideogram-homology-inputs'
        )
      )),
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children= sprintf('%s Stop-one', name_tag)),
        dccInput(
          id=sprintf('%s-stopone', id_tag),
          placeholder = 'Enter chromosomes',
          type = 'number',
          value = stopone,
          className = 'ideogram-homology-inputs'
        )
      )),
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children= sprintf('%s Start-two', name_tag)),
        dccInput(
          id=sprintf('%s-starttwo', id_tag),
          placeholder =sprintf('%s Starttwo', name_tag),
          type = 'number',
          value = starttwo,
          className = 'ideogram-homology-inputs'
        )
      )),
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children= sprintf('%s Stop-two', name_tag)),
        dccInput(
          id=sprintf('%s-stoptwo', id_tag),
          placeholder = 'Enter chromosomes',
          type = 'number',
          value = stoptwo,
          className = 'ideogram-homology-inputs'
        )
      ))
      
    )))
  }
  
  
  
  options = list(
    'custom' = list(
      htmlH4('Organism'),
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Species'),
        dccDropdown(
          className = 'ideogram-dropdown',
          id = 'organism-change',
          options = list(
            
            list('label' = 'Human', 'value' = 'human'),
            
            list('label' = 'Drosophila-Melanogaster', 'value' = 7227),
            
            list('label' = 'Zea-mays', 'value' = 4577),
            
            list('label' = 'Pan-troglodytes', 'value' = 9598),
            
            list('label' = 'Macaca-fascicularis', 'value' = 9541)
          ) , value = 'human'
        )
      )),
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Sex'),
        daqToggleSwitch(
          id='sex-switch',
          color = '#230047',
          label = list('m', 'f'),
          size = 35,
          labelPosition = 'bottom',
          value = FALSE
        )
      )),
      
      htmlDiv(id='ideogram-resolution-option', children = list(
        htmlDiv(className = 'app-controls-block', children = list(
          htmlDiv(className = 'app-controls-name', children = 'Resolution'),
          
          dccDropdown(
            className = 'ideogram-dropdown',
            id = 'resolution-select',
            options = list(
              list('label' = '500 bphs', 'value' = 500),
              list('label' = '550 bphs', 'value' = 550),
              list('label' = '650 bphs', 'value' = 850)
              
              
            ) , value = 550
          )
        ))
      )),
      
      htmlHr(),
      
      htmlH4('Labels'),
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Band'),
        daqToggleSwitch(
          id='bandlabel-switch',
          color = '#230047',
          label = list('off', 'on'),
          size = 35,
          labelPosition = 'bottom',
          value = TRUE
        )
      )),
      
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Chromosome'),
        daqToggleSwitch(
          id='chromlabel-switch',
          color = '#230047',
          label = list('off', 'on'),
          size = 35,
          labelPosition = 'bottom',
          value = TRUE
        )
      )),
      
      htmlHr(),
      
      htmlH4('Chromosome display'),
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Orientation'),
        dccDropdown(
          className = 'ideogram-dropdown',
          id = 'orientation-switch',
          options = list(
            list('label' = 'Vertical', 'value' = 'vertical'),
            list('label' = 'Horizontal', 'value' = 'horizontal')
          ), value = 'horizontal'
        )
      )),
      
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Rotatable'),
        daqToggleSwitch(
          id='rotatable-switch',
          color = '#230047',
          label = list('off', 'on'),
          size = 35,
          labelPosition = 'bottom',
          value = TRUE
        )
      )),
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Margin'),
        dccSlider(
          id = 'chr-margin-input',
          className = 'ideogram-slider',
          value = 10,
          min=2, max=50,
        ),
        
        htmlDiv(className = 'app-controls-name', children = 'Height'),
        dccSlider(
          id='chr-height-input',
          className = 'ideogram-slider',
          min = 100, max = 700,
          value = 300
        ),
        
        htmlDiv(className = 'app-controls-name', children = 'Width'),
        dccSlider(
          id='chr-width-input',
          className = 'ideogram-slider',
          min = 5, max = 30,
          value = 8
        )
        
        
      )),
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Fully Banded'),
        daqToggleSwitch(
          id='fullband-switch',
          color = '#230047',
          label = list('off', 'on'),
          size = 35,
          labelPosition = 'bottom',
          value = TRUE
        )
      )),
      
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Select Chromosomes'),
        dccDropdown(
          className = 'ideogram-dropdown',
          id = 'displayed-chromosomes',
          options = listOfOptions,
          multi = TRUE,
          value = listofchromosomes
        )
      ))
      
  
      
    ),
    
    'homology' = list(
      htmlDiv(className = 'app-controls-name', children = 'Chromosomes:' ),
      dccDropdown(
        className = 'ideogram-dropdown',
        id = 'chr-select-1',
        options = listOfOptions , value = '1'
      ),
      
      dccDropdown(
        className = 'ideogram-dropdown',
        id = 'chr-select-2',
        value = '2',
        options = listOfOptions
      ),
      
      chromosome_div(
        id_tag = 'chrone',
        name_tag = 'Chr 1',
        startone=50000,
        stopone=900000,
        starttwo=155701383,
        stoptwo=156030895
      ),
      
      chromosome_div(
        id_tag = 'chrtwo',
        name_tag = 'Chr 2',
        startone=10001,
        stopone=2781479,
        starttwo=56887903,
        stoptwo=57217415
      )
    ),
    
    'brush' = list(
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', id = 'brush-control-name', children = 'Chromosome:'),
        dccDropdown(
          className = 'ideogram-dropdown',
          id = 'chr-brush',
          options = listOfOptions, value = 'X'
        )
      )),
      
      htmlHr(),
      
      htmlDiv(
        id = 'brush-data',
        children = list(
          htmlH4('Brush Data'),
          'Start: ',
          htmlSpan(
            '',
            id='brush-print-start',
            style = list('color' = '#0D76Bf')
          ),
          htmlBr(),
          
          'Extent: ',
          
          htmlSpan(
            '',
            id='brush-print-extent',
            style = list('color' = '#0D76BF')
          ),
          
          htmlBr(),
          'End: ',
          htmlSpan(
            '',
            id='brush-print-end',
            style = list('color' = '#0D76BF')
          )
        ),
        className = 'ideogram-databox-parameters'
      )
    ),
    
    'annotations' = list(
      htmlH4('Hover data'),
      htmlDiv(children = list(
        htmlSpan(
          'None',
          id='annote-data',
          style = list('color' = '#0D76BF')
        )
      ), className = 'ideogram-databox-parameters'),
      
      htmlHr(),
      
      htmlH4('Appearance'),
      
      htmlDiv(className = 'app-controls-block', children = list(
        htmlDiv(className = 'app-controls-name', children = 'Type: '),
        dccDropdown(
          className = 'ideogram-dropdown',
          id= 'annotation-select',
          options = list(
            list('label' = 'Histogram', 'value' = 'histogram'),
            list('label' = 'Overlay-1', 'value' = 'overlay-1'),
            list('label' = 'Overlay-2', 'value' = 'overlay-2')
          ), value = 'histogram'
        )
      )),
      
      htmlDiv(id = 'ideogram-histo-options', children = list(
        htmlDiv(className = 'app-controls-block', children = list(
          htmlDiv(className = 'app-controls-name', children = 'Color:'),
          dccInput(
            id='color-input',
            className = 'ideogram-annot-inputs',
            type = 'text',
            value = '#FF0000'
          )
        )),
        
        htmlDiv(className = 'app-controls-block', children=list(
          htmlDiv(className = 'app-controls-name', children = 'Bar width:'),
          dccSlider(
            id='bar-input',
            className = 'ideogram-slider',
            value = 3,
            min = 1,
            max = 10
          )
        )),
        
        htmlDiv(className = 'app-controls-block', children=list(
          htmlDiv(className = 'app-controls-name', children = 'Height:'),
          dccSlider(
            id='height-input',
            className = 'ideogram-slider',
            value = 3,
            min = 1,
            max = 10
          )
        )),
        
        htmlDiv(className = 'app-controls-block', children = list(
          htmlDiv(className = 'app-controls-name', children = 'Orientation:'),
          dccDropdown(
            className = 'ideogram-dropdown',
            id= 'orientation-anote',
            options = list(
              list('label' = 'Vertical', 'value'='vertical'),
              list('label' = 'Horizontal', 'value' = 'horizontal')
            ), value = 'horizontal'
          )
        ))
        
        
      ))
    )
  )
  
  ideograms_initial = list(
    'custom' = list(
      id = 'ideo-custom',
      dataDir = 'https://unpkg.com/ideogram@1.3.0/dist/data/bands/native/',
      orientation = 'horizontal',
      organism = 'human',
      chrHeight = 300,
      chrMargin = 10,
      chrWidth = 8,
      rotatable = TRUE
    ),
    
    'homology' = list(
      id = 'ideo-homology',
      showBandLabels = TRUE,
      showChromosomeLabels = TRUE,
      showFullyBanded = TRUE,
      fullChromosomeLabels = TRUE,
      chrHeight = 400,
      chrMargin = 200,
      rotatable = FALSE,
      perspective = 'comparative'
    ),
    
    'brush' = list(
      id = 'brush-ideo',
      dataDir = 'https://unpkg.com/ideogram@1.3.0/dist/data/bands/native/',
      organism = 'human',
      chromosomes = list('1'),
      brush = 'chr1:1-10000000',
      chrHeight = 900,
      resolution = 550,
      orientation = 'horizontal'
    ),
    
    'annotations' = list(
      id='ideo-annotations',
      dataDir = 'https://unpkg.com/ideogram@1.3.0/dist/data/bands/native/',
      organism = 'human',
      assembly = 'GRCh37',
      showBandLabels = TRUE,
      chrHeight = 275,
      chrMargin = 28,
      rotatable = TRUE,
      filterable = TRUE,
      className = 'ideogram-cluster'
    )
  )
  
  #####################################################################################################
  
  #Layout
  
  
  ideolayout <- htmlDiv(list(
    htmlDiv(id = 'ideogram-body', className = 'app-body', children = list(
      dccLoading(className = 'dashbio-loading', children = htmlDiv(id = 'ideogram-container')),
      htmlDiv(className = 'control-tabs', children = list(
        dccTabs(id='ideogram-control-tabs', value = 'what-is', children=list(
          dccTab(
            label = 'About',
            value = 'what-is',
            children = htmlDiv(className = 'control-tab', children = list(
              htmlH4(className = 'what-is', children = 'What is Ideogram'),
              htmlP('Ideogram is a tool used to schematically represent
                      chromosomes. Bands on the chromosomes can show the locations
                      of specific genes.'),
              htmlP('In the "View" tab, you can choose to interact with several
                      different features of the Ideogram component. You can customize the 
                      appearance of the ideogram, as well as choose a different organism
                      to display, under the "Custom" option. The homology, brush, and 
                      annotation features are demonstrated under the corresponding options.', id = 'text')
            ))
          ),
          
          dccTab(
            label = 'View',
            value = 'view',
            children = htmlDiv(className = 'control-tab', children = list(
              htmlDiv(id='ideogram-feature-select', children = list(
                htmlDiv(className = 'app-controls-block', children = list(
                  htmlDiv(className = 'app-controls-name', children = 'View feature:'),
                  
                  dccDropdown(
                    className = 'ideogram-dropdown',
                    id='ideogram-feature-dropdown',
                    options = list(
                      list('label' = 'Customizability', 'value' = 'custom'),
                      list('label' = 'Homology', 'value' = 'homology'),
                      list('label' = 'Brush', 'value' = 'brush'),
                      list('label' = 'Annotations', 'value' = 'annotations')
                      
                    ), clearable = FALSE, value = 'custom'
                  )
                ))
              )),
              
              htmlHr(),
              htmlDiv(
                id = 'ideogram-feature-view-options'
              )
            ))
          )
        ))
      )),
      
      
  
      
      htmlDiv(do.call(dashbioIdeogram, ideograms_initial[['custom']]), style = list('display' = 'none')),
      
      
      
      dccStore(id='ideo-custom-data', data=ideograms_initial[['custom']]),
      dccStore(id='ideo-homology-data', data=ideograms_initial[['homology']]),
      dccStore(id='brush-ideo-data', data=ideograms_initial[['brush']]),    # Might be brush-ideo-data
      dccStore(id='ideo-annotations-data', data=ideograms_initial[['annotations']])
    )
    )
  ))
  
  header <- htmlDiv(children = list(
    htmlSpan(list(htmlA(
      htmlImg(src='assets/plotly-dash-bio-logo.png', height = '36', width = '180',
              style = list('top' = '10', 'margin-left' = '10px', 'margin-right' = '15px')),
      href='https://dash-bio.plotly.host/Portal/'))),
    
      "Ideogram",
  
      
      htmlA(href="http://github.com/plotly/dash-bio/blob/master/tests/dashbio_demos/app_ideogram.py",
            style = list("color" = "white",
                         "border" = "1px solid white",
                         "font-size" = "18px",
                         "position" = "relative",
                         "float" = "right",
                         "padding-top" = "10px",
                         "padding-left" = "20px",
                         "padding-right" = "20px",
                         "padding-bottom" = "10px",
                         "font-weight" = "100",
                         "margin-right" = "30px",
                         "top" = "17px",
                         "text-decoration" = "none"
                         
                         ),
            children = "View on Github"),
    
    
    htmlImg(src = "assets/github.png", width = '50', height = '50', style = list(
      "position" = "relative", "float" = "right",
      "margin-right" = "10px",
      "top" = "17px"
    ))
    ), className = 'ideogramheader')
  
  
  
  
  app$layout(htmlDiv(list(
      header,
      ideolayout,
      htmlDiv(list(
      "                "
      ), className = "blankspace"),
      daqToggleSwitch(
        id = 'daq-light-dark-theme',
        color = '#FDFCFF',
        labelPosition = 'bottom',
        label = list('Light', 'Dark'),
        style = list('width' = '250px', 'margin' = 'auto', 'display' = 'none'),
        value = FALSE
      )
    ), style = list('background-color' = '#FDFCFF')
  ))
  
  
  
  #####################################################################################################
  
  # Callbacks
  
  
  
  # Select which chromosomes to show for the customizable feature.
  
   app$callback(
      output = list(id = 'ideo-custom', property = 'chromosomes'),
      params = list(input(id = 'displayed-chromosomes', property = 'value')),
  
      update_ideogram <- function(value) {
        return(value)
      }
    )
  
  # Updates the plot that is shown in the ideogram container. 
  
  app$callback(
    output = list(id = 'ideogram-container', property = 'children'),
    params = list(
      input(id = 'ideo-custom-data', property = 'data'),
      input(id = 'ideo-homology-data', property = 'data'),
      input(id = 'brush-ideo-data', property = 'data'),
      input(id = 'ideo-annotations-data', property = 'data'),
      input(id = 'ideogram-feature-dropdown', property = 'value')
    ),
    
    update_ideogram <- function(
      ideo_custom,
      ideo_homology,
      ideo_brush,
      ideo_annotations,
      selected_ideo
    ) {
      
      ideograms = list(
        'custom' = ideo_custom,
        'homology' = ideo_homology,
        'brush' = ideo_brush,
        'annotations' = ideo_annotations
      )
      
      
      return(do.call(dashbioIdeogram, ideograms[[selected_ideo]]))
      
      
    }
  )
  
  
  
  # Updates which options are shown based on the selected feature.
  
  app$callback(output = list(id = 'ideogram-feature-view-options', property = 'children'),
               params = list(input(id = 'ideogram-feature-dropdown', property = 'value')),
               
               show_options <- function(feature) {
                 return(options[[feature]])
               }
  )
  
  
  # Change resolution for the plot. 
  
  app$callback(
    output = list(id = 'ideogram-resolution-option', property = 'style'),
    params = list(input(id = 'organism-change', 'value')),
    
    show_hide_resolution <- function(organism) {
      if (organism != 'human') {
        resolution <- list('display' = 'none')
      }
      
      else {
        resolution <- list('display' = 'block')
      }
      return(resolution)
    }
  )
  
  
  
  
  app$callback(
    output = list(id = 'annote-data', property = 'children'),
    params = list(
      input(id = 'ideo-annotations' , property = 'annotationsData')
    ),
    
    annote_data <- function(data) {
      if (exists('data') == FALSE) {
        data = 'None'
      }
      
    data = gsub('<br>', ' ', data)
    return(data)
    }
  )
  
  # Change ideogram annotation:
  
  app$callback(
    output = list(id = 'ideo-annotations-data', property = 'data'),
    params = list(
      input(id = 'annotation-select', property = 'value'),
      input(id = 'bar-input', property = 'value'),
      input(id = 'orientation-anote', property = 'value'),
      input(id = 'color-input', property = 'value'),
      input(id = 'height-input', property = 'value'),
      state(id = 'ideo-annotations-data', property = 'data')
    ),
    
    change_annotation_ideogram <- function(
      annotation_select,
      bar_input,
      orientation_input,
      color_input,
      height_input,
      current
    ) {
      if (!exists('current')) {
        current = ideograms_initial[['annotations']]
        
        annotations_layout = ""
        annotations_path = NULL
        annotations_assembly = NULL
        annotation_tracks = NULL
        annotation_height = height_input
        
        
        
        if (annotation_select == 'overlay-1') {
          annotations_layout = 'overlay'
          annotations_path = 'assets/10_virtual_cnvs.json'
        }
        
        else if (annotation_select == 'histogram') {
          annotations_path = 'assets/SRR562646.json'
          annotations_assembly = 'GRCh37'
        }  
        
        else if (annotations_select == 'overlay-2') {
          annotations_path = 'assets/1000_virtual_snvs.json'
          annotation_tracks  = list(
            list('id' = 'pathogenicTrack',
                 'displayName' = 'Pathogenic',
                 'color' = '#F00',
                 'shape' = 'triangle'
            ),
            list('id' = 'uncertainSignificanceTracks',
                 'displayName' = 'Uncertain Significance',
                 'color' = '#ccc',
                 'shape' = 'triangle'
            ),
            list('id' = 'benignTrack',
                 'displayName' = 'Benign',
                 'color' = '#8D4',
                 'shape' = 'triangle'
            )
          )
        }  
        
        current[['annotationsLayout']] = annotations_layout
        current[['barWidth']] = bar_input
        current[['annotationsPath']] = annotations_path
        current[['assembly']] = annotations_assembly
        current[['annotationTracks']] = annotation_tracks
        current[['annotationHeight']] = annotation_height
        current[['annotationsColor']] = color_input
        current[['orientation']] = orientation_input
      } 
      
      else {
        annotations_layout = ""
        annotations_path = NULL
        annotations_assembly = NULL
        annotation_tracks = NULL
        annotation_height = height_input
        
        
        
        if (annotation_select == 'overlay-1') {
          annotations_layout = 'overlay'
          annotations_path = 'https://eweitz.github.io/ideogram/data/annotations/10_virtual_cnvs.json'
        }
        
        else if (annotation_select == 'histogram') {
          annotations_path = 'https://eweitz.github.io/ideogram/data/annotations/SRR562646.json'
          annotations_assembly = 'GRCh37'
        }  
        
        else if (annotation_select == 'overlay-2') {
          annotations_path = 'https://eweitz.github.io/ideogram/data/annotations/1000_virtual_snvs.json'
          annotation_tracks  = list(
            list('id' = 'pathogenicTrack',
                 'displayName' = 'Pathogenic',
                 'color' = '#F00',
                 'shape' = 'triangle'
            ),
            list('id' = 'uncertainSignificanceTracks',
                 'displayName' = 'Uncertain Significance',
                 'color' = '#ccc',
                 'shape' = 'triangle'
            ),
            list('id' = 'benignTrack',
                 'displayName' = 'Benign',
                 'color' = '#8D4',
                 'shape' = 'triangle'
            )
          )
        }   
        
        
        
        current[['annotationsLayout']] = annotations_layout
        current[['barWidth']] = bar_input
        current[['annotationsPath']] = annotations_path
        current[['assembly']] = annotations_assembly
        current[['annotationTracks']] = annotation_tracks
        current[['annotationHeight']] = annotation_height
        current[['annotationsColor']] = color_input
        current[['orientation']] = orientation_input
      } 
        
     print(current)
      
     return(current)
    }
  )
  
  
  
  
  
  
  
  # Brush Feature Callbacks
  
  
  app$callback(
    output = list(id = 'brush-ideo-data', property = 'data'),
    params = list(
      input(id = 'chr-brush', property = 'value'),
      state(id = 'brush-ideo-data', property = 'data')
    ),
    
    
    change_brush_ideo <- function(brush_value, current) {
      if (is.null('current')) {
        current = ideograms_initial[['brush']]
        current['chromosomes'] = as.character(brush_value)
        current['brush'] = sprintf('chr%s:1-10000000', brush_value)
      }
      else {
        current['chromosomes'] = as.character(brush_value)
        current['brush'] = sprintf('chr%s:1-10000000', brush_value)
      }
      
      print(brush_value)
      
      print(current)
      
      return(current)
    }
    
  )
  
  
  
  app$callback(
    output = list(id = 'brush-print-start', property = 'children'),
    params = list(
      input(id = 'brush-ideo', 'brushData')
    ),
    brush_data_start <- function(brush_data) {
      answer = NA
      
      if (exists('brush_data')) {
        answer = brush_data[['start']]
      }
      
      return(answer)
    }
  )
  
  
  app$callback(
    output = list(id = 'brush-print-end', property = 'children'),
    params = list(
      input(id = 'brush-ideo', 'brushData')
    ),
    brush_data_end <- function(brush_data) {
      answer = NA
      
      if (exists('brush_data')) {
        answer = brush_data[['end']]
      }
      
      return(answer)
    }
  )
  
  
  app$callback(
    output = list(id = 'brush-print-extent', property = 'children'),
    params = list(
      input(id = 'brush-ideo', 'brushData')
    ),
    brush_data_extent <- function(brush_data) {
      answer = NA
      
      if (exists('brush_data')) {
        answer = brush_data[['extent']]
      }
      
      return(answer)
    }
  )

  
  app$callback(
    output = list(id = 'ideo-homology-data', property = 'data'),
    params = list(
      input(id = 'chr-select-1', property = 'value'),
      input(id = 'chr-select-2',  property = 'value'),
      input(id = 'chrone-startone',  property = 'value'),
      input(id = 'chrone-stopone',  property = 'value'),
      input(id = 'chrone-starttwo',  property = 'value'),
      input(id = 'chrone-stoptwo',  property = 'value'),
      input(id = 'chrtwo-startone',  property = 'value'),
      input(id = 'chrtwo-stopone',  property = 'value'),
      input(id = 'chrtwo-starttwo',  property = 'value'),
      input(id = 'chrtwo-stoptwo',  property = 'value'),
      state(id = 'ideo-homology-data', property = 'data')
    ),
    
    change_homology_ideogram <- function(
      chr_selected_1,
      chr_selected_2,
      start_one,
      stop_one,
      start_two,
      stop_two,
      start_one_a,
      stop_one_a,
      start_two_a,
      stop_two_a,
      current
    ) {
      
      current = ideograms_initial[['homology']]
        
      current[['chromosomes']] = list(chr_selected_1, chr_selected_2)
      
      print(current)
      
      current[['homology']] = list(
        'chrOne' = list(
          'organism' = 9606,
          'start' = list(start_one, start_two),
          'stop' = list(stop_one, stop_two)
        ),
        
        'chrTwo' = list(
          'organism' = 9606,
          'start' = list(start_one_a, start_two_a),
          'stop' = list(stop_one_a, stop_two_a)
        )
      )
        
      return(current)
      
    }
  )
  
  
  
  
  
  
  # Custom Callback
  
  
  app$callback(
    output = list(id = 'ideo-custom-data', property = 'data'),
    params = list(
      input(id = 'organism-change', property = 'value'),
      input(id = 'bandlabel-switch', property = 'value'),
      input(id = 'chromlabel-switch', property = 'value'),
      input(id = 'orientation-switch', property = 'value'),
      input(id = 'chr-width-input', property = 'value'),
      input(id = 'chr-height-input', property = 'value'),
      input(id = 'chr-margin-input', property = 'value'),
      input(id = 'rotatable-switch', property = 'value'),
      input(id = 'resolution-select', property = 'value'),
      input(id = 'sex-switch', property = 'value'),
      input(id = 'fullband-switch', property = 'value'),
      state(id = 'ideo-custom-data', property = 'data')
    ),
    
    change_custom_ideogram <- function(
      organism_sel,
      show_band_labels,
      show_chromosome_labels,
      orientation_value,
      chr_width,
      chr_height,
      chr_margin,
      rotatable_value,
      resolution_value,
      sex_value,
      show_banded,
      current) {
      
  
      if (is.null('current')) {
        current = ideograms_initial[['custom']]
        
        current['organism'] = organism_sel
        current['showBandLabels'] = show_band_labels
        current['showChromosomeLabels'] = show_chromosome_labels
        current['orientation'] = orientation_value
        current['chrMargin'] = chr_margin
        current['chrWidth'] = chr_width
        current['chrHeight'] = chr_height
        current['rotatable'] = rotatable_value
        current['showFullyBanded'] = show_banded
        current['resolution'] = resolution_value
        current['sex'] = ifelse(sex_value, 'female', 'male')
      }
      
      else {
        current['organism'] = organism_sel
        current['showBandLabels'] = show_band_labels
        current['showChromosomeLabels'] = show_chromosome_labels
        current['orientation'] = orientation_value
        current['chrWidth'] = chr_width
        current['chrHeight'] = chr_height
        current['chrMargin'] = chr_margin
        current['rotatable'] = rotatable_value
        current['showFullyBanded'] = show_banded
        current['resolution'] = resolution_value
        current['sex'] = ifelse(sex_value, 'female', 'male')
      }
  
      return(current)
    }
  )
  


if (appName != "") {
  
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
  
} else {
  
  app$run_server(showcase = TRUE)
  
}

  
