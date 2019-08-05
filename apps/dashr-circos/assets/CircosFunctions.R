# Dataframe for Table

df <- circos_dataframe[['chords']]

df <- flatten(df, recursive = TRUE)



description <- function() {
  return("'Vizualize and analyze similarities and differences between ' \
'genes in a single plot, using the powerful Circos graph.'")
}

# Header colors
header_color <- function() {
  return(list('bg_color' = '#262B3D', 'font_color' = '#FFF', 'light_logo' = TRUE))
}


# Empty Circos needed for circos graph callback

empty <- dashbioCircos(
  id = 'main-circos',
  selectEvent = list(),
  layout = list(),
  size = 700,
  config = list(),
  tracks = list(),
  enableZoomPan = TRUE,
  enableDownloadSVG = TRUE
)


circos_layout <- function() {
  return(
    htmlDiv(
      id='circos-body',
      className = 'app-body',
      children= list(
        htmlDiv(id = 'graph-container', className = 'circos-hold',children = list(
        dccLoading(className = 'dashbio-loading', children = htmlDiv(
          id='circos-hold',
          className = 'circos-hold',
          children = list(),
          style = list()
        )))),
        dccLoading(className = 'dashbio-loading', children = htmlDiv(
          id='chords-plot',
          className = 'circos-hold',
          children = list(),
          style = list()
        )),
        
      
      htmlDiv(id = 'circos-control-tabs', className = 'control-tabs', children=list(
        dccTabs(id='circos-tabs', value = 'what-is', children = list(
          dccTab(
            label = 'About',
            value = 'what-is',
            children = htmlDiv(
              id = 'control-tab', children = list(
                htmlH4(className = 'what-is', children = "What is Circos?"),
                
                  htmlP("Circos is a circular visualization of data, and can be 
                  used to highlight relationships between objects in a 
                  dataset (e.g., genes that are located on different chromosomes
                  in the genome of an organism."),
                  
                  htmlP("A Dash Circos graph consists of two main parts: the layout and
                  the tracks. The layout sets the basic parameters of the graph, 
                  such as radius, ticks, labels, etc; the tracks are graph layouts
                  that take in a series of data points to display."),
                  
                  htmlP("The visualizations supported by Dash Circos are: heatmaps, chords,
                  highlights, histograms, line, scatter, stack, and text graphs."),
                  
                  htmlP("In the 'Data' tab, you can opt to use preloaded datasets; additionally
                  you can download sample data you would use with a Dash Circos component,
                  upload that sample data, and render it with the 'Render' button."), 
                  
                  htmlP("In the 'Graph' tab, you can choose the type of Circos graph to 
                  display, control the size of the graph, and access data that are
                  generated upon hovering over parts of the graph."),
                  
                  htmlP("In the 'Table' tab, you can view the datasets that define the parameters
                  of the graph, such as the layout, the highlights, and the chords. 
                  You can interact with Circos through this table by selecting the 
                  'Chords' graph in the 'Graph' tab, then viewing the 'Chords' dataset
                  in the 'Table' tab."),
                
                
                htmlDiv(list(
                  'Reference:  ',
                  htmlA('Seminal paper', href='http://www.doi.org/10.1101/gr.092759.109')
                )),
                
                htmlBr(),
                
                htmlDiv(list(
                  'For a look into Circos and the Circos API, please visit the 
              original repository ',
                  htmlA('here', href='https://github.com/nicgirault/circosJS'),
                  '.'
                ))
              )
            )
          ),
          
          dccTab(
            label='Data',
            value = 'data',
            children = htmlDiv(className = 'control-tab', children = list(
              htmlDiv(className = 'app-controls-block', children = list(
                htmlDiv(className = 'app-controls-name', children = 'Data source'),
                dccDropdown(
                  id = 'circos-preloaded-uploaded',
                  options = list(
                    list('label' = 'Preloaded', 'value'='preloaded'),
                    list('label' ='Upload', 'value' = 'upload')
                    ), value = 'preloaded'
                  )
                )
              ),
              htmlHr(),
              htmlA(
                htmlButton(
                  id='circos-download-button',
                  className = 'control-download',
                  children = 'Download Sample Data'
                ),
                href ="assets/circos_graph_data.json",
                download = 'circos_graph_data.json'
              ),
              
              htmlDiv(id = 'circos-uploaded-data', children = list(htmlDiv(id = 'options', style = list(), list(
                dccUpload(
                  id='upload-data',
                  className = 'control-upload',
                  children = htmlDiv(list(
                    "Drag and Drop of click to import .csv file here!"
                  )) , multiple = TRUE
                ),
                
                htmlDiv(className = 'app-controls-block', children = list(
                  htmlDiv(className = 'app-controls-name', children = 'Select Upload Data'),
                  dccDropdown(
                    id = "circos-view-dataset-custom",
                    options = list(
                      list("label"="Layout", "value" = 0),
                      list("label" = "Track 1", "value" = 1),
                      list("label" = "Track 2", "value" = 2)
                    ),
                    value = 0
                  )
                )),
                htmlButton(
                  "Render uploaded dataset",
                  id = "render-button",
                  className = "control-download"
                )
              ))
              
            ))))
          ),
          
          dccTab(
            label = 'Graph',
            value = 'graph',
            children = htmlDiv(className = 'control-tab', children= list(
              htmlDiv(className = 'app-controls-block', children=list(
                htmlDiv(className = 'app-controls-name', children = 'Graph Type'),
                dccDropdown(
                  id='circos-graph-type',
                  options = list(
                    list('label' = 'Heatmap', 'value' = 'heatmap'),
                    list('label' = 'Histogram', 'value' = 'histogram'), 
                    list('label' = 'Chords', 'value' = 'chords'),
                    list('label' = 'Line', 'value' = 'line'),
                    list('label' = 'Highlight', 'value' = 'highlight'),
                    list('label' = 'Text', 'value' = 'text'),
                    list('label' = 'Stack', 'value' = 'stack'),
                    list('label' = 'Scatter', 'value' = 'scatter')
                    
                  ),
                  value = 'chords'
                ),
                htmlDiv(className = 'app-controls-desc', id = 'chords-text')
              )),
              htmlDiv(className = 'app-controls-block', children = list(
                htmlDiv(className = 'app-controls-name', children = 'Graph size'),
                dccSlider(
                  id='circos-size',
                  min = 500,
                  max = 800,
                  step = 10,
                  value = 650
                )
              )),
              htmlHr(),
              htmlH5('Hover Data'),
              htmlDiv(
                id = 'event-data-select',
                children = list()
              )
            ))
          ),
          
          dccTab(
            label = 'Table',
            value = 'table',
            children = htmlDiv(className = 'control-tab', children = list(
              htmlDiv(className = 'app-controls-name', children = list(
                htmlDiv(className = 'app-controls-name', children = 'View dataset'),
                dccDropdown(
                  id = 'circos-view-dataset',
                  options = list(
                    list('label' = 'Cytobands', 'value'= 'cytobands'),
                    list('label' = 'Chords', 'value'= 'chords'),
                    list('label' = 'Highlight' , 'value' = 'highlights')
                    
                  ), value = 'chords'
                )
              )),
              htmlDiv(id = 'circos-table-container',
                      children = list(
                        dashDataTable(
                          id = 'data-table',
                          columns = lapply(colnames(df), 
                                           function(colName){
                                             list(
                                               id = colName,
                                               name = colName
                                             )
                                           }),
                          data= df_to_list(df),
                          row_selectable='multi',
                          selected_rows = character(0),
                          # sorting = TRUE,
                          # filtering = TRUE,
                          # css = list(list(
                          #   'selector' = '.dash-cell div.dash-cell-value',
                          #   'rule' = 'display: inline;',
                          #            'white-space: inherit;',
                          #            'overflow: auto;',
                          #            'text-overflow: inherit;'
                          # )),
                          style_cell = list(
                            'whiteSpace' = 'no-wrap',
                            'overflow' = 'hidden',
                            'textOverflow' = 'ellipsis',
                            'maxWidth' = 100,
                            'fontWeight' = 100,
                            'fontSize' = '11pt',
                            'fontFamily' = 'Courier New',
                            'backgroundColor' = '#1F2132'
                          ),
                          
                          style_header = list(
                            'backgroundColor' = '#1F2132',
                            'textAlign' = 'center'
                          ),
                          
                          style_table = list(
                            'maxHeight' = '310px',
                            'width' = '320px',
                            'marginTop' = '5px',
                            'marginBottom' = '10px'
                          )
                        )
                      )),
              
              htmlDiv(id = 'expected-index')
            ))
          )
        ))
      )
        ) ,
      
      htmlDiv(list(
        htmlDiv(id='output-data-upload'),
        htmlDiv(id='event-data-store')
      ), className = 'circos-display-none')
        
        
      )
    )
  )
}




get_circos_graph <- function(key, size) {
  # if (is.null(data) == TRUE) {
  #   data = list(NULL,NULL,NULL)
  # }
  
  
  circos_graphs = list(
    # 'upload-custom-dataset' = dashbioCircos(
    #   id = 'main-circos',
    #   selectEvent = list('0' = 'both', '1' = 'both'),
    #   layout = data[[1]],
    #   config = list(
    #     'innerRadius' = size/2 - 80,
    #     'outerRadius' = size/2 - 40,
    #     'ticks' = list('display' = FALSE, 'labelDenominator' = 1000000),
    #     'labels' = list(
    #       'position' = 'center',
    #       'display' = FALSE,
    #       'size' = 12,
    #       'color' = '#fff',
    #       'radialOffset' = 70
    #     )
    #   ),
    # 
    #   tracks = list(
    #     list(
    #       'type' = 'HIGHLIGHT',
    #       'data' = data[[2]],
    #       'config' = list(
    #         'innerRadius' = size/2 -80,
    #         'outerRadius' = size/2 -40,
    #         'opacity' = 0.3,
    #         'tooltipContent' = list('name' = 'all'),
    #         'color' = list('name' = 'color')
    #       )
    #     ),
    #     list(
    #       'type' = 'HIGHLIGHT',
    #       'data' = data[[3]],
    #       'config' = list(
    #         'innerRadius' = size/2 -80,
    #         'outerRadius' = size/2 -40,
    #         'opacity' = 0.3,
    #         'tooltipContent' = list('name' = 'all'),
    #         'color' = list('name' = 'color')
    #       )
    #     )
    #   ), size = 700
    # ),
    # 
    # 
    # 'select-dataset-parser' = dashbioCircos(
    #   id = 'main-circos',
    #   selectEvent = list('0' = 'hover', '1' = 'click'),
    #   layout = parsed_layout,
    #   config = list(
    #     'innerRadius' = size/2 - 80,
    #     'outerRadius' = size/2 - 40,
    #     'ticks' = list('display' = FALSE, 'labelDenominator' = 1000000),
    #     'labels' = list(
    #       'position' = 'center',
    #       'display' = FALSE,
    #       'size' = 8,
    #       'color' = '#fff',
    #       'radialOffset' = 90
    #     )
    #   ),
    # 
    #   tracks = list(
    #     list(
    #       'type' = 'HIGHLIGHT',
    #       'data' = parsed_track_one,
    #       'config' = list(
    #         'innerRadius' = size/2 -80,
    #         'outerRadius' = size/2 -40,
    #         'opacity' = 0.3,
    #         'tooltipContent' = list('name' = 'block_id'),
    #         'color' = list('name' = 'color')
    #       )
    #     ),
    #     list(
    #       'type' = 'HIGHLIGHT',
    #       'data' = parsed_track_two,
    #       'config' = list(
    #         'innerRadius' = size/2 -80,
    #         'outerRadius' = size/2 -40,
    #         'opacity' = 0.3,
    #         'tooltipContent' = list('name' = 'block_id'),
    #         'color' = list('name' = 'color')
    #       )
    #     )
    #   ), size = 700
    # ),
    
    
    'heatmap' = dashbioCircos(
      id = 'main-circos',
      selectEvent = list('0' = 'hover', '1' = 'hover'),
      layout = circos_graph_data[['month_layout']],
      config = list(
        'innerRadius' = size/2 - 80,
        'outerRadius' = size/2 - 30,
        'ticks' = list('display' = FALSE),
        'labels' = list(
          'position' = 'center',
          'display' = TRUE,
          'size' = 14,
          'color' = '#fff',
          'radialOffset' = 15
        )
      ),
      
      tracks = list(
        list(
          'type' = 'HEATMAP',
          'data' = circos_graph_data[['heatmap']],
          'config' = list(
            'innerRadius' = 0.8,
            'outerRadius' = 0.98,
            'logScale' = FALSE,
            'tooltipContent' = list('name' = 'value'),
            'color' = 'YlOrRd'
          )
        ),
        list(
          'type' = 'HEATMAP',
          'data' = circos_graph_data[['heatmap']],
          'config' = list(
            'innerRadius' = 0.7,
            'outerRadius' = 0.79,
            'logScale' = FALSE,
            'tooltipContent' = list('name' = 'value'),
            'color' = 'Blues'
          )
        )
      ), size = 700
    ),
    
    
    'chords' = dashbioCircos(
      id = 'main-circos',
      selectEvent = list('0' = 'both', '1' = 'both'),
      layout = circos_graph_data[['GRCh37']],
      config = list(
        'innerRadius' = size/2 - 80,
        'outerRadius' = size/2 - 30,
        'ticks' = list('display' = FALSE, 'labelDenominator' = 1000000),
        'labels' = list(
          'position' = 'center',
          'display' = TRUE,
          'size' = 11,
          'color' = '#fff',
          'radialOffset' = 75
        )
      ),
      
      tracks = list(
        list(
          'type' = 'HIGHLIGHT',
          'data' = circos_graph_data[['cytobands']],
          'config' = list(
            'innerRadius' = size/2 - 80,
            'outerRadius' = size/2 - 40,
            'opacity' = 0.3,
            'tooltipContent' = list('name' = 'all'),
            'color' = list('name' = 'color')
          )
        ),
        list(
          'type' = 'CHORDS',
          'data' = circos_graph_data[['chords']],
          'config' = list(
            'opacity' = 0.7,
            'color' = list('name' = 'color'),
            'tooltipContent' = list(
              'source' = 'source',
              'sourceID' = 'id',
              'target' = 'target',
              'targetID' = 'id',
              'targetEnd' = 'end'
            )
          )
        )
      ), size = 700
    ),
    
    
    
    'highlight' = dashbioCircos(
      id = 'main-circos',
      selectEvent = list('0' = 'hover'),
      layout = circos_graph_data[['GRCh37']],
      config = list(
        'innerRadius' = size/2 - 100,
        'outerRadius' = size/2 - 50,
        'ticks' = list('display' = FALSE),
        'labels' = list('display' = FALSE)
      ),
      
      tracks = list(
        list(
          'type' = 'HIGHLIGHT',
          'data' = circos_graph_data[['cytobands']],
          'config' = list(
            'innerRadius' = size/2 - 100,
            'outerRadius' = size/2 - 50,
            'opacity' = 0.5,
            'tooltipContent' = list('name' = 'name'),
            'color' = list('name' = 'color')
          )
        )
      ), size = 700
    ),
    
    
    'histogram' = dashbioCircos(
      id = 'main-circos',
      selectEvent = list('0' = 'hover', '1' = 'hover'),
      layout = circos_graph_data[['GRCh37']],
      config = list(
        'innerRadius' = size/2 - 150,
        'outerRadius' = size/2 - 120,
        'ticks' = list('display' = FALSE, 'labelDenominator' = 1000000),
        'labels' = list('display' = FALSE)
      ),
      
      tracks = list(
        list(
          'type' = 'HIGHLIGHT',
          'data' = circos_graph_data[['cytobands']],
          'config' = list(
            'innerRadius' = size/2 - 150,
            'outerRadius' = size/2 - 120,
            'opacity' = 0.6,
            'tooltipContent' = list('name' = 'name'),
            'color' = list('name' = 'color')
          )
        ),
        list(
          'type' = 'HISTOGRAM',
          'data' = circos_graph_data[['histogram']],
          'config' = list(
            'innerRadius' = 1.01,
            'outerRadius' = 1.4,
            'color' = list('name' = 'color'),
            'tooltipContent' = list('name' = 'value')
          )
        )
      ), size = 700
    ),
    
    
    'line' = dashbioCircos(
      id = 'main-circos',
      selectEvent = list(
        '0' = 'both',
        '1' = 'both',
        '2' = 'both',
        '3' = 'both',
        '4' = 'both',
        '5' = 'both',
        '6' = 'both',
        '7' = 'both'
      ),
      layout = Filter(function(x){x$id %in% c("chr1", "chr2", "chr3")}, circos_graph_data$GRCh37),
      config = list(
        'innerRadius' = size/2 - 150,
        'outerRadius' = size/2 - 130,
        'ticks' = list('display' = FALSE, 'spacing' = 1000000, 'labelSuffix' = ''),
        'labels' = list('display' = FALSE,
                        'position' = 'center',
                        'size' = 14,
                        'color' = '#fff',
                        'radialOffset' = 30
        )
      ),
      
      tracks = list(
        list(
          'type' = 'HIGHLIGHT',
          'data' = Filter(function(x){x$id %in% c("chr1", "chr2", "chr3")}, circos_graph_data$cytobands),
          'config' = list(
            'innerRadius' = size/2 - 150,
            'outerRadius' = size/2 - 130,
            'opacity' = 0.3,
            'tooltipContent' = list('name' = 'name'),
            'color' = list('name' = 'color')
          )
        ),
        list(
          'type' = 'LINE',
          'data' = circos_graph_data[['snp250']],
          'config' = list(
            'innerRadius' = 0.5,
            'outerRadius' = 0.8,
            'color' = '#222222',
            'tooltipContent' = list('source' = 'block_id',
                                    'target' = 'position',
                                    'targetEnd' = 'value'
            ),
            'axes' = list(
              list(
                'spacing' = 0.001,
                'thickness' = 1,
                'color' = '#666666'
              )
            ),
            
            'backgrounds' = list(
              list(
                'start' = 0,
                'end' = 0.002,
                'color' = '#f44336',
                'opacity' = 0.5
              ),
              list(
                'start' = 0.006,
                'end' = 0.002,
                'color' = '#4caf50',
                'opacity' = 0.5
              )
            )
            
          )
        ),
        
        list(
          'type' = 'SCATTER',
          'data' = circos_graph_data[['snp250']],
          'config' = list(
            'innerRadius' = 0.5,
            'outerRadius' = 0.8,
            'min' = 0,
            'max' = 0.015,
            'fill' = FALSE,
            'strokeWidth' = 0,
            'tooltipContent' = list(
              'source' = 'block_id',
              'target' = 'position',
              'targetEnd' = 'value'
            )
          )
        ),
        
        list(
          'type' = 'LINE',
          'data' = circos_graph_data[['snp']],
          'config' = list(
            'innerRadius' = 1.01,
            'outerRadius' = 1.15,
            'maxGap' = 1000000,
            'min' = 0,
            'max' = 0.015,
            'color' = '#222222',
            'tooltipContent' = list('name' = 'value'),
            'axes' = list(
              list('position' = 0.002, 'color' = '#f44336'),
              list('position' = 0.006, 'color' = '#4caf50')
            )
          )
        ),
        
        
        list(
          'type' = 'LINE',
          'data' = circos_graph_data[['snp1m']],
          'config' = list(
            'innerRadius' = 1.01,
            'outerRadius' = 1.15,
            'maxGap' = 1000000,
            'min' = 0,
            'max' = 0.015,
            'color' = '#f44336',
            'tooltipContent' = list('name' = 'value')
          )
        ),
        
        
        list(
          'type' = 'LINE',
          'data' = circos_graph_data[['snp']],
          'config' = list(
            'innerRadius' = 0.85,
            'outerRadius' = 0.95,
            'maxGap' = 1000000,
            'direction' = 'in',
            'min' = 0,
            'max' = 0.015,
            'color' = '#222222',
            'axes' = list(
              list('position' = 0.01, 'color' = '#4caf50'),
              list('position' = 0.008, 'color' = '#4caf50'),
              list('position' = 0.006, 'color' = '#4caf50'),
              list('position' = 0.002, 'color' = '#4caf50')
              
            )
          )
        ),
        
        list(
          'type' = 'LINE',
          'data' = circos_graph_data[['snp1m']],
          'config' = list(
            'innerRadius' = 0.85,
            'outerRadius' = 0.95,
            'maxGap' = 1000000,
            'direction' = 'in',
            'min' = 0,
            'max' = 0.015,
            'color' = '#f44336',
            'tooltipContent' = list('name' = 'value')
          )
        )
      ), size = 700
    ),
    
    'scatter' = dashbioCircos(
      id = 'main-circos',
      selectEvent = list(
        '0' = 'both',
        '1' = 'both',
        '2' = 'both',
        '3' = 'both',
        '4' = 'both'
      ),
      
      layout = Filter(function(x){x$id %in% c("chr1", "chr2", "chr3")}, circos_graph_data$GRCh37),
      config = list(
        'innerRadius' = size/2 - 150,
        'outerRadius' = size/2 - 130,
        'ticks' = list('display' = FALSE, 'spacing' = 1000000, 'labelSuffix' = ''),
        'labels' = list('display' = FALSE)
      ),
      
      tracks = list(
        list(
          'type' = 'HIGHLIGHT',
          'data' = Filter(function(x){x$id %in% c("chr1", "chr2", "chr3")}, circos_graph_data$cytobands),
          'config' = list(
            'innerRadius' = size/2 - 150,
            'outerRadius' = size/2 - 130,
            'opacity' = 0.8,
            'tooltipContent' = list('name' = 'name'),
            'color' = list('name' = 'color')
          )
        ),
        list(
          'type' = 'SCATTER',
          'data' = Filter(function(x){x$value > 0.007}, circos_graph_data$snp250),
          'config' = list(
            'innerRadius' = 0.65,
            'outerRadius' = 0.95,
            'color' = list('colorData' = 'name'),
            'tooltipContent' = list('source' = 'block_id',
                                    'target' = 'position',
                                    'targetEnd' = 'value'
            ),
            'strokeColor' = 'grey',
            'strokeWidth' = 1,
            'shape' = 'circle',
            'size' = 14,
            'min' = 0,
            'max' = 0.013,
            
            'axes' = list(
              list(
                'spacing' = 0.001,
                'thickness' = 1,
                'start' = 0.006,
                'color' = '#4caf50',
                'opacity' = 0.3
              ),
              list(
                'spacing' = 0.002,
                'thickness' = 1,
                'start' = 0.006,
                'color' = '#4caf50',
                'opacity' = 0.5
              ),
              list(
                'spacing' = 0.002,
                'thickness' = 1,
                'start' = 0.002,
                'end' = 0.006,
                'color' = '#666',
                'opacity' = 0.5
              ),
              list(
                'spacing' = 0.001,
                'thickness' = 1,
                'end' = 0.002,
                'color' = '#f44336',
                'opacity' = 0.5
              )
            ),
            
            'backgrounds' = list(
              list(
                'start' = 0.006,
                'color' = '#4caf50',
                'opacity' = 0.1
              ),
              list(
                'start' = 0.002,
                'end' = 0.006,
                'color' = '#d3d3d3',
                'opacity' = 0.1
              ),
              list(
                'end' = 0.002,
                'color' = '#f44336',
                'opacity' = 0.1
              )
            )
            
          )
        ),
        
        list(
          'type' = 'SCATTER',
          'data' = circos_graph_data[['snp250']],
          'config' = list(
            'innerRadius' = 1.075,
            'outerRadius' = 1.175,
            'min' = 0.007,
            'max' = 0.013,
            'color' = '#4caf50',
            'strokeWidth' = 1,
            'shape' = 'rectangle',
            'size' = 10,
            'strokeColor' = 'green',
            'tooltipContent' = list(
              'source' = 'block_id',
              'target' = 'position',
              'targetEnd' = 'value'
            ),
            'axes' = list(
              list(
                'spacing' = 0.001,
                'thickness' = 1,
                'color' = '#4caf50',
                'opacity' = 0.3
              ),
              list(
                'spacing' = 0.002,
                'thickness' = 1,
                'color' = '#4caf50',
                'opacity' = 0.5
              )
            ),
            
            'backgrounds' = list(
              list('start' = 0.007, 'color' = '#4caf50', 'opacity' = 0.1),
              list('start' = 0.009, 'color' = '#4caf50', 'opacity' = 0.1),
              list('start' = 0.011, 'color' = '#4caf50', 'opacity' = 0.1),
              list('start' = 0.013, 'color' = '#4caf50', 'opacity' = 0.1)
            )
          )
        ),
        
        
        list(
          'type' = 'SCATTER',
          'data' = Filter(function(x){x$value < 0.002}, circos_graph_data$snp250),
          'config' = list(
            'innerRadius' = 0.35,
            'outerRadius' = 0.60,
            'min' = 0,
            'max' = 0.002,
            'color' = '#f44336',
            'strokeWidth' = 1,
            'shape' = 'triangle',
            'size' = 10,
            'strokeColor' = 'red',
            'tooltipContent' = list(
              'source' = 'block_id',
              'target' = 'position',
              'targetEnd' = 'value'
            ),
            'axes' = list(
              list(
                'spacing' = 0.001,
                'thickness' = 1,
                'color' = '#f44336',
                'opacity' = 0.3
              ),
              list(
                'spacing' = 0.002,
                'thickness' = 1,
                'color' = '#f44336',
                'opacity' = 0.5
              )
            ),
            
            'backgrounds' = list(
              list('start' = 0.0004, 'color' = '#f44336', 'opacity' = 0.1),
              list('start' = 0.0008, 'color' = '#f44336', 'opacity' = 0.1),
              list('start' = 0.0012, 'color' = '#f44336', 'opacity' = 0.1),
              list('start' = 0.0016, 'color' = '#f44336', 'opacity' = 0.1),
              list('start' = 0.002, 'color' = '#f44336', 'opacity' = 0.1)
            )
          )
        ),
        
        
        list(
          'type' = 'SCATTER',
          'data' = circos_graph_data[['snp250']],
          'config' = list(
            'innerRadius' = 0.65,
            'outerRadius' = 0.95,
            'min' = 0,
            'max' = 0.013,
            'strokeWidth' = 1,
            'shape' = 'circle',
            'size' = 14,
            'strokeColor' = 'grey',
            'tooltipContent' = list(
              'source' = 'block_id',
              'target' = 'position',
              'targetEnd' = 'value'
            ),
            'axes' = list(
              list(
                'spacing' = 0.001,
                'start' = 0.006,
                'thickness' = 1,
                'color' = '#4caf50',
                'opacity' = 0.3
              ),
              list(
                'spacing' = 0.002,
                'start' = 0.006,
                'thickness' = 1,
                'color' = '#4caf50',
                'opacity' = 0.5
              ),
              list(
                'spacing' = 0.002,
                'start' = 0.002,
                'end' = 0.006,
                'thickness' = 1,
                'color' = '#666',
                'opacity' = 0.5
              ),
              list(
                'spacing' = 0.001,
                'end' = 0.002,
                'thickness' = 1,
                'color' = '#f44336',
                'opacity' = 0.5
              )
            ),
            'backgrounds' = list(
              list('start' = 0.006, 'color' = '#4caf50', 'opacity' = 0.1),
              list('start' = 0.002, 'end' = 0.006,'color' = '#d3d3d3', 'opacity' = 0.1),
              list('start' = 0.002, 'color' = '#f44336', 'opacity' = 0.1)
            )
          )
        )
        
        
      ), size = 700
    ),
    
    
    'stack' = dashbioCircos(
      id = 'main-circos',
      selectEvent = list('0' = 'hover'),
      layout = list(
        list(
          'id'= 'chr9',
          'len' = 8000000,
          'label' = 'chr9',
          'color' = '#FFCC00'
        )
      ),
      config = list(
        'innerRadius' = size/2 - 80,
        'outerRadius' = size/2 - 30,
        'ticks' = list('display' = FALSE, 'labels' = FALSE, 'spacing' = 10000),
        'labels' = list('display' = FALSE, 'labels' = FALSE, 'spacing' = 10000)
      ),
      
      tracks = list(
        list(
          'type' = 'STACK',
          'data' = circos_graph_data[['stack']],
          'config' = list(
            'innerRadius' = 0.7,
            'outerRadius' = 1,
            'thickness' = 4,
            'margin' = 0.01*8000000,
            'direction' = 'out',
            'strokeWidth' = 0,
            'opacity' = 0.5,
            'tooltipContent' = list('name' = 'chr'),
            'color' = list(
              'conditional' = list(
                'end' = 'end',
                'start' = 'start',
                'value' = list(150000, 120000, 90000, 60000, 30000),
                'color' = list(
                  'red',
                  'black',
                  '#fff',
                  '#999',
                  '#BBB'
                )
              )
            )
          )
        )
      ), size = 700
    ),
    
    
    'text' = dashbioCircos(
      id = 'main-circos',
      selectEvent = list('0' = 'hover', '1' = 'both'),
      layout = circos_graph_data[['GRCh37']][1],
      config = list(
        list(
          'innerRadius' = size/2 - 100,
          'outerRadius' = size/2 - 80,
          'ticks' = list('display' = FALSE),
          'labels' = list('display' = FALSE)
        )
      ),
      
      tracks = list(
        list(
          'type' = 'HIGHLIGHT',
          'data' = Filter(function(x){x$block_id == circos_graph_data$GRCh37[1]['id']}, circos_graph_data$cytobands),
          'config' = list(
            'innerRadius' = size / 2 - 100,
            'outerRadius' = size / 2 - 80,
            'opacity' = 0.7,
            'tooltipContent' = list('name' = 'name'),
            'color' = list('name' = 'color')
          )
        )
      ), size = 700
    )
    
    
    
  )
  
  return(circos_graphs[[key]])
}









