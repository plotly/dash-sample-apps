appName <- Sys.getenv("DASH_APP_NAME")

if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

library(dashBio)
library(dashDaq)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(readr)
library(jsonlite)
library(RCurl)

app <- Dash$new()

dataset1 <- read_file("data/alignment_viewer_sample.fasta")
dataset2 <- read_file("data/alignment_viewer_p53.fasta")
dataset3 <- read_file("data/alignment_viewer_p53_clustalo.fasta")

DATASETS <- list("dataset1"=dataset1,
                 "dataset2"=dataset2,
                 "dataset3"=dataset3)

COLORSCALES_DICT <- list(
  list('value'= 'buried', 'label'= 'Buried'),
  list('value'= 'cinema', 'label'= 'Cinema'),
  list('value'= 'clustal2', 'label'= 'Clustal2'),
  list('value'= 'clustal', 'label'= 'Clustal'),
  list('value'= 'helix', 'label'= 'Helix'),
  list('value'= 'hydro', 'label'= 'Hydrophobicity'),
  list('value'= 'lesk', 'label'= 'Lesk'),
  list('value'= 'mae', 'label'= 'Mae'),
  list('value'= 'nucleotide', 'label'= 'Nucleotide'),
  list('value'= 'purine', 'label'= 'Purine'),
  list('value'= 'strand', 'label'= 'Strand'),
  list('value'= 'taylor', 'label'= 'Taylor'),
  list('value'= 'turn', 'label'= 'Turn'),
  list('value'= 'zappo', 'label'= 'Zappo')
)

GAP_COLORS_OPT = list(
  'black',
  'grey',
  'white',
  'turquoise',
  'blue',
  'green',
  'red',
  'purple'
)

CONSERVATION_COLORS_OPT = list(
  'Blackbody',
  'Bluered',
  'Blues',
  'Earth',
  'Electric',
  'Greens',
  'Greys',
  'Hot',
  'Jet',
  'Picnic',
  'Portland',
  'Rainbow',
  'RdBu',
  'Reds',
  'Viridis',
  'YlGnBu',
  'YlOrRd'
)

alignlayout <- htmlDiv(id="alignment-body", className="app-body",children=list(
  htmlDiv(list(
    htmlDiv(id='alignment-control-tabs', className='control-tabs', children=list(
      dccTabs(id='alignment-tabs', value='what-is',children=list(
        dccTab(label="About",
               value="what-is",
               children=htmlDiv(
                 className="control-tab", children=list(
                   htmlH4(className="what-is",
                          children="What is Alignment Viewer?"
                   ),
                   htmlP("The Alignment Viewer (MSA) component is used to align
                    multiple genomic or proteomic sequences from a FASTA or
                    Clustal file. Among its extensive set of features,
                    the multiple sequence alignment viewer can display
                    multiple subplots showing gap and conservation info,
                    alongside industry standard colorscale support and
                    consensus sequence. No matter what size your alignment
                    is, Alignment Viewer is able to display your genes or
                    proteins snappily thanks to the underlying WebGL
                    architecture powering the component. You can quickly
                    scroll through your long sequence with a slider or a
                    heatmap overview."
                   ),
                   htmlP("Note that the AlignmentChart only returns a chart of
                    the sequence, while AlignmentViewer has integrated
                    controls for colorscale, heatmaps, and subplots allowing
                    you to interactively control your sequences."
                   ),
                   htmlP("Read more about the component here:
                      https://github.com/plotly/react-alignment-viewer"
                   )
                 )
               )
        ),
        dccTab(
          label='Data',
          value='alignment-tab-select',
          children=htmlDiv(className='control-tab',children=list(
            htmlDiv(className="app-controls-block",children=list(
              htmlDiv(
                className='fullwidth-app-controls-name',
                children="Select preloaded dataset"
              ),
              dccDropdown(
                id="alignment-dropdown",
                options=list(
                  list(label= 'Sample.fasta', value= 'dataset1'),
                  list(label= 'P53.fasta naive', value= 'dataset2'),
                  list(label= 'P53.fasta aligned (ClustalW)', value= 'dataset3')
                ),
                value="dataset3"
              )
            )),
            htmlDiv(className="app-controls-block",children=list(
              htmlDiv(className="fullwidth-app-controls-name",
                      children="Upload your own dataset"),
              
              htmlA(
                htmlButton(
                  "Download sample data",
                  className="control-download"
                ),
                href="data/alignment_viewer_p53_clustalo.fasta",
                download="p53_clustalo.fasta"
              ),
              
              htmlDiv(id='alignment-file-upload-container', children=list(
                dccUpload(
                  id='alignment-file-upload',
                  className='control-upload',
                  children=htmlDiv(list(
                    "Drag and drop FASTA files or select files."
                  ))
                )
              ))
            ))
          ))
        ),
        dccTab(
          label='Interactions',
          value='control-tab-select2',
          children=htmlDiv(className='control-tab', children=list(
            htmlDiv(
              className='app-controls-name',
              children='Event Metadata'
            ),
            htmlP('Hover or click on data to see it here.'),
            htmlDiv(
              id='alignment-events'
            )
          )
          
          )
        ), #end Interactions tab
        dccTab(
          label="Graph",
          value="control-tab-customize",
          children=htmlDiv(
            className="control-tab",children=list(
              htmlDiv(list(
                htmlH3("General",className="alignment-settings-section"),
                htmlDiv(
                  className="app-controls-block",
                  children=list(
                    htmlDiv(
                      className='app-controls-name',
                      children="Colorscale",
                    ),
                    dccDropdown(
                      id='alignment-colorscale-dropdown',
                      className='app-controls-block-dropdown',
                      options=COLORSCALES_DICT,
                      value='clustal2'
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children='Choose the color theme of the viewer.'
                    )
                  )
                ),
                htmlDiv(
                  className='app-controls-block',
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='Overview'),
                    dccDropdown(
                      id='alignment-overview-dropdown',
                      className='app-controls-block-dropdown',
                      options=list(
                        list('label'='Heatmap', 'value'='heatmap'),
                        list('label'='Slider', 'value'='slider'),
                        list('label'='None', 'value'='none')
                      ),
                      value='heatmap'
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children='Show slider, heatmap or no overview.'
                    )
                  )
                ),
                htmlDiv(
                  className='app-controls-block',
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='Consensus'),
                    dccRadioItems(
                      id='alignment-showconsensus-radio',
                      className='alignment-radio',
                      options=list(
                        list("label"="Show","value"=TRUE),
                        list("label"="Hide","value"=FALSE)
                      ),
                      value=TRUE,
                      labelStyle=list(
                        "display"="inline-block",
                        "margin-right"="8px"
                      )
                    ), #end Consensus dccRadioItems
                    htmlDiv(
                      className='app-controls-desc',
                      children='Toggle the consensus (most frequent) sequence.'
                    )
                  )
                ),
                htmlDiv(
                  className="app-controls-block",
                  children=list(
                    htmlDiv(
                      className="app-controls-name",
                      children="Text size"
                    ),
                    dccSlider(
                      className='control-slider',
                      id='alignment-textsize-slider',
                      value=10,
                      min=8,
                      max=12,
                      marks= as.list(setNames(c(paste((8:12))), ((8:12)))),
                      step=1,
                      updatemode = "drag"
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children='Adjust the font size (in px) of viewer text.'
                    )
                  )
                )
              )),
              htmlHr(),
              htmlDiv(list(
                htmlH3("Conservation",className='alignment-settings-section'),
                htmlDiv(
                  className='app-controls-block',
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='Barplot'),
                    dccRadioItems(
                      id='alignment-showconservation-radio',
                      className='alignment-radio',
                      options=list(
                        list("label"="Show","value"=TRUE),
                        list("label"="Hide","value"=FALSE)
                      ),
                      value=TRUE,
                      labelStyle=list(
                        "display"="inline-block",
                        "margin-right"="8px"
                      )
                    ), #end Conservations dccRadioItems
                    htmlDiv(
                      className='app-controls-desc',
                      children='Show or hide the conservation barplot.'
                    )
                  )
                ),
                htmlDiv(
                  className='app-controls-block',
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='Colorscale'),
                    dccDropdown(
                      id='alignment-conservationcolorscale-dropdown',
                      className='app-controls-block-dropdown',
                      options= lapply(CONSERVATION_COLORS_OPT, function(x){
                        list(label = x, value = x)
                      }),
                      value="Viridis"
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children='Change the colorscale for the conservation barplot.'
                    )
                  )
                ),
                htmlDiv(
                  className='app-controls-block',
                  children=list(
                    htmlDiv(
                      className='app-controls-name',
                      children='Method'
                    ),
                    dccDropdown(
                      id='alignment-conservationmethod-dropdown',
                      className='app-controls-block-dropdown',
                      options=list(
                        list("label"="Entropy","value"="entropy"),
                        list("label"="Conservation","value"="conservation")
                      ),
                      value="entropy"
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children="Conservation (MLE) or normalized entropy."
                    )
                  )
                )
              )), # end Conservation block
              htmlHr(),
              htmlDiv(list(
                htmlH3("Conservation gap", 
                       className='alignment-settings-section'),
                htmlDiv(
                  className='app-controls-block',
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='Colorscale'),
                    dccRadioItems(
                      id='alignment-correctgap-radio',
                      className='alignment-radio',
                      options=list(
                        list('label'='Yes', 'value'=TRUE),
                        list('label'='No', 'value'=FALSE)
                      ),
                      value=TRUE,
                      labelStyle=list(
                        "display"="inline-block",
                        "margin-right"="8px"
                      )
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children="Lowers conservation of high gap sequences."
                    )
                  )
                ), #end app-controls block
                htmlDiv(
                  className="app-controls-block",
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='Gap'),
                    dccRadioItems(
                      id='alignment-showgap-radio',
                      className='alignment-radio',
                      options=list(
                        list('label'='Show', 'value'=TRUE),
                        list('label'='Hide', 'value'=FALSE)
                      ),
                      value=TRUE,
                      labelStyle=list(
                        "display"="inline-block",
                        "margin-right"="8px"
                      )
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children="Show/hide the gap barplot."
                    )
                  )
                ), #end app-controls-block
                htmlDiv(
                  className="app-controls-block",
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='Color'),
                    dccDropdown(
                      id='alignment-gapcolor-dropdown',
                      className='app-controls-block-dropdown',
                      options=lapply(GAP_COLORS_OPT, function(x){
                        list(label = x, value = x)
                      }),
                      value='grey',
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children="Set the color of the traces that represent the gap."
                    )
                  )
                ), #end Color app-controls-block
                htmlDiv(
                  className="app-controls-block",
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='Group'),
                    dccRadioItems(
                      id='alignment-groupbars-radio',
                      className='alignment-radio',
                      options=list(
                        list('label'='Yes', 'value'=TRUE),
                        list('label'='No', 'value'=FALSE)
                      ),
                      value=FALSE,
                      labelStyle=list(
                        "display"="inline-block",
                        "margin-right"="8px"
                      )
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children="Group gap and conservation bars."
                    )
                  )
                )
              )),#end Conservation Gap Block
              htmlHr(),
              htmlDiv(list(
                htmlH3("Layout", className='alignment-settings-section'),
                htmlDiv(
                  className='app-controls-block',
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='Labels'),
                    dccRadioItems(
                      id='alignment-showlabel-radio',
                      className='alignment-radio',
                      options=list(
                        list('label'='Show', 'value'=TRUE),
                        list('label'='Hide', 'value'=FALSE)
                      ),
                      value=TRUE,
                      labelStyle=list(
                        "display"="inline-block",
                        "margin-right"="8px"
                      )
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children="Show track labels on the left."
                    )
                  )
                ),
                htmlDiv(
                  className='alignment-settings',
                  children=list(
                    htmlDiv(className='app-controls-name',
                            children='IDs'),
                    dccRadioItems(
                      id='alignment-showid-radio',
                      className='alignment-radio',
                      options=list(
                        list('label'='Show', 'value'=TRUE),
                        list('label'='Hide', 'value'=FALSE)
                      ),
                      value=TRUE,
                      labelStyle=list(
                        "display"="inline-block",
                        "margin-right"="8px"
                      )
                    ),
                    htmlDiv(
                      className='app-controls-desc',
                      children='Show track IDs on the left.'
                    )
                  )
                )
              )) #end Layout Block
              
            ))
        ) #end Graph Tab
      ))
    ))
  )),
  dccLoading(
    className='dashbio-loading',
    children=htmlDiv(list(
      dashbioAlignmentChart(
        id='alignment-chart',
        height=725,
        data=dataset3,
      )
    ))
  ), #end dccLoading
  dccStore(id='alignment-data-store')
))

###HEADER###
header <- htmlDiv(children = list(
  htmlSpan(list(htmlA(
    children = list(
      htmlImg(src='assets/plotly-dash-bio-logo.png', height = '36', width = '180',
              style = list('top' = '10', 'margin-left' = '10px', 'margin-top' = '20px', 'margin-right' = '20px'))
    ),
    href='https://dash-bio.plotly.host/Portal/'
  ))),
  "Alignment Viewer",
  
  htmlA(href="https://github.com/plotly/dash-sample-apps/blob/master/apps/dashr-alignment-viewer/app.R",
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
  
  htmlImg(src = "assets/github.png", 
          width = '50', 
          height = '50', 
          style = list(
            "position" = "relative", "float" = "right",
            "margin-right" = "10px",
            "top" = "17px"
          )
  )
), className = 'alignmentviewerheader')

app$layout(htmlDiv(
  list(
    header,
    alignlayout,
    htmlDiv(list(" "), className = "blankspace"),
    daqToggleSwitch(
      id = 'daq-light-dark-theme',
      color = '#FDFCFF',
      labelPosition = 'bottom',
      label = list('Light', 'Dark'),
      style = list(
        'width' = '250px',
        'margin' = 'auto',
        'display' = 'none'
      ),
      value = FALSE
    )
  ),
  style = list('background-color' = '#FDFCFF')
))

app$callback(
  output('alignment-data-store', 'data'),
  list(
    input('alignment-dropdown', 'value'),
    input('alignment-file-upload', 'contents'),
    input('alignment-file-upload', 'filename')
  ),

  update_storage <- function(dropdown, contents, filename) {

    if ((!is.null(contents)) && grepl("fasta",filename)) {
      content_type <- strsplit(contents, ",")[[1]][1]
      content_string <- strsplit(contents, ",")[[1]][2]
      content <- base64Decode(content_string)
    } else {
      content <- DATASETS[[dropdown]]
    }
    return(content)
  }
)

# Handle event data
app$callback(
  output("alignment-events", "children"),
  list(input("alignment-chart", "eventDatum")),

  event_data_select <- function(data) {
    if (!exists("data")) {
      return("No data.")
    } else {
      dj <- fromJSON(data)
    }

    if(length(names(dj))==0){
      return('No event data to display.')
    } else {
      l<- paste0("- ",names(dj),": ",dj)
      q <- list()
      for(i in 1:length(l)){
        q[[i]]<- htmlDiv(l[i])
      }
      return(q)
    }

  }
)

# Render main chart
app$callback(
  output("alignment-chart",'data'),
  list(input('alignment-data-store',"data")),

  update_chart <- function(input_data){
    return(input_data)
  }
)

#Customization callbacks
app$callback(
  output('alignment-chart', 'colorscale'),
  list(input('alignment-colorscale-dropdown', 'value')),

  customize_colorscale <- function(val) {
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'overview'),
  list(input('alignment-overview-dropdown', 'value')),

  customize_overview <- function(val) {
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'showconsensus'),
  list(input('alignment-showconsensus-radio', 'value')),

  customize_showconsensus <- function(val){
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'textsize'),
  list(input('alignment-textsize-slider', 'value')),

  customize_textsize <- function(val){
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'showconservation'),
  list(input('alignment-showconservation-radio', 'value')),

  customize_showconservation <- function(val){
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'conservationcolorscale'),
  list(input('alignment-conservationcolorscale-dropdown', 'value')),

  customize_conservationcolorscale <- function(val){
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'conservationmethod'),
  list(input('alignment-conservationmethod-dropdown', 'value')),

  customize_conservationmethod <- function(val){
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'correctgap'),
  list(input('alignment-correctgap-radio', 'value')),

  customize_correctgap <- function(val){
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'showgap'),
  list(input('alignment-showgap-radio', 'value')),

  customize_showgap <- function(val){
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'gapcolor'),
  list(input('alignment-gapcolor-dropdown', 'value')),

  customize_gapcolor <- function(val){
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'groupbars'),
  list(input('alignment-groupbars-radio', 'value')),

  customize_groupbars <- function(val)
  return(val)

)

app$callback(
  output('alignment-chart', 'showlabel'),
  list(input('alignment-showlabel-radio', 'value')),

  customize_showlabe <- function(val){
    return(val)
  }
)

app$callback(
  output('alignment-chart', 'showid'),
  list(input('alignment-showid-radio', 'value')),

  customize_showid <- function(val){
    return(val)
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(showcase = TRUE)
}
