appName <- Sys.getenv("DASH_APP_NAME")


if (appName != ""){
  
  pathPrefix <- sprintf("/%s/", appName)
  
  
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  
  
  setwd(sprintf("/app/apps/%s", appName))
  
}

library(dashBio)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(readr)
library(jsonlite)
library(Biostrings)
library(RCurl)
library(sets)
source("protein_reader.R")

app <- Dash$new()


initialCov <- list(
  list('start'= 26, 'end'= 29, 'color'= 'rgb(255,255,255)', 
       'bgcolor'= 'rgb(0,0,255)', 'tooltip'= 'Beta strand', 'underscore'= TRUE),
  list('start'= 33, 'end'= 43, 'color'= 'rgb(0,0,0)',
       'bgcolor'= 'rgb(100,100,200)', 'tooltip'= 'Helix', 'underscore'= TRUE),
  list('start'= 44, 'end'= 46, 'color'= 'rgb(0,0,0)',
       'bgcolor'= 'rgb(100,100,200)', 'tooltip'= 'Helix', 'underscore'= TRUE),
  list('start'= 48, 'end'= 50, 'color'= 'rgb(255,255,255)',
       'bgcolor'= 'rgb(0,0,255)', 'tooltip'= 'Beta strand', 'underscore'= TRUE),
  list('start'= 56, 'end'= 58, 'color'= 'rgb(255,255,255)',
       'bgcolor'= 'rgb(0,0,255)', 'tooltip'= 'Beta strand', 'underscore'= TRUE),
  list('start'= 59, 'end'= 66, 'color'= 'rgb(0,0,200)',
       'bgcolor'= 'rgb(200,200,0)', 'tooltip'= 'Turn', 'underscore'= FALSE),
  list('start'= 74, 'end'= 76, 'color'= 'rgb(255,255,255)',
       'bgcolor'= 'rgb(0,0,255)', 'tooltip'= 'Beta strand', 'underscore'= TRUE),
  list('start'= 79, 'end'= 81, 'color'= 'rgb(0,0,0)',
       'bgcolor'= 'rgb(100,100,200)', 'tooltip'= 'Helix', 'underscore'= TRUE),
  list('start'= 84, 'end'= 86, 'color'= 'rgb(0,0,200)',
       'bgcolor'= 'rgb(200,200,0)', 'tooltip'= 'Turn', 'underscore'= FALSE),
  list('start'= 91, 'end'= 97, 'color'= 'rgb(0,0,0)',
       'bgcolor'= 'rgb(100,100,200)', 'tooltip'= 'Helix', 'underscore'= TRUE),
  list('start'= 98, 'end'= 101, 'color'= 'rgb(255,255,255)',
       'bgcolor'= 'rgb(0,0,255)', 'tooltip'= 'Beta strand', 'underscore'= TRUE),
  list('start'= 102, 'end'= 106, 'color'= 'rgb(0,0,0)',
       'bgcolor'= 'rgb(100,100,200)', 'tooltip'= 'Helix', 'underscore'= TRUE),
  list('start'= 107, 'end'= 109, 'color'= 'rgb(0,0,200)',
       'bgcolor'= 'rgb(200,200,0)', 'tooltip'= 'Turn', 'underscore'= FALSE)
)

header_colors <- list(
  bg_color = "#C50063",
  font_color = "white"
)


app$layout(htmlDiv(children = list(
  htmlDiv(
    id = "app-page-header",
    style = list(
      width = "100%",
      background = header_colors[["bg_color"]],
      color = header_colors[["font_color"]]
    ),
    children = list(
      htmlA(
        id = "dashbio-logo",
        children = list(htmlImg(src = "assets/plotly-dash-bio-logo.png")),
        href = "/Portal"
      ),
      htmlH2("Dash Sequence Viewer"),
      htmlA(
        id = "gh-link",
        children = list("View on GitHub"),
        href = "https://github.com/plotly/dash-bio/blob/master/tests/dashbio_demos/app_sequence_viewer.py",
        style = list(color = "white", border = "solid 1px white")
      ),
      htmlImg(src = "assets/github.png")
    )
  ),
  htmlBr(),
  htmlDiv(
    id = "seq-view-body",
    className = 'app-body',
    children = list(
      htmlDiv(id = 'seq-view-container',
              children = list(
                htmlDiv(id = 'seq-view-component-container',
                        children = list(
                          dashbioSequenceViewer(id = "sequence-viewer",
                                                sequenceMaxHeight = '100px')
                        )),
                htmlDiv(
                  id = 'seq-view-info-container',
                  children = htmlDiv(
                    id = 'seq-view-info',
                    children = list(
                      htmlDiv(id = 'seq-view-info-desc',
                              children = list(
                                htmlSpan("Description",
                                         className = 'seq-view-info-element-title'),
                                htmlDiv(id = 'desc-info',
                                        children = list())
                              )),
                      
                      htmlBr(),
                      
                      htmlDiv(id = 'seq-view-info-aa-comp',
                              children = list(
                                htmlSpan("Amino acid composition",
                                         className = 'seq-view-info-element-title'),
                                htmlDiv(id = 'test-selection')
                              )),
                      
                      htmlBr(),
                      
                      htmlDiv(id = 'seq-view-info-coverage-clicked',
                              children = list(
                                htmlSpan("Coverage entry clicked",
                                         className = 'seq-view-info-element-title'),
                                htmlDiv(id = 'test-coverage-clicked')
                              )),
                      
                      htmlBr(),
                      htmlDiv(id = 'seq-view-info-mouse-selection',
                              children = list(
                                htmlSpan("Mouse selection",
                                         className = 'seq-view-info-element-title'),
                                htmlDiv(id = 'test-mouse-selection')
                              )),
                      
                      htmlBr(),
                      
                      htmlDiv(id = 'seq-view-info-subpart-sel',
                              children = list(
                                htmlSpan("Subpart selected",
                                         className = 'seq-view-info-element-title'),
                                htmlDiv(id = 'test-subpart-selection')
                              ))
                    )
                  )
                )
              )),
      #end seq-view-container
      
      htmlDiv(
        id = 'seq-view-control-tabs',
        className = 'control-tabs',
        children = list(dccTabs(
          id = 'seq-view-tabs',
          value = 'what-is',
          children = list(
            dccTab(
              label = 'About',
              value = 'what-is',
              children = htmlDiv(
                className = 'control-tab',
                children = list(
                  htmlH4(className = 'what-is', children = 'What is Sequence Viewer?'),
                  htmlP(
                    'Sequence Viewer is a component that allows you
                to display genomic and proteomic sequences. In
                this app, you can choose to view one of the preloaded
                data sets or upload your own FASTA file in the "Data"
                tab. For FASTA files with multiple entries, the entry
                to display in the component can be selected in the
                "Sequence" tab.'
                  ),
                  htmlP(
                    'In the "Sequence" tab, you can also select a region
                of the sequence to highlight and view its amino
                acid composition in the box under the component.'
                  ),
                  htmlP(
                    'You can additionally create a sequence coverage
                (i.e., a collection of subsequences to highlight
                and annotate). These subsequences can be extracted
                from your mouse selection, or from the results of the
                search that is shown in the component. Upon clicking on
                a coverage entry, you will be able to see the
                annotation that you have provided.'
                  )
                )
              )
            ),
            dccTab(
              label = 'Data',
              value = 'data',
              children = htmlDiv(
                className = 'control-tab',
                children = list(
                  htmlDiv(
                    id = 'preloaded-and-uploaded-alert',
                    className = 'app-controls-desc',
                    children = list(
                      'You have uploaded your own data. In order
                to view it, please ensure that the "preloaded
                sequences" dropdown has been cleared.'
                    ),
                    style = list('display' = 'none')
                  ),
                  
                  htmlDiv(
                    className = 'app-controls-block',
                    children = list(
                      htmlDiv("Preloaded sequences",
                              className = 'app-controls-name'),
                      dccDropdown(
                        className = 'app-dropdown',
                        id = 'preloaded-sequences',
                        options = list(
                          list(label = "insulin", value = "./data/sequence_viewer_P01308.fasta"),
                          list(label = "keratin", value = "./data/sequence_viewer_P04264.fasta"),
                          list(label = "albumin", value = "./data/sequence_viewer_NX_P02768.fasta"),
                          list(label = "myosin (gene)", value = "./data/sequence_viewer_myosin.fasta"),
                          list(label = "HflX (gene)", value = "./data/sequence_viewer_hflx.fasta")
                        ),
                        value = "./data/sequence_viewer_P01308.fasta"
                      ),
                      htmlDiv(id = "testdiv11")
                    )
                  ),
                  
                  htmlDiv(id = 'seq-view-fasta-upload',
                          children = list(
                            dccUpload(
                              id = 'upload-fasta-data',
                              className = 'control-upload',
                              children = htmlDiv("Drag and drop or click to upload a file.")
                            )
                          )),
                  
                  htmlA(
                    children = htmlButton(
                      "Download sample FASTA data",
                      id = 'seq-view-download-sample-data',
                      className = 'control-download',
                    ),
                    href = "./assets/sample_data/sequence_viewer_tubulin.fasta",
                    download = "tubulin.fasta"
                  )
                )
              )
            ),
            #end tab
            
            dccTab(
              label = 'Sequence',
              value = 'sequence',
              children = htmlDiv(
                className = 'control-tab',
                children = list(
                  htmlDiv(
                    id = 'seq-view-entry-dropdown-container',
                    className = 'app-controls-block',
                    children = list(
                      htmlDiv("View entry:",
                              className = 'app-controls-name'),
                      
                      dccDropdown(
                        className = 'app-dropdown',
                        id = 'fasta-entry-dropdown',
                        options = list(list('label' = 1, 'value' = 1)),
                        value = 1
                      ),
                      
                      htmlDiv(id = 'seq-view-number-entries',
                              className = 'app-controls-desc')
                    )
                  ),
                  
                  htmlBr(),
                  
                  htmlDiv(id = 'seq-view-sel-or-cov-container',
                          children = list(
                            htmlDiv("Selection or coverage:",
                                    className = 'app-controls-name'),
                            dccRadioItems(
                              id = 'selection-or-coverage',
                              options = list(
                                list('label' = 'selection', 'value' = 'sel'),
                                list('label' = 'coverage', 'value' = 'cov')
                              ),
                              value = 'sel'
                            )
                          )),
                  
                  htmlHr(),
                  htmlDiv(
                    id = 'cov-options',
                    children = list(
                      htmlDiv(
                        className = 'app-controls-block',
                        children = list(
                          htmlDiv("Add coverage from selection by:",
                                  className = 'app-controls-name'),
                          
                          dccRadioItems(
                            id = 'mouse-sel-or-subpart-sel',
                            options = list(
                              list('label' = 'mouse', 'value' = 'mouse'),
                              list('label' = 'search', 'value' = 'subpart')
                            ),
                            value = 'mouse'
                          )
                        )
                      ),
                      
                      htmlDiv(
                        className = "app-controls-block",
                        children = list(
                          htmlDiv("Text color:",
                                  className = 'app-controls-name'),
                          
                          dccInput(id = 'coverage-color',
                                   type = 'text',
                                   value = 'rgb(255, 0, 0)')
                        )
                      ),
                      
                      htmlDiv(
                        className = "app-controls-block",
                        children = list(
                          htmlDiv("Background color:",
                                  className = 'app-controls-name'),
                          
                          dccInput(id = 'coverage-bg-color',
                                   type = 'text',
                                   value = 'rgb(0, 0, 255)')
                        )
                      ),
                      
                      htmlDiv(
                        className = "app-controls-block",
                        children = list(
                          htmlDiv("Tooltip:",
                                  className = 'app-controls-name'),
                          
                          dccInput(
                            id = 'coverage-tooltip',
                            type = 'text',
                            value = '',
                            placeholder = 'hover text'
                          )
                        )
                      ),
                      
                      htmlDiv(
                        className = "app-controls-block",
                        children = list(
                          htmlDiv("Underscore text: ",
                                  className = 'app-controls-name'),
                          
                          dccChecklist(
                            id = 'coverage-underscore',
                            options = list(list('label' = '', 'value' =
                                                  'underscore')),
                            value = list()
                          )
                        )
                      ),
                      
                      htmlDiv(
                        className = "app-controls-block",
                        children = list(
                          htmlButton(id = 'coverage-submit',
                                     children = 'Submit'),
                          htmlButton(id = 'coverage-reset',
                                     children = 'Reset')
                        )
                      )
                    )
                  ),
                  #end cov-options div
                  htmlDiv(
                    id = 'seq-view-sel-slider-container',
                    children = list(
                      htmlDiv(
                        className = 'app-controls-block',
                        children = list(
                          htmlDiv(className = 'app-controls-name',
                                  children = "Selection region:"),
                          dccRadioItems(
                            id = 'sel-slider-or-input',
                            options = list(
                              list('label' = 'slider', 'value' = 'slider'),
                              list('label' = 'input', 'value' = 'input')
                            ),
                            value = 'slider'
                          )
                        )
                      ),
                      
                      htmlDiv(className = 'app-controls-block', children =
                                list(
                                  dccRangeSlider(
                                    id = 'sel-slider',
                                    min = 0,
                                    max = 0,
                                    step = 1,
                                    value = list(0, 0)
                                  )
                                )),
                      
                      htmlDiv(className = 'app-controls-block', children =
                                list(
                                  htmlDiv(
                                    id = 'sel-region-inputs',
                                    children = list(
                                      "From: ",
                                      dccInput(
                                        id = 'sel-region-low',
                                        type = 'number',
                                        min = 0,
                                        max = 0,
                                        placeholder = "low"
                                      ),
                                      
                                      "To: ",
                                      dccInput(
                                        id = 'sel-region-high',
                                        type = 'number',
                                        min = 0,
                                        max = 0,
                                        placeholder = "high"
                                      )
                                    ),
                                    style = list('display' = 'none')
                                  )
                                )),
                      
                      htmlDiv(id = 'seq-view-dna-or-protein-container',
                              children = list(
                                htmlDiv(
                                  className = 'app-controls-block',
                                  children = list(
                                    htmlDiv(className = 'app-controls-name',
                                            children = "Translate selection from:"),
                                    dccDropdown(
                                      id = 'translation-alphabet',
                                      options = list(
                                        list('label' = 'DNA', 'value' = 'dna'),
                                        list('label' = 'RNA', 'value' = 'rna')
                                      ),
                                      value = NULL
                                    )
                                  )
                                )
                              )),
                      htmlDiv(className = 'app-controls-name',
                              children = "Selection highlight color:"),
                      
                      dccDropdown(
                        className = 'app-dropdown',
                        id = 'sel-color',
                        options = list(
                          list('label' = 'violet', 'value' = 'violet'),
                          list('label' = 'indigo', 'value' = 'indigo'),
                          list('label' = 'blue', 'value' = 'blue'),
                          list('label' = 'green', 'value' = 'green'),
                          list('label' = 'yellow', 'value' = 'yellow'),
                          list('label' = 'orange', 'value' = 'orange'),
                          list('label' = 'red', 'value' = 'red')
                        ),
                        value = 'indigo'
                      )
                    )
                  ) #end seq-view-sel-slider-container div
                )
              )
            )
          )
        ))
      ),
      dccStore(id = 'coverage-storage',
               data = initialCov),
      dccStore(id = 'clear-coverage',
               data = 0),
      dccStore(id = 'current-sequence',
               data = 0)
    )
  )
  
)))

app$callback(
  output('preloaded-sequences', 'value'),
  list(input('upload-fasta-data', 'contents'),
       state('preloaded-sequences', 'value')),
  
  remove_preloaded <- function(contents,current){
    if(!is.null(contents[[1]])){
      return(NULL)
    } else{
      return(current)
    }
  }
)

app$callback(
  output('preloaded-and-uploaded-alert', 'style'),
  list(input('preloaded-sequences','value'),
       state('upload-fasta-data', 'contents')
  ),
  display_preloaded_uploaded_warning <- function(preloaded, contents){
    answer <- list('display'='none')
    if(!is.null(contents[[1]]) && !is.null(preloaded[[1]])){
      answer <- list('display'='inline-block')
    }
    return(answer)
  }
)

# sequence viewer sequence

app$callback(
  output('sequence-viewer', 'sequence'),
  list(input('upload-fasta-data', 'contents'),
       input('fasta-entry-dropdown', 'value'),
       input('preloaded-sequences', 'value')),
  
  update_sequence <- function(upload_contents, entry, preloaded){
    if(is.null(entry)){
      return("")
    }
    
    if(!is.null(preloaded[[1]])){
      protein= read_fasta(preloaded)[[entry]]
    }
    else if(!is.null(upload_contents) && is.null(preloaded[[1]])){
      data=""
      
      data <- tryCatch({
        content_type <- strsplit(upload_contents, ",")[[1]][1]
        content_string <- strsplit(upload_contents, ",")[[1]][2]
        data <- base64Decode(content_string)
      },
      error=function(e){
        data=""
      })
      
      if(data==""){
        return('-')
      }
      #protein= read_fasta(data,is_datafile=FALSE)[[entry]]
      protein= read_fasta(data,is_datafile=FALSE)[[entry]]
    } 
    else {
      return('-')
    }
    return(protein$sequence)
  }
)

# coverage

#a way of getting the timestamp for a dropdown change
app$callback(
  output('current-sequence', 'data'),
  list(input('sequence-viewer', 'sequence'),
       state('current-sequence', 'data')),
  
  signal_sequence_updated <- function(seq, current){
    if(is.null(current[[1]])){
      return(0)
    } else{ 
      return(current[[1]]+1) 
    }
  }
)

#whether or not to clear the coverage, based on a change in the dropdown
app$callback(
  output('clear-coverage', 'data'),
  list(input('current-sequence', 'modified_timestamp'),
       input('coverage-submit', 'n_clicks_timestamp'),
       state('current-sequence', 'modified_timestamp')),
  
  signal_clear_coverage<-function(c_seq, c_timestamp, s_timestamp){
    #if the coverage has been modified at all and it was modified more recently than the sequence, keep the current coverage
    answer <- 1
    if(!is.null(c_timestamp[[1]])){
      if(c_timestamp[[1]] > s_timestamp){
        answer <- 0
      } else{
        answer <- 1
      }
    }
    # if the coverage has not yet been modified, we can clear the coverage
    return(answer)
  }
)



#clear the subpart selected
app$callback(
  output('sequence-viewer', 'subpartSelected'),
  list(input('sequence-viewer', 'sequence')),
  
  clear_subpart_sel <- function(seq){
    return(list())
  }
)

# check if there is an overlap with an existing coverage item
overlaps_coverage <- function(current_cov, cov_item){
  if(length(current_cov)==0){
    return(FALSE)
  } else{
    for(i in 1:length(current_cov)){
      current_range <- c(current_cov[[i]]$start:current_cov[[i]]$end)
      item_range <- c(cov_item$start:cov_item$end)
      if(length(intersect(current_range,item_range)) > 0){
        return(TRUE) 
      } else{
        return(FALSE)
      }
    }
  }
}

app$callback(
  output('coverage-storage', 'data'),
  list(input('coverage-submit', 'n_clicks'),
       input('coverage-reset', 'n_clicks'),
       input('clear-coverage', 'data'),
       state('preloaded-sequences', 'value'),
       state('coverage-storage', 'data'),
       state('mouse-sel-or-subpart-sel', 'value'),
       state('sequence-viewer', 'mouseSelection'),
       state('sequence-viewer', 'subpartSelected'),
       state('coverage-color', 'value'),
       state('coverage-bg-color', 'value'),
       state('coverage-underscore', 'value'),
       state('coverage-tooltip', 'value'),
       state('coverage-submit', 'n_clicks_timestamp'),
       state('coverage-reset', 'n_clicks_timestamp')),
  
  edit_coverage<- function(
    s_nclicks,
    r_nclicks,
    clear_coverage,
    preloaded,
    current_cov,
    mouse_subpart,
    mouse_sel,
    subpart_sel,
    color,
    bgcolor,
    underscore,
    tooltip,
    s_timestamp,
    r_timestamp
  ){
    
    # if the coverage hasn't been updated by resetting or
    # adding, and the sequence hasn't changed from the
    # initial one, then we return the initial coverage
    if( is.null(s_nclicks[[1]]) && is.null(r_nclicks[[1]]) && (grepl('P01308',preloaded)) ){
      return(initialCov)
    }
    
    if(!is.null(r_timestamp[[1]]) && (is.null(s_timestamp[[1]]) || (s_timestamp[[1]] < r_timestamp[[1]]) )){
      return(list())
    }
    
    if(clear_coverage == 1){
      return(list())
    }
    if(mouse_subpart=='mouse'){
      if(!is.null(mouse_sel[[1]]) && !is.null(color) ){
        if(!overlaps_coverage(current_cov,mouse_sel)){
          current_cov <- append(current_cov,
                                list(list('start'=mouse_sel$start - 1,
                                          'end'=mouse_sel$end,
                                          'color'=color,
                                          'bgcolor'=bgcolor,
                                          'underscore'=ifelse(length(underscore)>0,TRUE,FALSE),
                                          'tooltip'=tooltip
                                )))
        }
      }
    }
    else if(mouse_subpart=='subpart'){
      cov_items = list()
      for(subpart in subpart_sel){
        cov_items <- append(cov_items,
                            list(list('start'=subpart$start-1,
                                      'end'=subpart$end,
                                      'color'=color,
                                      'bgcolor'=bgcolor,
                                      'underscore'=ifelse(length(underscore)>0,TRUE,FALSE),
                                      'tooltip'=tooltip
                            )))
      }
      for(cov_item in cov_items){
        if(!overlaps_coverage(current_cov,cov_item)){
          current_cov <- append(current_cov, cov_item)
        }
      }
    }
    #sort so that the tooltips can match up
    current_cov[order((sapply(current_cov,function(x) x$start)),decreasing=TRUE)]
    
    return(current_cov)
  }
)

app$callback(
  output("sequence-viewer","coverage"),
  list(input('coverage-storage','data'),
       input('selection-or-coverage','value')),
  
  apply_coverage <- function(coverage_stored, sel_or_cov){
    if(sel_or_cov != "cov"){
      return(list())
    }
    return(coverage_stored)
  }
)

#selection
app$callback(
  output('sel-slider','value'),
  list(input('sequence-viewer','sequence')),
  
  reset_selection <- function(seq){
    return(list(0,0))
  }
)

app$callback(
  output("sel-region-low",'value'),
  list(input('sequence-viewer','sequence')),
  
  reset_selection_low <- function(low){
    return(0)
  }
)

app$callback(
  output("sel-region-high","value"),
  list(input("sequence-viewer",'sequence')),
  
  reset_selection_high <- function(high){
    return(0)
  }
)


app$callback(
  output('sequence-viewer', 'selection'),
  list(input('sel-slider', 'value'),
       input('sel-region-low', 'value'),
       input('sel-region-high', 'value'),
       input('selection-or-coverage', 'value'),
       input('sel-color', 'value'),
       state('sel-slider-or-input', 'value')),
  
  update_sel <- function(slider_value, sel_low, sel_high, sel_or_cov, color, slider_input){
    answer <- list()
    if(sel_low > sel_high){
      sel_low <- 0
      sel_high <- 0
    }
    
    if(sel_or_cov != "sel"){
      answer<-list()
    } else{
      if(is.null(color)){
        color<-'blue'
      }
      if(slider_input == "slider"){
        answer <- list(slider_value[[1]],slider_value[[2]],color)
      }
      else if(slider_input == "input"){
        answer <- list(sel_low,sel_high,color)
      }
    }
    return(answer)
  }
)

# clear mouse selection

app$callback(
  output('sequence-viewer', 'mouseSelection'),
  list(input('sequence-viewer', 'sequence')),
  clear_mouse_selection <- function(seq){
    return(NULL)
  }
)

# controls

app$callback(
  output('sel-slider', 'style'),
  list(input('sel-slider-or-input', 'value')),
  
  show_hide_slider <- function(slider_input){
    if(slider_input=='slider'){
      return(list('display'='block'))
    } else{
      return(list('display'='none'))
    }
  }
)

app$callback(
  output('sel-region-inputs', 'style'),
  list(input('sel-slider-or-input', 'value')),
  
  show_hide_inputs <- function(slider_input){
    if(slider_input=='input'){
      return(list('display'='block'))
    } else{
      return(list('display'='none'))
    }
  }
)

app$callback(
  output('cov-options', 'style'),
  list(input('selection-or-coverage', 'value')),
  
  show_cov_options <- function(sel_or_cov){
    answer <- list('display'='none')
    if(sel_or_cov == 'cov'){
      answer <- list('display'='inline-block')
    }
    return(answer)
  }
)

app$callback(
  output('seq-view-sel-slider-container', 'style'),
  list(input('selection-or-coverage', 'value')),
  
  enable_disable_slider <- function(sel_or_cov){
    answer <- list('display'='none')
    if(sel_or_cov == 'sel'){
      answer <- list('display'='inline-block')
    }
    return(answer)
  }
)

app$callback(
  output('seq-view-number-entries', 'children'),
  list(input('fasta-entry-dropdown', 'options')),
  
  update_num_entries <- function(entries){
    return(paste0("Number of entries: ", length(entries)))
  }
)

app$callback(
  output('fasta-entry-dropdown', 'options'),
  list(input('upload-fasta-data', 'contents'),
       input('preloaded-sequences', 'value')),
  
  update_protein_options <- function(upload_contents, preloaded){
    dropdown_options = list(list("label"=1,'value'=1))
    
    #proteins= read_fasta("./data/sequence_viewer_NX_P02768.fasta")
    #proteins = read_fasta("./data/sequence_viewer_P01308.fasta")
    #proteins = read_fasta(preloaded)[[1]]
    #return(as.character(proteins$accession))
    
    if(!is.null(preloaded[[1]])){
      proteins = read_fasta(preloaded)
    }
    else if(!is.null(upload_contents) && is.null(preloaded[[1]])){
      data= ""
      
      data <- tryCatch({
        content_type <- strsplit(upload_contents, ",")[[1]][1]
        content_string <- strsplit(upload_contents, ",")[[1]][2]
        data <- base64Decode(content_string)
      },
      error=function(e){
        data <-""
      })
      # content_type <- strsplit(upload_contents, ",")[[1]][1]
      # content_string <- strsplit(upload_contents, ",")[[1]][2]
      # data <- base64Decode(content_string)
      
      proteins= read_fasta(data,is_datafile=FALSE)[[1]]
    }
    else{
      return(dropdown_options)
    }
    
    if(class(proteins)=='list'){
      dropdown_options <- list()
      for(i in 1:(length(proteins))){
        dropdown_options <- append(dropdown_options,list(list("label"=i,"value"=i)))
      }
      return(dropdown_options)
    }
  }
)

app$callback(
  output('sel-slider', 'max'),
  list(input('sequence-viewer', 'sequence')),
  
  update_slider_value <- function(seq){
    if(is.null(seq)){
      seq=""
    }
    return(nchar(seq))
  }
)

app$callback(
  output('sel-region-high', 'max'),
  list(input('sequence-viewer', 'sequence')),
  
  update_sel_low_max <- function(seq){
    if(is.null(seq)){
      seq=""
    }
    return(nchar(seq))
  }
)

app$callback(
  output('sel-region-low', 'max'),
  list(input('sequence-viewer', 'sequence')),
  
  update_sel_high_max <- function(seq){
    if(is.null(seq)){
      seq=""
    }
    return(nchar(seq))
  }
)

app$callback(
  output('sequence-viewer', 'title'),
  list(input('sequence-viewer', 'sequence'),
       input('fasta-entry-dropdown', 'value'),
       input('preloaded-sequences', 'value'),
       state('upload-fasta-data', 'contents')),
  
  update_sequence_title <- function(seq, entry, preloaded, upload_contents){
    if (is.null(entry)){
      return("")
    }
    
    if(!is.null(preloaded[[1]])){
      protein= read_fasta(preloaded)[[entry]]
    }
    else if(!is.null(upload_contents) && is.null(preloaded[[1]])){
      data=""
      data <- tryCatch({
        content_type <- strsplit(upload_contents, ",")[[1]][1]
        content_string <- strsplit(upload_contents, ",")[[1]][2]
        data <- base64Decode(content_string)
      },
      error=function(e){
        data <- ""
      })
      
      if(data==""){
        return('')
      }
      
      protein= read_fasta(data,is_datafile=FALSE)[[entry]]
    }
    else {
      return('')
    }
    
    titles = list('Name', 'Entry Name', 'Protein Name', 'Identifier', 'Desc-0')
    
    for(t in titles){
      ifelse(is.null(protein$description[[t]]), next, return(protein$description[[t]]))
    }
    return('')
  }
)


# app$callback(
#   output('testdiv11', 'children'),
#   list(input('sequence-viewer', 'selection')),
# 
#   testq<- function(dom){
#     
#     if(is.null(dom) || length(dom) < 2){
#       print("ok")
#     }
#     print(dom)
#     print(".")
#     print(dom[[1]])
#     print("..")
#     print(dom[[2]])
#     return(as.character(dom))
#   }
# )


# info display
app$callback( 
  output('test-selection', 'children'),
  list(input('sequence-viewer', 'selection'),
       input('translation-alphabet', 'value'),
       input('sequence-viewer', 'sequence')),

  get_aa_comp <- function(v, alphabet, seq){
    answer <- ""
    break_and_return = FALSE

    if(is.null(v) || length(v) < 2){
      warning("pass")
    } else {
      res <- tryCatch({
        substr(seq,v[[1]],v[[2]])
      },
      error=function(e){
        "error"
      })

      if(res!="error"){
        subsequence <- res
      } else{
        answer <- htmlTable(list())
        break_and_return <- TRUE
      }

      if(break_and_return==FALSE){

        # default - file represents a protein
        aa_string <- subsequence

        if(grepl("dna", alphabet)){
          print(alphabet)
          # remove partial codons
          subsequence <- substr(subsequence,1,(nchar(subsequence) - nchar(subsequence)%%3))
          print("subsequence ok")

          if(nchar(subsequence)%%3 != 0){
            print("OK")
          } else{
            subsequence <- subsequence
            print(subsequence)
          }

          #ok <-as.character(translate(DNAString(x="AATAATCATATTATCAGACAAGGTAATACTGTACACAATGAATATAAAACACAGTTACATGCAACTTCATGGAAAAATCTCACAACAATAATTTTGGGCAAAAGAAGTAAGATACTAAAAAATAAATTCAGTTGTATTCACATGAAGTTCAAAACCAGCAAAACTAAACTAGGGATGCATATACAAATTTTAAAACATTAAGAAATGCAGAGAAATAATTATCATAAAAATAAATTCATAGTAACCTTTT"),if.fuzzy.codon="solve"))
          #nchar(substr(ok,1,(nchar(ok) - nchar(ok)%%3)))

          res_dna <- tryCatch({
            as.character(translate(DNAString(x=subsequence),if.fuzzy.codon = "solve"))
          },
          error=function(e){
            "Sequence does not represent DNA"
          },
          warning=function(e){
            as.character(translate(DNAString(x=subsequence),if.fuzzy.codon = "solve"))
          })

          if(grepl("Sequence does not represent DNA",res_dna)){
            answer <- "Sequence does not represent DNA"
            break_and_return <- TRUE
          } else{
            aa_string <- res_dna
          }
        } #end DNA if 

        else if(alphabet == 'rna'){
          # remove partial codons
          subsequence <- substr(subsequence,1,(nchar(subsequence) - nchar(subsequence)%%3))

          if(nchar(subsequence)%%3 != 0){
            print("OK")
          } else{
            subsequence <- subsequence
          }

          #In Biopython, generic_rna and translate method automatically substitute thymines (T) for uracil (U).
          #In Biostrings, it doesn't do that, so to add that functionality, a if statement is added.
          if(grepl("T",subsequence)){
            res_rna <- tryCatch({
              as.character(translate(RNAString(DNAString(x=subsequence)),if.fuzzy.codon = "solve"))
            },
            error=function(e){
              "Sequence does not represent RNA"
            },
            warning=function(e){
              as.character(translate(RNAString(DNAString(x=subsequence)),if.fuzzy.codon = "solve"))
            })
          } else {
            res_rna <- tryCatch({
              as.character(translate(RNAString(x=subsequence),if.fuzzy.codon = "solve"))
            },
            error=function(e){
              "Sequence does not represent RNA"
            },
            warning=function(e){
              as.character(translate(RNAString(x=subsequence),if.fuzzy.codon = "solve"))
            })
          }

          if(grepl("Sequence does not represent RNA",res_rna)){
            answer <- "Sequence does not represent RNA"
            break_and_return <- TRUE
          } else{
            aa_string <- res_rna
          }
        } #end RNA if


        if(break_and_return==FALSE){
          # all unique amino acids
          amino_acids <- unique(strsplit(aa_string, "")[[1]])
    
          #In Python version, the seq3() function automatically maps "*" to "TER" but the AMINO_ACID_CODE from BioStrings in R doesn't do that, so we edit the function to do that:
          custom_AMINO_ACID_CODE <- AMINO_ACID_CODE
          custom_AMINO_ACID_CODE[[27]] <- "Ter"
          names(custom_AMINO_ACID_CODE)[27] <- "\\*"
          
          aa_counts <-lapply(amino_acids, function(x){
            if(grepl("*",x,fixed=TRUE)){
              x <- "\\*"
            }
            list('aa' = custom_AMINO_ACID_CODE[x][[1]], 'count' = as.numeric(str_count(aa_string,x)))
          })

          # sort by most common AA in sequence
          aa_counts <- aa_counts[order((sapply(aa_counts,function(x) x[[2]])),decreasing=TRUE)]


          summary <- lapply(aa_counts, function(x){
            htmlTr(list(htmlTd(x$aa),
                        htmlTd(as.character(x$count))))
          })

          if((grepl("dna|rna", alphabet)) && (length(summary)>0)){
            answer <- list(sprintf("Protein translated from %s: %s",toupper(alphabet),aa_string),htmlTable(summary))
          } else{
            answer <- htmlTable(summary)
          }
        }
      } #end if
    } #end else
    return(answer)
  }
)

app$callback(
  output('test-coverage-clicked', 'children'),
  list(input('sequence-viewer', 'coverageClicked'),
       state('coverage-storage', 'data')),
  update_cov_clicked<-function(index, current_cov){
    if (is.null(index[[1]]) || is.null(current_cov)){
      return("")
    }
    index <- index+1
    return(sprintf("Index: %d Tooltip: %s",index,current_cov[[index]]$tooltip))
  }
)

app$callback(
  output('test-mouse-selection', 'children'),
  list(input('sequence-viewer', 'mouseSelection')),

  update_mouse_sel<-function(v){
    if(!is.null(v)){
      return(v$selection)
    }
    return("")
  }
)

app$callback(
  output('fasta-entry-dropdown', 'value'),
  list(input('preloaded-sequences', 'value'),
       input('upload-fasta-data', 'contents')),

  update_dropdown_value<-function(v, c){
    return(1)
  }
)

app$callback(
  output('desc-info','children'),
  list(input('upload-fasta-data','contents'),
       input('fasta-entry-dropdown','value'),
       input('preloaded-sequences', 'value')),

  update_desc_info <- function(upload_contents, entry, preloaded){
    if(is.null(entry)){
      return('Please select an entry.')
    }

    if(!is.null(preloaded[[1]])){
      protein= read_fasta(preloaded)[[entry]]
    }
    else if(!is.null(upload_contents) && is.null(preloaded[[1]])){
      data=""

      data <- tryCatch({
        content_type <- strsplit(upload_contents, ",")[[1]][1]
        content_string <- strsplit(upload_contents, ",")[[1]][2]
        data <- base64Decode(content_string)
      },
      error=function(e){
        data <- ""
      })


      if(data==""){
        return(list())
      }
      protein <- tryCatch({
        read_fasta(data,is_datafile=FALSE)[[entry]]
      },
      error=function(e){
        return("NA")
      })

    }
    else {
      return('Please either upload a file or select one from the dropdown.')
    }

    desc <- list()
    for(i in names(protein$description)){

      if(!grepl("desc-",i)){
        tmp <- sprintf("%s:",i)
      }else{
        tmp <- "-"
      }
      tmp <- paste(tmp, protein$description[[i]],sep= " ")
      desc <- append(desc, tmp)
      desc <- append(desc, list(htmlBr()))
    }
    return(desc)
  }
)

app$callback(
  output('test-subpart-selection', 'children'),
  list(input('sequence-viewer', 'subpartSelected')),

  update_subpart_sel <- function(v){
    if(is.null(v)){
      return('')
    }
    test <- list()
    for(sel in v){
      if(nchar(sel$sequence)==0){
        next
      }
      test <- append(test, sprintf("Start: %d ",sel$start))
      test <- append(test, sprintf("End: %d ",sel$end))
      test <- append(test, sprintf("Sequence: %s",sel$sequence))
      test <- append(test, list(htmlBr()))
    }
    return(test)
  }
)

if (appName != "") {
  
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
  
} else {
  
  app$run_server(showcase = TRUE)
  
}
