appName <- Sys.getenv("DASH_APP_NAME")
if (appName != "") {
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

# set working directory
setwd("/app")

library(dash)
library(dashBio)
library(dashCoreComponents)
library(dashHtmlComponents)
library(Biostrings)
library(readr)
library(ape)
library(stringr)
library(rentrez)
library(dashTable)
library(plotly)

source("assets/helpers.R")

app <- Dash$new()

# Functions and Data Loading

FASTA_DATA <- readBStringSet(filepath = "data/alignment_viewer_p53_clustalo.fasta")

blank_dataframe <- data.frame(matrix(ncol=0,nrow=0))


# App Layout Elements

#Header

header <- htmlDiv(
  id = "app-page-header",
  style = list(
    width = "100%",
    background = "#262B3D" ,
    color = "#E2EFFA"
  ),
  children = list(
    htmlA(
      id = "dashbio-logo",
      children = list(
        htmlImg(src='assets/plotly-dash-bio-logo.png', height = '36', width = '180',
                style = list('top' = '10', 'margin' = '10px'))
      ),
      href = "/Portal"
    ),
    htmlH2("NCBI Nucleotide DB Explorer"),
    htmlA(
      id = "gh-link",
      children = list("View on GitHub"),
      href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-ncbi-explorer",
      style = list(color = "white", border = "solid 1px white")
    ),
    htmlImg(
      src = "assets/GitHub-Mark-Light-64px.png"
    )
  )
)

#First Row 

first_row <- htmlDiv(list(
  htmlDiv(children = list(
    htmlH3("Search for a Dataset:", style = list("margin" = 10)),
    htmlP("Use keywords to search the Nucleotide database for a gene or organism of interest.
                    This will retrieve the top 10 related datasets' Accession IDs along with their descriptions.
                    ", style = list("margin" = 10)),
    htmlP('Example searches: "Basiliscus basiliscus[Organism]" or "BRCA1[Gene]"
                    ', style = list("margin" = 10)),
    htmlP('Select the table cell containing your desired Accession ID to add this dataset
    to the alignment and generate a new figure Then, filter through datasets to analyze the CpG patterns
    and nucleotide compositions for your sequence of interest.
                    ', style = list("margin" = 10)),
    dccInput(
      id = 'search-input',
      placeholder = "Enter a species or gene...",
      type = "text",
      value = "",
      style = list("margin" = "10px")
    ),
    htmlButton(
      id = 'search-button', n_clicks = 0, children = 'Search for Datasets'
    ),
    htmlDiv(id = "search-output"),
    dashDataTable(
      id = "table",
      columns = lapply(colnames(blank_dataframe), 
                       function(colName){
                         list(
                           id = colName,
                           name = colName
                         )
                       }),
      data = df_to_list(blank_dataframe)
    )
  ),className = 'item-c'),
  
  htmlDiv(dccTabs(id = 'circos-tabs', value = 'what-is', children = list(
    dccTab(
      label = 'About',
      value = 'what-is',
      children = htmlDiv(
        id = 'control-tab', children = list(
          htmlH4(className = 'what-is', children = "What is NCBI Nucleotide DB Explorer?"),
          
          htmlP("The NCBI Nucleotide sequence database is an open access, annotated 
                    collection of all publicly available nucleotide sequences and
                    their protein translations. This database is produced and 
                    maintained by the National Center for Biotechnology Information.
                    ", style = list("padding" = "5px")),
          htmlP("This app will allow you to explore, analyze and compare datasets on the
          Nucleotide database. Load datasets using the GenBank Accession ID's directly in
          the 'Explorer Controls' tab, or find the top datasets corresponding to a species or gene of interest
          in the 'Database Search'.
                    ", style = list("padding" = "5px")),
          htmlP("Generate and analyze the alignment, sequence composition, and CpG island
          distribution of the datasets, or export the data as a FASTA file for downstream analysis.
                    ", style = list("padding" = "5px")),
          dccMarkdown("The NCBI Nucleotide Database Explorer incorporates data retrieved from NCBI. 
          Please consult the [NCBI Website and Data Usage Policies and Disclaimers](https://www.ncbi.nlm.nih.gov/home/about/policies/) for further information.
                    ", style = list("padding-left" = "5px"))
          )
      )
    ),
    
    dccTab(
      label = 'Explorer Controls',
      value = 'data',
      children = htmlDiv(className = 'control-tab', children = list(
        htmlDiv(className = 'app-controls-block', children = list(
          htmlP("Enter a single GenBank accession ID or multiple ID's separated by commas
                    and select the button to generate an alignment and the sequence of the dataset.
                    Generating the alignment will also create a FASTA file which you can download.
                    Enter  additional datasets to add these sequences to the alignment. Clear your FASTA file with the button below.", style = list("margin" = 10)),
          htmlP('Example: Single Dataset - NR_108049', style = list("margin" = 10)),
          htmlP(' Multiple Datasets - HM161150, FJ356743, 
                    JQ073190, GU457971, FJ356741', style = list("margin" = 10)),
          dccInput(
            id = 'genbank-input',
            placeholder = "Enter an Accession ID...",
            type = "text",
            value = "",
            style = list("margin" = 10)
          ),
          htmlButton(
            id = 'submit-button', n_clicks = 0, children = 'Generate Sequence'
          ),
          htmlButton(
            id = 'submit-button-2', n_clicks = 0, children = 'Generate Alignment'
          ),
          htmlButton(
            id = 'reset-file', n_clicks = 0, children = 'Clear FASTA File'
          ),
          htmlA(
            htmlButton(
              "Download FASTA data",
              className="control-download"
            ),
            href="data/random_fasta.fasta",
            download="random_fasta.fasta"
          ),
          
          htmlDiv(id='alignment-file-upload-container', children=list(
            dccUpload(
              id='alignment-file-upload',
              className='control-upload',
              children=htmlDiv(list(
                "Drag and drop FASTA files or select files."
              ))
            )
          )),
          
          dccStore(
            id = 'genbank-sequence'),
          
          dccStore(
            id = 'fasta-file'
          ),
          
          dccStore(
            id = 'dataframe-store'
          )
        ))
      ))
    )
  )), className = 'item-a'),
  
  htmlDiv(list(
    htmlP("Multiple Sequence Alignment Chart", 
          style = list("font-family" = "Open Sans", "font-size" = "22px", "color" = "#262B3D", "text-align" = "center")),
    dccLoading(type = "dot", children = htmlDiv(
      id = 'alignment-container',
      children = list(),
      style = list("overflow-x" = "scroll", "display" = "flex", "justify-content" = "center")
    ))), className = 'item-b'),
  
  htmlDiv(children = list(
    htmlP("Sequence Viewer", style = list("font-family" = "Open Sans", "font-size" = "22px", "color" = "#262B3D", "text-align" = "center")),
    dashbioSequenceViewer(
      id = 'sequence-viewer',
      sequence = as.character(FASTA_DATA$`sp|Q9W678|P53_BARBU Cellular tumor antigen p53 OS=Barbus barbus GN=tp53 PE=2 SV=1`),
      selection = list(10,20, "green"),
      charsPerLine = 40
    )), className = 'item-d'),
  
  htmlDiv(list(
    dccGraph(id = "gc_graph")
  ), className = 'item-e'),
  
  htmlDiv(dccGraph(
    id = 'sequence-pie-chart'
  ), className = 'item-f'),
  
  htmlDiv(children = list(
    htmlP("Selected Dataset", style = list("font-family" = "Open Sans", "font-size" = "14px", "color" = "#262B3D", "text-align" = "center")),
    dccRadioItems(
      id = "cpg-options",
      options = list(list(label = "Default", value = 1)),
      value = 1),
    htmlButton(
      id = 'reset-button', n_clicks = 0, children = 'Refresh',
      style = list("width" = "100%", "padding" = "0px", "margin" = "0px", "margin-top" = "15px")
      )
    )
    ,className = 'item-g')
  
), className = "container")

app$layout(
  htmlDiv(
    list(
      header,
      first_row
    )
  )
)

# Callbacks for Sequence Viewer Output

app$callback(
  output(id = 'genbank-sequence', property = 'data'),
  params = list(
    input(id = 'submit-button', property = 'n_clicks'),
    state(id = 'genbank-input', property = 'value')
  ),
  
  update_data <- function(n_clicks, accession_id) {
    accession_id <- gsub("\\s", "", accession_id)
    accession_id <- gsub('"', "", accession_id)
    
    if (n_clicks > 0) {
      sequence_info = list()
      genbank_file <- read.GenBank(access.nb = accession_id, as.character = TRUE)
      sequence = toupper(paste(genbank_file[[1]], collapse = ''))
      sequence_info <- c(sequence_info, "sequence" = sequence)
      return(sequence_info)
    }
  }
)

app$callback(
  output(id = 'sequence-viewer', property = 'sequence'),
  params = list(
    input(id = 'genbank-sequence', property = 'data'),
    input(id ='cpg-options', property = 'value')
  ),
  
  update_sequence <- function(data, value) {
    if (value != 1) {
      dna_data <- readDNAStringSet("data/random_fasta.fasta", "fasta")
      return(as.character(dna_data[[value]]))
    }
    else{
      return(data[["sequence"]])
    }
  }
)

# Callback for Alignment-Viewer Output

app$callback(
  output(id = 'alignment-container', property = 'children'),
  params = list(
    input(id = 'submit-button-2', property = 'n_clicks'),
    input(id = 'dataframe-store', property = "data"),
    input(id = "table", property = "active_cell"),
    state(id = 'genbank-input', property = 'value')
  ),
  
  update_fasta <- function(n_clicks, search_data, cell, accession_id) {

    accession_id <- gsub("\\s", "", accession_id)
    
    accession_id <- gsub('"', "", accession_id)
    
    accession_id <- strsplit(accession_id, ",")
    
    accession_id <- unlist(accession_id)
    
    ctx = app$callback_context()
    
    if (ctx$triggered$prop_id == "submit-button-2.n_clicks") {
      bio_file <- read.GenBank(access.nb = accession_id)
      write.dna(bio_file, file ="data/random_fasta.fasta", format = "fasta", append = TRUE, nbcol = 6, colsep = "", colw = 10)
      fasta_file <- toupper(read_file("data/random_fasta.fasta"))
      alignment_chart <- dashbioAlignmentChart(id = 'alignment-chart', data = fasta_file, height = 750, width = 1150)
    }
    
    else if (!is.null(search_data[[1]]) && !is.null(cell[[1]])) {
      search_data  <-  as.data.frame(matrix(unlist(search_data), nrow=length(unlist(search_data[1]))), stringsAsFactors = FALSE)
      accession_id <- search_data[1, cell$row + 1]
      accession_id <- gsub('"', "", accession_id)
      accession_id <- gsub(" ", "", accession_id)
      
      bio_file <- read.GenBank(access.nb = accession_id)
      write.dna(bio_file, file ="data/random_fasta.fasta", format = "fasta", append = TRUE, nbcol = 6, colsep = "", colw = 10)
      fasta_file <- toupper(read_file("data/random_fasta.fasta"))
      alignment_chart <- dashbioAlignmentChart(id = 'alignment-chart', data = fasta_file, height = 750, width = 1150)
    }
    
    
    else if (n_clicks < 1) {
      alignment_chart <- dashbioAlignmentChart(
        id = 'alignment-chart',
        data = read_file("data/alignment_viewer_p53_clustalo.fasta"),
        height = 750,
        width = 1150,
        opacity = 0.5
      )
    }
    return(alignment_chart)
  }
)

#Callbacks for database search output.

app$callback(
  output(id = "search-output", property = "children"),
  params = list(
    input(id = "search-button", property = "n_clicks"),
    state(id = "search-input", property = "value")
  ),
  
  update_search_results <- function(n_clicks, search_term) {
    if (n_clicks > 0) {
      sequence_dataframe <- results_dataframe(search_term)
      
      results_table <- dashDataTable(
        id = "table",
        style_table = list('overflowX' = 'scroll', 'overflowY' = 'scroll', 'maxHeight' = '100'),
        style_data = list(
          whiteSpace = "normal"
        ),
        css = list(
          list(
            selector = '.dash-cell div.dash-cell-value',
            rule = 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'
          )
        ),
        columns = lapply(colnames(sequence_dataframe),
                         function(colName){
                           list("name" = colName,
                                "id" = colName)
                         }),
        data = df_to_list(sequence_dataframe)
      )
      
      return(results_table)
    }
  }
)

app$callback(
  output(id = "dataframe-store", property = "data"),
  params = list(
    input(id = "search-button", property = "n_clicks"),
    state(id = "search-input", property = "value")
  ),
  
  update_search_results <- function(n_clicks, search_term) {
    if (n_clicks > 0) {
      sequence_dataframe <- results_dataframe(search_term)
      
      return(sequence_dataframe)
    }
  }
)

#Callbacks for Nucleotide Pie Chart

app$callback(
  output(id = "sequence-pie-chart", property = "figure"),
  params = list(
    input(id = 'submit-button', property = 'n_clicks'),
    input(id = "sequence-viewer", property = "sequence")),
  
  update_pie_chart <- function(n_clicks, sequence) {
    if (n_clicks > 0 | !is.null(sequence[[1]])) {
      sequence_string <- as.character(sequence)
      chart_sequence <- strsplit(sequence_string, split = "")
      p <- plot_ly(as.data.frame(table(chart_sequence)), labels = ~chart_sequence, values = ~Freq, type = 'pie', hole = 0.6) %>%
        layout(paper_bgcolor="#FFFFFF", title = "Nucleotide Base Composition", margin = 10,
               title = list(
                 family = "Open Sans",
                 size = 22,
                 color = '#262B3D'),
               font = list(
                 color = '#262B3D'
               )
        )
      
      
    }
    else {
      p <- (plot_ly(as.data.frame(table(strsplit(as.character(FASTA_DATA$`sp|Q9W678|P53_BARBU Cellular tumor antigen p53 OS=Barbus barbus GN=tp53 PE=2 SV=1`) 
                                                 , split = ""))), labels = ~Var1, values = ~Freq, type = 'pie', hole = 0.6) %>%
              layout(paper_bgcolor="#ffffff", title = "Nucleotide Base Composition", margin = 10,
                     title = list(
                       family = "Open Sans",
                       size = 22,
                       color = '#262B3D'),
                     font = list(
                       color = '#262B3D'
                     )
              )
      )
    }
    return(p)
  }
)

#Callbacks for CpG Chart, Options, and Sequence Selection

app$callback(
  output(id = "cpg-options", property = "options"),
  params = list(
    input(id = "submit-button-2", property = "n_clicks"),
    input(id = 'search-button', property = 'n_clicks'),
    input(id = 'reset-button', property = 'n_clicks')
  ),
  
  select_trace <- function(n_clicks, search_clicks, reset_clicks) {
    if (n_clicks > 0 | search_clicks > 0 | reset_clicks > 0) {
      dna_data = readDNAStringSet("data/random_fasta.fasta", "fasta")
      trace_names <- unique(names(dna_data))
      
      list_options <- lapply(trace_names, function(x) {
        list(label = x, value = match(x, trace_names))
      })
    }
    else {
      list_options <- list(list(label = "Default", value = 0))
    }
    return(list_options)
  }
)

app$callback(
  output(id = "gc_graph", property = "figure"),
  params = list(
    input(id ='submit-button-2', property = 'n_clicks'),
    input(id ='cpg-options', property = 'value')
  ),
  update_gc <- function(n_clicks, option) {
    if (n_clicks > 0) {
      dna_data = readDNAStringSet("data/random_fasta.fasta", "fasta")
      window = 100
      gc = rowSums(letterFrequencyInSlidingView(dna_data[[option]], window, c("G", "C")))/window
      g = as.data.frame(gc)
      g$index = as.numeric(rownames(g))
      m <- g[which.max(g$gc), ]
      
      p <- plot_ly(as.data.frame(g), x=~index, y=~gc, type = "scatter", mode = "line", line = list(color="#262B3D")) %>%
        layout(paper_bgcolor="#ffffff", plot_bgcolor="#ffffff", title = "CpG Island Distribution", margin = 20,
               title = list(
                 family = "Open Sans",
                 size = 22,
                 color = '#262B3D'),
               font = list(
                 color = 'rgb(50,50,50)'
               ),
               annotations = list(
                 x = m$index,
                 y = m$gc,
                 text = "Select a peak to highlight the sequence above.",
                 xref = "x",
                 yref = "y",
                 showarrow = TRUE,
                 arrowhead = 7,
                 ax = 20,
                 ay = -40
               ),
               xaxis = list(title = "Sequence Index"),
               yaxis = list(title = "GC Content over 100 BP windows")
        )
    }
    
    else {
      dna_data = readDNAStringSet("data/dnastring.fasta", "fasta")
      window = 100
      gc = rowSums(letterFrequencyInSlidingView(dna_data[[option]], window, c("G", "C")))/window
      g = as.data.frame(gc)
      g$index = as.numeric(rownames(g))
      m <- g[which.max(g$gc), ]
      
      p <- plot_ly(as.data.frame(g), x=~index, y=~gc, type = "scatter", mode = "line", line = list(color="#262B3D")) %>%
        layout(paper_bgcolor="#ffffff", plot_bgcolor="#ffffff", title = "CpG Island Distribution", margin = 20,
               title = list(
                 family = "Open Sans",
                 size = 22,
                 color = '#262B3D'),
               font = list(
                 color = 'rgb(50,50,50)'
               ),
               annotations = list(
                 x = m$index,
                 y = m$gc,
                 text = "Select a peak to highlight the sequence above.",
                 xref = "x",
                 yref = "y",
                 showarrow = TRUE,
                 arrowhead = 7,
                 ax = 20,
                 ay = -40
               ),
               xaxis = list(title = "Sequence Index"),
               yaxis = list(title = "GC Content over 100 BP windows")
        )
    }
    return(p)
  }
)


app$callback(
  output(id = 'genbank-input', property = "value"),
  params = list(
    input(id = 'submit-button-2', property = "n_clicks")
  ),
  
  update_string <- function(n_clicks) {
    if(n_clicks > 0) {
      return("")
    }
  }
)


app$callback(
  output(id = "sequence-viewer", property = "selection"),
  params = list(
    input(id = "gc_graph", property = "clickData")
  ),
  update_selection <- function(hoverdata) {
    point <- hoverdata$points[[1]]$pointNumber
    selection_fixed <- list(as.numeric(point-30), as.numeric(point+30), "#262B3D")
    return(selection_fixed)
  }
)

app$callback(
  output(id = 'fasta-file', property = 'data'),
  params = list(
    input(id = 'reset-file', property = "n_clicks")
  ),
  reset_fasta <- function(n_clicks) {
    if (n_clicks > 0) {
      write.dna(" ", file ="data/random_fasta.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = "", colw = 10)
      print("FASTA CLEARED")
    }
  }
)


if(appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
} else {
  app$run_server(showcase = TRUE)
}
