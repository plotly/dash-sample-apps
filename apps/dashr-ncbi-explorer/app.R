appName <- Sys.getenv("DASH_APP_NAME")


if (appName != ""){
  
  pathPrefix <- sprintf("/%s/", appName)
  
  
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  
  
  
  setwd(sprintf("/app/apps/%s", appName))
  
}

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

app <- Dash$new()

# Functions and Data Loading

FASTA_DATA <- readBStringSet(filepath = "data/alignment_viewer_p53_clustalo.FASTA")

blank_dataframe <- data.frame(matrix(ncol=0,nrow=0))

generate_table <- function(df, nrows=20)
  
  # function generates a dash table from a supplied data frame (df)
  
  #  and number of rows (nrows) to display
{
  
  n <- min(nrows, nrow(df))
  
  rows <- lapply(seq(1, n), function(i) {
    
    htmlTr(children = lapply(as.character(df[i,]), htmlTd))
    
  })
  
  header <- htmlTr(children = lapply(names(df), htmlTh))
  
  htmlTable(
    
    children = c(list(header), rows)
  )
}

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
      href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-clustergram",
      style = list(color = "white", border = "solid 1px white")
    ),
    htmlImg(
      src = "assets/GitHub-Mark-Light-64px.png"
    )
  )
)


#First Row 

first_row <- htmlDiv(list(
  htmlDiv(dccTabs(id = 'circos-tabs', value = 'what-is', children = list(
    dccTab(
      label = 'About',
      value = 'what-is',
      children = htmlDiv(
        id = 'control-tab', children = list(
          htmlH4(className = 'what-is', children = "What is NCBI Explorer?"),
          
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
                    ", style = list("padding" = "5px"))
          
        )
      )
    ),
    
    dccTab(
      label = 'Explorer Controls',
      value = 'data',
      children = htmlDiv(className = 'control-tab', children = list(
        htmlDiv(className = 'app-controls-block', children = list(
          htmlP("Enter a single GenBank accession ID or multiple ID's seperated by commas
                    and select the button to generate an alignment or the sequence of the dataset.
                    Enter  addtional datasets to add these sequences to the alignment.", style = list("margin" = 10)),
          htmlP('Example: Single Dataset - NR_108049', style = list("margin" = 10)),
          htmlP(' Multiple Datasets - "JF806202", "HM161150", "FJ356743", "JF80620", 
                    "JQ073190", "GU457971", "FJ356741", "JF806"', style = list("margin" = 10)),
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
      children = list()
    ))), className = 'item-b'),
  
  htmlDiv(children = list(
    htmlH3("Search for a Dataset:", style = list("margin" = 10)),
    htmlP("Use keywords to search the Nucleotide database for a gene or organism of interest.
                    This will retrieve the top 10 related datasets' Accession IDs along with their descriptions.
                    ", style = list("margin" = 10)),
    htmlP('Example searches: "Basiliscus basiliscus[Organism]" or "BRCA1[Gene]
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
  
  htmlDiv(children = list(
    htmlP("Sequence Viewer", style = list("font-family" = "Open Sans", "font-size" = "22px", "color" = "#262B3D", "text-align" = "center")),
    dashbioSequenceViewer(
      id = 'sequence-viewer',
      sequence = as.character(FASTA_DATA$`sp|Q9W678|P53_BARBU Cellular tumor antigen p53 OS=Barbus barbus GN=tp53 PE=2 SV=1`),
      selection = list(10,20, "green"),
      charsPerLine = 90
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
      value = 1)),className = 'item-g')
  
), className = "container")


app$layout(
  header,
  first_row
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
    if (value != 0) {
      dna_data <- readDNAStringSet("data/random_fasta.FASTA", "fasta")
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
    
    if (n_clicks > 0) {
      bio_file <- read.GenBank(access.nb = accession_id)
      write.dna(bio_file, file ="data/random_fasta.FASTA", format = "fasta", append = TRUE, nbcol = 6, colsep = "", colw = 10)
      fasta_file <- toupper(read_file("data/random_fasta.FASTA"))
      alignment_chart <- dashbioAlignmentChart(id = 'alignment-chart', data = fasta_file, height = 750, width = 850)
    }
    
    else if (!is.null(search_data[[1]]) && !is.null(cell[[1]])) {
      search_data  <-  as.data.frame(matrix(unlist(search_data), nrow=length(unlist(search_data[1]))), stringsAsFactors = FALSE)
      accession_id <- search_data[1, cell$row + 1]
      accession_id <- gsub('"', "", accession_id)
      accession_id <- gsub(" ", "", accession_id)
      
      bio_file <- read.GenBank(access.nb = accession_id)
      write.dna(bio_file, file ="data/random_fasta.FASTA", format = "fasta", append = TRUE, nbcol = 6, colsep = "", colw = 10)
      fasta_file <- toupper(read_file("data/random_fasta.FASTA"))
      alignment_chart <- dashbioAlignmentChart(id = 'alignment-chart', data = fasta_file, height = 750, width = 850)
    }
    
    
    else if (n_clicks < 1) {
      alignment_chart <- dashbioAlignmentChart(
        id = 'alignment-chart',
        data = read_file("data/alignment_viewer_p53_clustalo.FASTA"),
        height = 750,
        width = 850,
        opacity = 0.5
      )
    }
    return(alignment_chart)
  }
)