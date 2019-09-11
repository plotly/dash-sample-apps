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

df <- read.csv(
  file = "https://raw.githubusercontent.com/plotly/datasets/master/solar.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


FASTA_DATA <- readBStringSet(filepath = "data/alignment_viewer_p53_clustalo.FASTA")


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

# Layout Elements


second_row <- htmlDiv(children =list(
  htmlDiv(dashbioSequenceViewer(
    id = 'sequence-viewer',
    sequence = as.character(FASTA_DATA$`sp|Q9W678|P53_BARBU Cellular tumor antigen p53 OS=Barbus barbus GN=tp53 PE=2 SV=1`),
    selection = list(10,20, "green"),
    charsPerLine = 90
  )),
  htmlDiv(dccGraph(
    id = 'sequence-pie-chart'
  )),
  htmlDiv(list(
    dccGraph(id = "gc_graph")
  ),style = list("flex" = "-moz-max-content")
  ),
  dccRadioItems(
    id = "cpg-options",
    options = list(list(label = "Default", value = 1)),
    value = 1,
    style = list("margin" = "10px")
  )
), style = list("display" = "flex", "width" = "98%", "margin" = "2%", "flex-wrap" = "wrap"))


options <- htmlDiv(
  id = 'circos-body',
  className = 'app-body',
  children = list(
    dccLoading(className = 'dashbio-loading', type = "dot", children = htmlDiv(
      id = 'alignment-container',
      children = list()
    )),
    
    htmlDiv(id = 'circos-control-tabs', className = 'control-tabs', children = list(
      dccTabs(id = 'circos-tabs', value = 'what-is', children = list(
        dccTab(
          label = 'About',
          value = 'what-is',
          children = htmlDiv(
            id = 'control-tab', children = list(
              htmlH4(className = 'what-is', children = "What is NCBI Explorer?"),
              
              htmlP("The GenBank sequence database is an open access, annotated 
                    collection of all publicly available nucleotide sequences and
                    their protein translations. This database is produced and 
                    maintained by the National Center for Biotechnology Information
                    as part of the International Nucleotide Sequence Database Collaboration."),
              htmlBr(),
              htmlP("This app will allow you to explore, analyze and compare datasets using
                    the GenBank Accession ID or by searching up organisms or gene's of interest
                    using keyword entries in the 'Database Search'. You can analyze the alignment, sequence
                    composition, and CpG island distribution of any of the datasets. 
                    ")
            )
          )
        ),
        
        dccTab(
          label = 'Explorer Controls',
          value = 'data',
          children = htmlDiv(className = 'control-tab', children = list(
            htmlDiv(className = 'app-controls-block', children = list(
              htmlDiv(className = 'app-controls-name', children = "GenBank Access:"),
              htmlP("Enter a single GenBank accession ID or multiple ID's seperated by commas
                    and select the button to generate an alignment or the sequence of the dataset.
                    Enter additional datasets to add additional sequences to the alignment."),
              htmlP('Example: Single Dataset - NR_108049', style = list("margin" = 10)),
              htmlP(' Multiple Datasets - "JF806202", "HM161150", "FJ356743", "JF80620", 
                    "JQ073190", "GU457971", "FJ356741", "JF806"', style = list("margin" = 10)),
              dccInput(
                id = 'genbank-input',
                placeholder = "Enter a GenBank Accession Number...",
                type = "text",
                value = "",
                style = list("margin" = "10px")
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
              )
            ))
          ))
        ),
        
        dccTab(
          label = 'Database Search',
          children = htmlDiv(className = 'control-tab', children = list(
            htmlDiv(className = 'app-controls-block', children = list(
              htmlH3("Dataset Search:"),
              htmlP('Search for an organism or gene from the Nucleotide database,
                    to retrieve the top 10 related dataset Accession IDs and their descriptions.
                    Example search: "Basiliscus basiliscus[Organism]"
                    
                    ', style = list("margin" = 10)),
              dccInput(
                id = 'search-input',
                placeholder = "Search for a species or genes...",
                type = "text",
                value = "",
                style = list("margin" = "10px")
              ),
              htmlButton(
                id = 'search-button', n_clicks = 0, children = 'Search for Datasets'
              ),
              htmlDiv(id = "search-output", children = dashDataTable(
                id = "table",
                columns = lapply(colnames(df), 
                                 function(colName){
                                   list(
                                     id = colName,
                                     name = colName
                                   )
                                 }),
                data = df_to_list(df)
              ))
            ))
          ))
        )
      ))
    )),
    second_row
  )
)



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
    htmlH2("NCBI Explorer"),
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



# App Layout

app$layout(htmlDiv(list(
  header,
  htmlBr(),
  options
)))


# Callbacks

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
    input(id = 'genbank-sequence', property = 'data')
  ),
  
  update_sequence <- function(data) {
    return(data[["sequence"]])
  }
)



app$callback(
  output(id = 'alignment-container', property = 'children'),
  params = list(
    input(id = 'submit-button-2', property = 'n_clicks'),
    state(id = 'genbank-input', property = 'value')
  ),
  
  update_fasta <- function(n_clicks, accession_id) {
    accession_id <- gsub("\\s", "", accession_id)
    
    accession_id <- gsub('"', "", accession_id)
    
    accession_id <- strsplit(accession_id, ",")
    
    accession_id <- unlist(accession_id)
    
    if (n_clicks > 0) {
      bio_file <- read.GenBank(access.nb = accession_id)
      write.dna(bio_file, file ="data/random_fasta.FASTA", format = "fasta", append = TRUE, nbcol = 6, colsep = "", colw = 10)
      fasta_file <- toupper(read_file("C:/Users/hamma/Documents/Gene Expression/data/random_fasta.FASTA"))
      return(
        dashbioAlignmentChart(id = 'alignment-chart', data = fasta_file, height = 600, width = 950)
      )
    }
    
    else if (n_clicks < 1) {
      return(
        dashbioAlignmentChart(
          id = 'alignment-chart',
          data = read_file("data/alignment_viewer_p53_clustalo.FASTA"),
          height = 600,
          width = 950
        )
      )
    }
  }
)


app$callback(
  output(id = "search-output", property = "children"),
  params = list(
    input(id = "search-button", property = "n_clicks"),
    state(id = "search-input", property = "value")
  ),
  
  update_search_results <- function(n_clicks, search_term) {
    if (n_clicks > 0) {
      search_term <- gsub('"', "", search_term)
      
      search <- entrez_search(db = "nuccore", term = search_term, retmax=10)
      
      search_seqs <- entrez_fetch(db = "nuccore", id = search$ids, rettype = "fasta")
      
      sequence_titles <- as.list(unlist(str_extract_all(search_seqs, ">.*\\n")))
      
      sequence_dataframe <- do.call(rbind.data.frame, sequence_titles)
      
      names(sequence_dataframe)[1] <- "Dataset Accession ID and Description"
      
      return(generate_table(sequence_dataframe))
    }
  }
)


app$callback(
  output(id = "sequence-pie-chart", property = "figure"),
  params = list(
    input(id = 'submit-button', property = 'n_clicks'),
    state(id = "sequence-viewer", property = "sequence")),
  
  update_pie_chart <- function(n_clicks, sequence) {
    if (n_clicks > 0) {
      sequence_string <- as.character(sequence)
      chart_sequence <- strsplit(sequence_string, split = "")
      p <- plot_ly(as.data.frame(table(chart_sequence)), labels = ~chart_sequence, values = ~Freq, type = 'pie', hole = 0.6,
                   width = 600, height = 375) %>%
        layout(paper_bgcolor="#262B3D", title = "Nucleotide or Amino Acid Composition", margin = 20,
               titlefont = list(
                 family = "Agency FB",
                 size = 35,
                 color = '#ffffff'),
               font = list(
                 color = '#ffffff'
               )
        )
      
      
    }
    else {
      
      p <- (plot_ly(as.data.frame(table(strsplit(as.character(FASTA_DATA$`sp|Q9W678|P53_BARBU Cellular tumor antigen p53 OS=Barbus barbus GN=tp53 PE=2 SV=1`)
                                                 , split = ""))), labels = ~Var1, values = ~Freq, type = 'pie', hole = 0.6,
                    width = 600, height = 375) %>%
              layout(paper_bgcolor="#262B3D", title = "Nucleotide or Amino Acid Composition", margin = 20,
                     titlefont = list(
                       family = "Agency FB",
                       size = 35,
                       color = '#ffffff'),
                     font = list(
                       color = '#ffffff'
                     )
              )
      )
    }
    return(p)
  }
)


app$callback(
  output(id = "cpg-options", property = "options"),
  params = list(
    input(id = "submit-button-2", property = "n_clicks")
  ),
  
  select_trace <- function(n_clicks) {
    if (n_clicks > 0) {
      dna_data = readDNAStringSet("data/random_fasta.FASTA", "fasta")
      trace_names <- unique(names(dna_data))
      
      list_options <- lapply(trace_names, function(x) {
        list(label = x, value = match(x, trace_names))
      })
    }
    else {
      list_options <- list(list(label = "Default", value = 1))
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
      dna_data = readDNAStringSet("data/random_fasta.FASTA", "fasta")
      window = 100
      gc = rowSums(letterFrequencyInSlidingView(dna_data[[option]], window, c("G", "C")))/window
      g = as.data.frame(gc)
      g$index = as.numeric(rownames(g))
      m <- g[which.max(g$gc), ]
      
      p <- plot_ly(as.data.frame(g), x=~index, y=~gc, type = "scatter", mode = "line", line = list(color="#ffffff")) %>%
        layout(paper_bgcolor="#262B3D", plot_bgcolor="#262B3D", title = "CpG Island Distribution", margin = 20,
               titlefont = list(
                 family = "Agency FB",
                 size = 35,
                 color = '#ffffff'),
               font = list(
                 color = '#ffffff'
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
               )
        )
    }
    
    else {
      dna_data = readDNAStringSet("data/dnastring.FASTA", "fasta")
      window = 100
      gc = rowSums(letterFrequencyInSlidingView(dna_data[[option]], window, c("G", "C")))/window
      g = as.data.frame(gc)
      g$index = as.numeric(rownames(g))
      m <- g[which.max(g$gc), ]
      
      p <- plot_ly(as.data.frame(g), x=~index, y=~gc, type = "scatter", mode = "line", line = list(color="#ffffff")) %>%
        layout(paper_bgcolor="#262B3D", plot_bgcolor="#262B3D", title = "CpG Island Distribution", margin = 20,
               titlefont = list(
                 family = "Agency FB",
                 size = 35,
                 color = '#ffffff'),
               font = list(
                 color = '#ffffff'
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
               )
        )
    }
    return(p)
  }
)

app$callback(
  output(id = "sequence-viewer", property = "selection"),
  params = list(
    input(id = "gc_graph", property = "clickData")
  ),
  update_selection <- function(hoverdata) {
    point <- hoverdata$points[[1]]$pointNumber
    selection_fixed <- list(point, as.numeric(point+20), "green")
    str(selection_fixed)
    return(selection_fixed)
  }
)

# app$callback(
#   output(id = 'genbank-sequence', property = "data"),
#   params = list(
#     input(id ='submit-button-2', property = 'n_clicks'),
#     state(id ='cpg-options', property = 'value')
#   ),
#   update_sequence_2 <- function(n_clicks, option) {
#     if (n_clicks > 0) {
#       dna_data = readDNAStringSet("data/random_fasta.FASTA", "fasta")
#       return(as.character(dna_data[[option]]))
#     }
#   }
# )



if (appName != "") {
  
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
  
} else {
  
  app$run_server(showcase = TRUE)
  
}

