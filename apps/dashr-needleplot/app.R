library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashBio)
library(jsonlite)
source("utils/utils.R")

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}


DATAPATH <- "sample_data/needle_" 

# Data used for the default demo plot
DATA_URL <- "https://raw.githubusercontent.com/bbglab/muts-needle-plot/master/snippets/data/"

DEMO_DATA <- list(
  list("mutData" = "TP53.json", label = "TP53"),
  list("mutData" = "ACVR1.json", label = "ACVR1"),
  list("mutData" = "SMARCA4.json", label = "SMARCA4"),
  list("mutData" = "ENTS00000557334.json", label = "ENTS00000557334"),
  list("mutData" = "PIK3CA.json", label = "PIK3CA"),
  list("mutData" = "ATRX.json", label = "ATRX")
)

# Values of a the load dataset dropdown
DEMO_KEY <- "demo"
FILE_KEY <- "file"
DATABASE_KEY <- "db"

# Values of the data download dropdown
FULL_KEY <- "needleplot"
MUT_KEY <- "mutations"
DOM_KEY <- "protein_domains"

# Value of the cheklist to load the protein domain individually
INDIV_DOMS_KEY <- "indivual_domain_loading"
UNIPROT_DOMS_KEY <- "uniprot_domain_loading"

DB_LAST_QUERY_KEY <- sprintf("%s-last-query", DATABASE_KEY)

HEAD_COLORS <- list(
  "#e41a1c",
  "#377eb8",
  "#4daf4a",
  "#984ea3",
  "#ff7f00",
  "#ffff33",
  "#a65628",
  "#f781bf",
  "#999999"
)

HEAD_SYMBOLS <- list(
	"circle",
	"square",
	"triangle-up",
	"diamond"
)

STEM_COLOR <- list(
  "grey",
	"black",
  "white"
)

DOMAIN_COLORS <- list(
  "#8dd3c7",
  "#ffffb3",
  "#bebada",
  "#fb8072",
  "#80b1d3",
  "#fdb462",
  "#b3de69",
  "#fccde5",
  "#d9d9d9",
  "#bc80bd",
  "#ccebc5",
  "#ffed6f"
)

description <- function(){
  paste(
    "Display gene mutation of the genome thanks to this needle plot.",
    "Also known under the lollipop plot name."
  )
}

header_colors <- function(){
  list(
    bg_color = "#0D1A51",
    font_color = "#FFFFFF"
  )
}

app <- Dash$new()

app$layout(
  htmlDiv(
    children = list(
      htmlDiv(
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
          htmlH2("Needle Plot"),
          htmlA(
            id = "gh-link",
            children = list("View on GitHub"),
            href = "https://github.com/plotly/dash-sample-apps/tree/master/apps/dashr-needleplot",
            style = list(color = "white", border = "solid 1px white")
          ),
          htmlImg(
            src = "assets/GitHub-Mark-Light-64px.png"
          )
        )
      ),
      htmlBr(),
      htmlDiv(
        id = "needleplot-body",
        className = "app-body",
        children = list(
          dccLoading(
            className = "dashbio-loading",
            children = htmlDiv(
              id = "needleplot-wrapper",
              children = dashbioNeedlePlot(
                id = "needle-plot",
                rangeSlider = TRUE
              ) 
            )
          ),
          htmlDiv(
            id = "needleplot-control-tabs",
            className = "control-tabs",
            children = list(
              dccTabs(
                id = "needleplot-tabs",
                value = "what-is",
                children = list(
                  dccTab(
                    label = "About",
                    value = "what-is",
                    children = htmlDiv(
                      className = "control-tab",
                      children = list(
                          htmlH4(
                          className = "what-is", 
                          children = "What is Needle Plot?"
                        ),
                        htmlP(
                          paste0(
                            "Needle Plot allows you to display mutations in",
                            "a genome. Due to its similarity to both a barplot",
                            "and a scatter plot, it can be used to plot",
                            "datasets that have too many mutations for a",
                            "barplot to be meaningful."
                          )
                        ),
                        htmlP(
                          paste0(
                            "In the \"Data\" tab, you can choose from preloaded",
                            "datasets, as well as upload your own. You can",
                            "additionally search the UniProt database for",
                            "data to plot. If you wish to save the data that are",
                            "plotted, you can choose to download all of it, or",
                            "just the data corresponding to mutations or domains."
                          ) 
                        ),
                        htmlP(
                          paste0(
                            "In the \"Graph\" tab, you can change the aesthetics of",
                            "the data points by customizing colors, marker shapes,",
                            "and more."
                          ) 
                        )
                      )
                    )
                  ),
                  dccTab(
                    label = "Data",
                    value = "datasets",
										children = htmlDiv(
                      className = "control-tab",
                      children = list(
                        htmlDiv(
                          className = "app-controls-block",
                          title = paste(
                            "\"Demo dataset\" choice will allow you to play",
                            "with the options.\n",
                            "\"UniProt dataset\" choice will retrieve protein",
                            "domain as well as mutation data from UniProt",
                            "database.\n\"Upload dataset\" choice will let",
                            "you choose your own mutation data with the",
                            "option to load the protein domains from pfam",
                            "database."
                          ),
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Dataset source:"
                            ),
                            dccDropdown(
                              id = "needle-dataset-select-dropdown",
                              options = list(
                                list(
                                  label = "Demo dataset", 
                                  value = DEMO_KEY
                                ),
                                list(
                                  label = "Upload dataset", 
                                  value = FILE_KEY
                                ),
                                list(
                                  label = "UniProt dataset", 
                                  value = DATABASE_KEY
                                )
                              ),
                              value = DEMO_KEY
                            )
                          )
                        ),
                        htmlHr(),
                        htmlDiv(
                          id = sprintf("needle-%s-div", DEMO_KEY),
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Dataset:"
                            ),
                            dccDropdown(
                              id = "needle-dataset-dropdown",
                              options = lapply(
                                1:length(DEMO_DATA),
                                function(i){
                                  list(
                                    label = DEMO_DATA[[i]][["label"]],
                                    value = i
                                  )
                                }
                              ),
                              value = 1
                            )
                          )
                        ),
                        htmlDiv(
                          id = "needle-protein-domains-select-div",
                          className = "app-controls-block",
                          title = paste(
                            "Check this box to enable loading of",
                            "mutation data such as the protein coordinate",
                            "(x), mutation number (y) and mutation type",
                            "(mutationGroups), individually from the protein",
                            "domains"
                          ),
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Load protein domains:"
                            ),
                            dccChecklist(
                              id = "needle-protein-domains-select-checklist",
                              className = "needle-checklist",
                              options = list(
                                list(
                                  label = "Individually", 
                                  value = INDIV_DOMS_KEY
                                ),
                                list(
                                  label = "From UniProt only", 
                                  value = UNIPROT_DOMS_KEY
                                )
                              ),
                              value = list()
                            )
                          )
                        ),
                        htmlDiv(
                          id = sprintf("needle-%s-div", DATABASE_KEY),
                          className = "app-controls-block",
                          children = list(
                            htmlH5("Search UniProt"),
                            htmlDiv(
                              title = paste(
                                "Enter the UniProt accession key",
                                "of the gene you want to display.\n",
                                "More information on https://www.uniprot.org/"
                              ),
                              children = list(
                                dccInput(
                                  id = "needle-sequence-input",
                                  value = "",
                                  type = "text",
                                  placeholder = "TP53, DDX3X, SMARCA4, ..."
                                ),
                                htmlButton(
                                  id = "needle-search-sequence-button",
                                  children = "submit",
                                  n_clicks = 0,
                                  n_clicks_timestamp = 0
                                )
                              )
                            ),
                            htmlDiv(id = "needle-uniprot-div")
                          )
                        ),
                        htmlDiv(
                          id = sprintf("needle-%s-div", FILE_KEY),
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              id = "needle-mutdata-file-div",
                              title = paste(
                                "Mutation data files are JSON files containing",
                                "the following fields :\n",
                                "- \"x\" (protein coordinate of the mutation); \n",
                                "- \"y\" (number of recorded mutations); \n",
                                "- \"mutationGroups\" (type of mutations); \n",
                                "- \"domains\" (protein domains).\n",
                                "\"x\", \"y\", and \"mutationGroups\" are arrays, they",
                                "must have the same length.",
                                "\"x\" is required. \"domains\" is an array of",
                                "JSON objects with required fields \"name\"",
                                "and \"coord\"; \"name\" can be any string (think",
                                "of it as a label), whereas \"coord\" should be",
                                "a string formatted like",
                                "\"<start_coord>-<stop_coord>\" where integer",
                                "<start_coord> is less than integer <stop_coord>",
                                "(e.g., \"23-34\"), giving the domain in protein",
                                "coordinates."
                              ),
                              children = list(
                                htmlH5("Upload mutation data JSON file"),
                                dccUpload(
                                  id = "needle-mutdata-file-upload",
                                  className = "control-upload",
                                  children = htmlDiv(
                                    list(
                                      "Drag and drop or ",
                                      htmlA("select files")
                                    )
                                  )
                                ),
                                htmlDiv(
                                  id = "needle-mutdata-file-info-div"
                                )
                              )
                            ),
                            htmlDiv(
                              id = "needle-domain-file-div",
                              title = paste(
                                "Protein data files accepted here can be of",
                                "two types: \n",
                                "- an array of JSON objects with the required",
                                "fields, i.e., \"name\" (any string) and \"coord\"",
                                "(a string specifying domains in protein",
                                "coordinates (e.g., \"23-34\");",
                                "- a JSON file with the same structure as: ",
                                "http://pfam.xfam.org/protein/P04637/graphic."
                              ),
                              children = list(
                                htmlH5("Upload protein data JSON file"),
                                dccUpload(
                                  id = "needle-domains-file-upload",
                                  className = "needle-upload",
                                  children = htmlDiv(
                                    list(
                                      "Drag and drop or ",
                                      htmlA("select files")
                                    )
                                  )
                                ),
                                htmlDiv(
                                  id = "needle-domains-file-info-div"
                                )
                              )
                            ),
                            htmlDiv(
                              id = "needle-domain-query-info-div"
                            )
                          )
                        ),
                        htmlHr(),
                        htmlBr(),
                        htmlDiv(
                          id = "needle-download-data-div",
                          children = list(
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Download data:"
                                ),
                                dccDropdown(
                                  id = "needle-download-data-dropdown",
                                  options = list(
                                    list(label = "All", value = FULL_KEY),
                                    list(label = "Mutations", value = MUT_KEY),
                                    list(label = "Domains", value = DOM_KEY)
                                  ),
                                  value = FULL_KEY
                                )
                              )
                            ),
                            htmlA(
                              id = "needle-download-data-button-link",
                              children = htmlButton(
                                id = "needle-download-data-button",
                                className = "control-download",
                                children = "Download graph data",
                                n_clicks = 0,
                                n_clicks_timestamp = 0
                              ),
                              href = "",
                              download = ""
                            )
                          ) 
                        )
                      )
                    )
                  ),
                  dccTab(
                    label = "Graph",
                    value = "graph",
                    children = htmlDiv(
                      className = "control-tab",
                      children = list(
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Stem thickness"
                            ),
                            dccSlider(
                              id = "needle-stem-thick-input",
                              value = 2,
                              min = 1,
                              max = 10,
                              step = 1
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Needle head size"
                            ),
                            dccSlider(
                              id = "needle-head-size-input",
                              value = 7,
                              min = 1,
                              max = 10,
                              step = 1
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Stem color()"
                            ),
                            dccDropdown(
                              id = "needle-stem-color-dropdown",
                              options = lapply(
                                STEM_COLOR,
                                function(col){
                                  list(label = col, value = col)
                                }
                              ),
                              value = STEM_COLOR[[1]]
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Head color(s)"
                            ),
                            dccDropdown(
                              id = "needle-head-color-dropdown",
                              options = lapply(
                                HEAD_COLORS,
                                function(col){
                                  list(label = col, value = col)
                                }
                              ),
                              value = HEAD_COLORS[1:4],
                              multi = TRUE
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Head symbol(s)"
                            ),
                            dccDropdown(
                              id = "needle-head-symbol-dropdown",
                              options = lapply(
                                HEAD_SYMBOLS,
                                function(sym){
                                  list(label = sym, value = sym)
                                }
                              ),
                              value = HEAD_SYMBOLS[1:4],
                              multi = TRUE
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Constant height needles"
                            ),
                            dccRadioItems(
                              id = "needle-stem-height-radioitems",
                              className = "needle-radio",
                              options = list(
                                list(label = "On", value = TRUE),
                                list(label = "Off", value = FALSE)
                              ),
                              value = FALSE
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Rangeslider display"
                            ),
                            dccRadioItems(
                              id = "needle-rangeslider-radioitems",
                              className = "needle-radio",
                              options = list(
                                list(label = "On", value = TRUE),
                                list(label = "Off", value = FALSE)
                              ),
                              value = TRUE 
                            )
                          )
                        ),
                        htmlDiv(
                          className = "app-controls-block",
                          children = list(
                            htmlDiv(
                              className = "app-controls-name",
                              children = "Small domains color(s)"
                            ),
                            dccDropdown(
                              id = "needle-domains-color-dropdown",
                              options = lapply(
                                DOMAIN_COLORS,
                                function(col){
                                  list(label = col, value = col)
                                }
                              ),
                              value = DOMAIN_COLORS[1:4],
                              multi = TRUE
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
        )
      ),
    dccStore(id = "needle-store")
    )
  )
)

app$callback(
  output(sprintf("needle-%s-div", DATABASE_KEY), "style"),
  list(
    input("needle-dataset-select-dropdown", "value"),
    state(sprintf("needle-%s-div", DATABASE_KEY), "style")
  ),
  function(load_choice, div_style){
    # updates what the user can use to load data to the graph,
    if (is.null(unlist(div_style))){
      div_style = list(display = "none")
    }
    if (load_choice == DATABASE_KEY){
      div_style[["display"]] <- "inherit"
    } else {
      div_style[["display"]] <- "none"
    }
    div_style
  }
)

app$callback(
  output(sprintf("needle-%s-div", DEMO_KEY), "style"),
  list(
    input("needle-dataset-select-dropdown", "value"),
    state(sprintf("needle-%s-div", DEMO_KEY), "style")
  ),
  function(load_choice, div_style){
    # updates what the user can use to load data to the graph,
    if (is.null(unlist(div_style))){
      div_style = list(display = "none")
    }
    if (load_choice == DEMO_KEY){
      div_style[["display"]] <- "block"
    } else {
      div_style[["display"]] <- "none"
    }
    div_style
  }
)

app$callback(
  output(sprintf("needle-%s-div", FILE_KEY), "style"),
  list(
    input("needle-dataset-select-dropdown", "value"),
    state(sprintf("needle-%s-div", FILE_KEY), "style")
  ),
  function(load_choice, div_style){
    # updates what the user can use to load data to the graph,
    if (is.null(unlist(div_style))){
      div_style = list(display = "none")
    }
    if (load_choice == FILE_KEY){
      div_style[["display"]] <- "inherit"
    } else {
      div_style[["display"]] <- "none"
    }
    div_style
  }
)

app$callback(
  output("needle-uniprot-div", "style"),
  list(
    input("needle-dataset-select-dropdown", "value"),
    state("needle-uniprot-div", "style")
  ),
  function(load_choice, div_style){
    # updates what the user can use to load data to the graph,
    if (is.null(unlist(div_style))){
      div_style = list(display = "none")
    }
    if (load_choice == DATABASE_KEY){
      div_style[["display"]] <- "inherit"
    } else {
      div_style[["display"]] <- "none"
    }
    div_style
  }
)

app$callback(
  output("needle-uniprot-div", "children"),
  list(
    input("needle-store", "data"),
    state("needle-dataset-select-dropdown", "value")
  ),
  function(stored_data, load_choice){
    # displays information about the query to the UniProt database 
    div <- htmlDiv()
    if (load_choice == DATABASE_KEY){
      if (stored_data[["info"]][["is_same_key"]] == TRUE){
        title <- "Query"
      } else {
        title <- "Last query"
      }
      #last_query <- stored_data[["info"]][[load_choice]]
      last_query <- stored_data[["info"]][[DATABASE_KEY]]
      if (last_query[[1]] != ""){
        div <- htmlDiv(
          list(
            htmlH5(title),
            htmlP(last_query)
          )
        )
      }
    }
    div
  }
)

app$callback(
  output("needle-sequence-input", "value"),
  list(
    input("needle-store", "data"),
    state("needle-dataset-select-dropdown", "value")
  ),
  function(stored_data, load_choice){
    # resets the last query if the user changed dataset loading option
    if (load_choice == DATABASE_KEY){
      answer <- stored_data[["info"]][[DB_LAST_QUERY_KEY]]
    } else {
      answer = ""
    }
    answer
  }
)

app$callback(
  output("needle-domain-file-div", "style"),
  list(
    input("needle-protein-domains-select-checklist", "value"),
    state("needle-domain-file-div", "value")
  ),
  function(domains_opt, div_style){
    # toggles the view of the domain upload div
    if (is.null(unlist(div_style))){
      div_style = list(display = "none")
    }
    if (
      (INDIV_DOMS_KEY %in% domains_opt) & 
        (!UNIPROT_DOMS_KEY %in% domains_opt)){
      div_style[["display"]] <- "inherit"
    } else {
      div_style[["display"]] <- "none"
    }
    div_style
  }
)

app$callback(
  output("needle-domain-query-info-div", "style"),
  list(
    input("needle-protein-domains-select-checklist", "value"),
    state("needle-domain-query-info-div", "style"),
    state("needle-dataset-select-dropdown", "value")
  ),
  function(domains_opt, div_style, load_choice){
    # toggles the view of the domain-query-info div which
    # displays information from the UniProt query
    if (is.null(unlist(div_style))){
      div_style = list(display = "none")
    }
    div_style[["display"]] <- "none"
    if (load_choice == FILE_KEY){
      if (UNIPROT_DOMS_KEY %in% domains_opt){
        div_style[["display"]] <- "inherit"
      }
    }
    div_style
  }
)

app$callback(
  output("needle-protein-domains-select-div", "style"),
  list(
    input("needle-dataset-select-dropdown", "value"),
    state("needle-protein-domains-select-div", "style")
  ),
  function(load_choice, div_style){
    # toggles the view of the protein domain select div which displays
    # the ability to load protein domains independently from mutation data
    if (is.null(unlist(div_style))){
      div_style = list(display = "none")
    }
    if (load_choice == FILE_KEY){
      div_style[["display"]] <- "inherit"
    } else {
      div_style[["display"]] <- "none"
    }
    div_style
  }
)


# UPLOAD RELATED CALLBACKS========
app$callback(
  output("needle-mutdata-file-upload", "contents"),
  list(
    input("needle-dataset-select-dropdown", "value"),
    state("needle-mutdata-file-upload", "contents")
  ),
  function(load_choice, mut_contents){
    # reset the content of the mutation upload
    if (load_choice == DEMO_KEY){
      answer <- NULL
    } else {
      answer <- mut_contents
    }
    answer
  }
)

app$callback(
  output("needle-domains-file-upload", "contents"),
  list(
    input("needle-dataset-select-dropdown", "value"),
    state("needle-domains-file-upload", "contents")
  ),
  function(load_choice, dom_contents){
    # reset the content of the protein domains upload
    if (load_choice == DEMO_KEY){
      answer <- NULL
    } else {
      answer <- dom_contents
    }
    answer
  }
)

app$callback(
  output("needle-mutdata-file-info-div", "children"),
  list(
    input("needle-mutdata-file-upload", "contents"),
    input("needle-mutdata-file-upload", "filename")
  ),
  function(mut_contents, mut_fname){
    # display the info about the source of the protein mutations data
    if (!is.null(unlist(mut_contents))){
      answer <- sprintf("Loaded from : %s", mut_fname)
    } else {
      answer <- list()
    }
    answer
  }
)

app$callback(
  output("needle-domains-file-info-div", "children"),
  list(
    input("needle-domains-file-upload", "contents"),
    input("needle-domains-file-upload", "filename")
  ),
  function(dom_contents, dom_fname){
    # display the info about the source of the protein domains data
    if (!is.null(unlist(dom_contents))){
      answer <- sprintf("Loaded from : %s", dom_fname)
    } else {
      answer <- list()
    }
    answer
  }
)

app$callback(
  output("needle-domain-query-info-div", "children"),
  list(
    input("needle-store", "data"),
    state("needle-protein-domains-select-checklist", "values")
  ),
  function(stored_data, domains_opt){
    # display the info about the source of the protein domains data
    div <- list()
    accession <- stored_data[[INDIV_DOMS_KEY]][["accession"]]
    if ((UNIPROT_DOMS_KEY %in% domains_opt) & (!is.null(accession))){
      div <- htmlDiv(
        htmlH5("Protein domains loaded from "),
        htmlA(
          sprintf("http://pfam.xfam.org/protein/%s/graphic", accession),
          href = sprintf("http://pfam.xfam.org/protein/%s/graphic", accession)
        )
      ) 
    }
    div
  }
)

# DATA RELATED CALLBACKS==========
app$callback(
  output("needle-download-data-button-link", "download"),
  list(
    input("needle-plot", "mutationData"),
    input("needle-download-data-dropdown", "value")
  ),
  function(x_, dl_data_choice){
    # change the file name of downloadable data based on user choice
    sprintf("%s_data.json", dl_data_choice)
  }
)

app$callback(
  output("needle-download-data-button-link", "href"),
  list(
    input("needle-plot", "mutationData"),
    input("needle-download-data-dropdown", "value"),
    state("needle-store", "data")
  ),
  function(x_, dl_data_choice, stored_data){
    # change the link to the downloadable data
    fname <- sprintf("%s_data.json", dl_data_choice)
    fpath <- sprintf("%s%s", stored_data[["info"]][["dl_data_path"]], fname)
    fpath
  }
)

app$callback(
  output("needle-download-data-div", "style"),
  list(
    input("needle-plot", "mutationData"),
    state("needle-download-data-div", "style")
  ),
  function(plotted_data, div_style){
    # disable download-data-div if there is no data plotted
    if (is.null(unlist(div_style))){
      div_style <- list(display = "block")
    }
    div_style <- list(display = "block")
    if ((is.null(plotted_data[["domains"]])) & (is.null(plotted_data[["x"]]))){
      div_style[["display"]] <- "none"
    }
    div_style
  }
)

app$callback(
  output("needle-plot", "mutationData"),
  list(
    input("needle-store", "data"),
    state("needle-download-data-dropdown", "value")
  ),
  function(stored_data, dl_data_choice){
    # reloads the dataset of thegraph if the stored data changed
    fname <- sprintf("%s_data.json", dl_data_choice)
    write_json(
      stored_data[["plot"]], 
      path = sprintf("%s%s", stored_data[["info"]][["dl_data_path"]], fname),
      auto_unbox = TRUE
    )
    stored_data[["plot"]]
  }
)


app$callback(
  output("needle-store", "data"),
  list(
    input("needle-search-sequence-button", "n_clicks"),
    input("needle-dataset-select-dropdown", "value"),
    input("needle-dataset-dropdown", "value"),
    input("needle-protein-domains-select-checklist", "value"),
    input("needle-mutdata-file-upload", "contents"),
    input("needle-domains-file-upload", "contents"),
    state("needle-mutdata-file-upload", "filename"),
    state("needle-domains-file-upload", "filename"),
    state("needle-sequence-input", "value"),
    state("needle-store", "data")
  ),
  function(
    x_,
    load_choice,
    demo_choice,
    domains_opt,
    mut_contents,
    dom_contents,
    mut_fname,
    dom_fname,
    query,
    stored_data
  ){
    if (is.null(unlist(stored_data))){
      stored_data <- list(
        plot = EMPTY_MUT_DATA,
        info = list(
          previous_key = load_choice,
          is_same_key = FALSE,
          previous_choice = "",
          dl_data_path = DATAPATH
        )
      )
      stored_data[["info"]][[DB_LAST_QUERY_KEY]] <- ""
      stored_data[["info"]][[DATABASE_KEY]] <- ""
      stored_data[[INDIV_DOMS_KEY]] <- list(
        domains = list(),
        accession = ""
      )
    }
    if (load_choice == DEMO_KEY){
      # loads datasets from a github repo for demo purposes
      fpath <- stored_data[["info"]][["dl_data_path"]]
      fname <- DEMO_DATA[[demo_choice]][["mutData"]]
      stored_data[["plot"]] <- load_mutation_data(
        sprintf("%s%s", fpath, fname)
      )
    }
    if (load_choice == DATABASE_KEY){
      # loading uniprot data set
      # Performs a simple query in the UniProt database to find an
      # accession number
      if (query != ""){
        gene_search <- query_into_dataframe(
          query, 
          fields = list(
            revieved = "yes", 
            database = "pfam"
          ),
          parameters = list(
            limit = 1, 
            columns = "id,entry name,length,genes,organism",
            sort = "score",
            format = "tab"
          )
        )
        stored_data[["info"]][[DATABASE_KEY]] <- capture.output(
          print(gene_search)
        ) 
        
        accession <- gene_search[1, "Entry"]
        domains <- load_protein_domains(accession = accession)

        stored_data[[INDIV_DOMS_KEY]] <- list(
          "domains" = domains,
          "accession" = accession
        )
        gff_data <- query_into_dataframe(
          query,
          fields = list(
            revieved = "yes",
            database = "pfam",
            accession = as.character(accession)
          ),
          parameters = list(
            format = "gff"
          ),
          names = c(
            "name", "db", "mut", "start", "end", "x1", "x2", "x3", "note"
          )
        )
        formatted_data <- parse_mutations_uniprot_data(gff_data = gff_data)
        
        stored_data[["plot"]][["x"]] <- formatted_data[["x"]]
        stored_data[["plot"]][["y"]] <- formatted_data[["y"]]
        stored_data[["plot"]][["mutationGroups"]] <- 
          formatted_data[["mutationGroups"]]
        stored_data[["plot"]][["domains"]] <- domains
      }
    }
  if (load_choice == FILE_KEY){
    stored_data[["plot"]] <- parse_mutation_upload_file(mut_contents, mut_fname) 
    if (INDIV_DOMS_KEY %in% domains_opt){
      stored_data[["plot"]][["domains"]] <- list()
      if (!UNIPROT_DOMS_KEY %in% domains_opt){
        if (!is.null(unlist(dom_contents))){
          stored_data[["plot"]][["domains"]] <- parse_domain_upload_file(
            dom_contents, dom_fname
          )
        }
      } else {
        if (INDIV_DOMS_KEY %in% stored_data){
          stored_data[["plot"]][["domains"]] <- 
            stored_data[[INDIV_DOMS_KEY]][["domains"]]
        }
      }
    }
  }
  if (load_choice != stored_data[["info"]][["previous_key"]]){
    stored_data[["info"]][["is_same_key"]] <- FALSE
  } else {
    stored_data[["info"]][["is_same_key"]] <- TRUE
  }
  stored_data
  }
)

# GRAPH OPTIONS CALLBACKS=========
app$callback(
  output("needle-plot", "rangeSlider"),
  list(input("needle-rangeslider-radioitems", "value")),
  function(val){
    return(val)
  }
)

app$callback(
  output("needle-plot", "needleStyle"),
  list(
    input("needle-stem-height-radioitems", "value"),
    input("needle-stem-thick-input", "value"),
    input("needle-stem-color-dropdown", "value"),
    input("needle-head-size-input", "value"),
    input("needle-head-color-dropdown", "value"),
    input("needle-head-symbol-dropdown", "value"),
    state("needle-plot", "needleStyle")
  ),
  function(
    const_height,
    stem_thick,
    stem_color,
    head_size,
    head_colors,
    head_symbols,
    needle_sty){
    if (is.null(unlist(needle_sty))){
      needle_sty <- list()
    }
    needle_sty[["stemConstHeight"]] <- const_height
    needle_sty[["stemThickness"]] <- stem_thick 
    needle_sty[["stemColor"]] <- stem_color 
    needle_sty[["headSize"]] <- head_size 
    needle_sty[["headColor"]] <- head_colors 
    needle_sty[["headSymbol"]] <- head_symbols 
    needle_sty
  }
)

app$callback(
  output("needle-plot", "domainStyle"),
  list(
    input("needle-domains-color-dropdown", "value"),
    state("needle-plot", "domainStyle")
  ),
  function(small_domains_colors, domains_sty){
    if (is.null(unlist(domains_sty))){
      domains_sty <- list()
    }
    domains_sty[["domainColor"]] <- small_domains_colors
    domains_sty
  }
)

if (appName != "") {
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050)) 
} else {
  app$run_server()
}
