library(jsonvalidate)
library(jsonlite)
library(httr)
# to make loading gff files easier:
library(ape)
library(plyr)
library(dplyr)

EMPTY_MUT_DATA <- list(
  x = list(),
  y = list(),
  mutationGroups = list(),
  domains = list()
)

PFAM_DOM_SCHEMA <-'{
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "region": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "text": {"type": "string"},
                        "start": {"type": "string"},
                        "end": {"type": "string"},
                    },
                    "required": ["text", "start", "end"]
                }
            },
        },
        "required": ["regions"]
    }
}'


PROT_DOM_SCHEMA <- '{
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "name": {"type": "string"},
            "coord": {"type": "string"},
        },
        "required": ["name", "coord"]
    }
}'

ARR_SCHEMA <- '{
    "type": "array",
    "items": {"type": ["string", "number"]}
}'

MUT_DATA_SCHEMA <- sprintf('{
  "type": "object",
  "properties": {
    "x": %s,
    "y": %s,
    "mutationGroups": %s,
    "domains": %s,
  },
  "required": ["x"]}', 
  ARR_SCHEMA, ARR_SCHEMA, ARR_SCHEMA, PROT_DOM_SCHEMA
)

parse_mutation_data <- function(mutation_data){
  # take a json object and extract the mutation data based on the schema
  # EMPTY_MUT_DATA. Could potentially use toJSON(fromJSON(mutation_data))?
  data <- EMPTY_MUT_DATA 
  #json_validate(mutation_data, MUT_DATA_SCHEMA)
  mutation_data <- fromJSON(mutation_data)
  for (k in names(data)){
    data[k] <- mutation_data[k]
  }
  data
}

load_mutation_data <- function(json_fname = NULL){
  #take a json object and extract the mutation data based one the schema EMPTY_MUT_DATA
  if (!is.null(json_fname)){
    mutationData <- readLines(json_fname)
    return(parse_mutation_data(mutationData))
  } else {
    return(EMPTY_MUT_DATA)
  }
}

query_fields <- list(
  "accession",
  "active",
  "annotation",
  "author",
  "cdantigen",
  "citation",
  "cluster",
  "count",
  "created",
  "database",
  "ec",
  "evidence",
  "existence",
  "family",
  "fragment",
  "gene",
  "gene_exact",
  "goa",
  "host",
  "id",
  "inn",
  "interactor",
  "keyword",
  "length",
  "lineage",
  "mass",
  "method",
  "mnemonic",
  "modified",
  "name",
  "organelle",
  "organism",
  "plasmid",
  "proteome",
  "proteomecomponent",
  "replaces",
  "reviewed",
  "scope",
  "sequence",
  "sequence_modified",
  "source",
  "strain",
  "taxonomy",
  "tissue",
  "web"
)

query_parameters <- list(
  format = list(
    "html",
    "tab",
    "xls",
    "fasta",
    "gff",
    "txt",
    "xml",
    "rdf",
    "list",
    "rss"
  ),
  columns = list(
		"citation",
		"clusters",
		"comments",
		"domains",
		"domain",
		"ec",
		"id",
		"entry name",
		"existence",
		"families",
		"features",
		"genes",
		"go",
		"go-id",
		"interactor",
		"keywords",
		"last-modified",
		"length",
		"organism",
		"organism-id",
		"pathway",
		"protein names",
		"reviewed",
		"sequence",
		"3d",
		"version",
		"virus hosts"
  ),
	sort = list("score"),
	include = list("yes", "no"),
	compress = list("yes", "no"),
	limit = "int",
	offset = "int"
)

validate_query_parameters <- function(parameters = NULL){
  if (is.null(parameters)){
    parameters = list()
  }
  validated_parameters <- list()
  for (param in names(parameters)){
    if (param %in% names(query_parameters)){
      if (param %in% list("limit", "offset")){
        if (is.numeric(parameters[[param]])){
          validated_parameters[[param]] <- as.character(parameters[[param]])
        }
      } else if (param == "format"){
        if (parameters[[param]] %in% query_parameters[[param]]){
          validated_parameters[[param]] <- parameters[[param]] 
        }
      } else if (param == "columns"){
        column_entry <- unlist(
          strsplit(
            gsub("+", "", parameters[[param]], fixed = TRUE), ","
          )
        )
        set_a <- unlist(unique(query_parameters[[param]]))
        set_b <- unlist(unique(column_entry))
        validated_items <- intersect(set_a, set_b)
        validated_column_entry <- list()
        for (i in 1:length(column_entry)){
          if (column_entry[[i]] %in% validated_items){
            validated_column_entry[[i]] <- column_entry[[i]]
          }
        }
        validated_parameters[[param]] <- paste(
          unlist(validated_column_entry), collapse = ","
        )
      } else if (param %in% list("include", "compress", "sort")){
        if (parameters[[param]] %in% query_parameters[[param]]){
          validated_parameters[[param]]  <- parameters[[param]]
        }
      }
    }
  }
  validated_parameters
}

# Functions apdapted from the UniprotQueryBuilder class defined in dash_bio_utils:
build_query <- function(query, fields = NULL, parameters = NULL){
  base_url <- "https://www.uniprot.org/uniprot/"
  base_query <- "?query=%s"
  field_separator <- "+AND+"
  parameter_separator <- "&"
  if (is.null(fields)){
    fields <- list()
  }
  if (is.null(parameters)){
    parameters <- list()
  }
  URL <- paste0(base_url, sprintf(base_query, query))
  for (fieldname in names(fields)){
    if (fieldname %in% query_fields){
      if (is.character(fields[[fieldname]])){
        URL <- paste0(URL, field_separator, fieldname, ":", fields[[fieldname]])
      } else {
        stop(
          sprintf("The value of the field %s is not a string format", fieldname)
        )
      }
    }
  }
  validated_parameters <- validate_query_parameters(parameters)
  for (param in names(validated_parameters)){
    URL <- paste0(
      URL, parameter_separator, param, "=", validated_parameters[[param]]
    )
  }
  URL
}

query_into_dataframe <- function(
  query, 
  fields = NULL, 
  parameters = NULL, 
  names = NULL
  ){
  target_url <- build_query(query, fields = fields, parameters = parameters)
  col_id <- "columns"
  col_names <- NULL
  if (is.null(names)){
    db <- read.csv(URLencode(target_url), sep = "\t")
  } else {
    db <- ape::read.gff(file = URLencode(target_url))
    colnames(db) <- names
  }
  db
}

pfam_domain_parser <- function(accession){
  URL <- sprintf("http://pfam.xfam.org/protein/%s/graphic", accession)  
  r <- GET(URL)
  jsonData <- content(r, "parsed")
  toJSON(jsonData)
}

parse_protein_domains_data <- function(domain_data){
  region_key <- "regions"
  region_name_key <- "text"
  region_start_key <- "start"
  region_stop_key <- "end"
  formatted_data <- list()

  if (json_validate(domain_data, PFAM_DOM_SCHEMA)){
    regionlist <- fromJSON(
      domain_data, simplifyVector = FALSE
    )[[1]][[region_key]]
    for (i in 1:length(regionlist)){
      formatted_data[[i]] <- list(
        name = unlist(regionlist[[i]][[region_name_key]]),
        coord = sprintf(
          "%s-%s",
          regionlist[[i]][[region_start_key]],
          regionlist[[i]][[region_stop_key]]
        )
      )
    }
  } else if (json_validate(domain_data, PROT_DOM_SCHEMA)){
    formatted_data <- domain_data
  }
  formatted_data
}

load_protein_domains <- function(accession){
  domain_data <- pfam_domain_parser(accession)
  parse_protein_domains_data(domain_data)
}

parse_mutations_uniprot_data <- function(
  gff_data, 
  start = "start", 
  stop = "end", 
  mut_types_to_skip = NULL
  ){
  if (is.null(mut_types_to_skip)){
    mut_types_to_skip <- list(
      "Chain",
      "Region"
    )
  }
  if (!"Chain" %in% mut_types_to_skip){
    mut_types_to_skip <- c(mut_types_to_skip, list("Chain"))
  }
  # Selects the various mutations types in the dataset, except types contained in the above list
  mut_types <- as.character(
    unique(gff_data[!gff_data$mut %in% mut_types_to_skip, "mut"])
  )
  x <- list()
  y <- list()
  mutationgroups <- list()

  for (mut in mut_types){
    data_coord <- gff_data[gff_data$mut == mut, c(start, stop)]
    # split between single and multi-site coordinates
    single_sites <- data_coord[data_coord[,start] == data_coord[,stop],]
    multi_sites <- data_coord[data_coord[,start] != data_coord[,stop],]
    multi_sites[,start] <- paste(
      as.character(multi_sites[,start]),
      as.character(multi_sites[,stop]),
      sep = "-"
    )
    sorted_data <- plyr::count(c(single_sites[,start], multi_sites[,start])) %>% 
      mutate(mut = mut)

    x <- c(x, as.character(sorted_data[,"x"])) 
    y <- c(y, as.character(sorted_data[,"freq"]))
    mutationgroups <- c(mutationgroups, sorted_data[, "mut"])
  }

  # order the results by frequency of occurrence
  # There is an issue with the ordering of the mutation groups:
  # For the time being, the data look correct, but the mutation groups they are 
  # attached to are sometimes re-ordered in an unexpected manner
  # For an example, source DDX3X from the UniProt database 
  # The "mutagenesis" and "helix" labels are somehow reversed...
  order_df <- as.data.frame(
    sort(
      table(unlist(mutationgroups)), 
      decreasing = TRUE
    )
  )
  order_df <- merge(
    data.frame(
      mut = unlist(mutationgroups),
      x = unlist(x),
      y = unlist(y)
    ),
    order_df,
    by.x = "mut", by.y = "Var1" 
  )
  order_df <- order_df[order(order_df[, "Freq"], decreasing = TRUE),]
  
  formatted_data = list(
    "x" = as.character(order_df$x),
    "y" = as.numeric(order_df$y),
    "mutationGroups" = as.character(order_df$mut),
    domains = list()
  )
  formatted_data
}

parse_mutation_upload_file <- function(contents, fname){
  data <- EMPTY_MUT_DATA
  if (endsWith(fname, ".json"))
    content_string <- unlist(strsplit(contents, ","))[2]
    decoded <- base64_dec(content_string)
    data <- fromJSON(rawToChar(decoded))
  data
}

parse_domain_upload_file <- function(contents, fname){
  data <- list()
  if (endsWith(fname, ".json")){
    content_string <- unlist(strsplit(contents, ","))[2]
    decoded <- base64_dec(content_string)
    data <- fromJSON(rawToChar(decoded))
  }
}

