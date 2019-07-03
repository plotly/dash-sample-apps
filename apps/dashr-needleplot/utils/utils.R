library(jsonvalidate)
library(jsonlite)

EMPTY_MUT_DATA <- list(
  x = list(),
  y = list(),
  mutationGroups = list(),
  domains = list()
)

PFAM_DOM_SCHEMAr<-'{
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


#mutationData <- readLines("../sample_data/needle_needleplot_data.json")
#mutationData <- readLines("../sample_data/needle_ATRX.json")

#json_validate(mutationData, MUT_DATA_SCHEMA, error = TRUE)
#json_validate(mutationData, PROT_DOM_SCHEMA, error = TRUE)

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

}

