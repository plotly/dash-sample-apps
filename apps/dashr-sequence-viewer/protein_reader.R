library(Biostrings)
library(stringr)

DATABASES = list(
  'gb'= list('Accession', 'Locus'),
  'emb'= list('Accession', 'Locus'),
  'dbj'= list('Accession', 'Locus'),
  'pir'= list('Entry'),
  'prf'= list('Name'),
  'sp'= list('Accession', 'Entry Name', 'Protein Name', 'Organism Name',
         'Organism Identifier', 'Gene Name', 'Protein Existence',
         'Sequence Version'),
  'tr'= list('Accession', 'Entry Name', 'Protein Name', 'Organism Name',
             'Organism Identifier', 'Gene Name', 'Protein Existence',
             'Sequence Version'),
  'pdb'= list('Entry', 'Chain'),
  'pat'= list('Country', 'Number'),
  'bbs'= list('Number'),
  'gnl'= list('Database', 'Identifier'),
  'ref'= list('Accession', 'Locus'),
  'lcl'= list('Identifier'),
  'nxp'= list('Identifier', 'Gene Name', 'Protein Name', 'Isoform Name')
)

ReadFasta <- function(datapath_or_datastring, is_datafile = TRUE) {
  # Read a file in FASTA format, either from a file or from a string of raw
  # data.
  # :param (string) datapath_or_datastring: Either the path to the FASTA file (can be relative
  #                                         or absolute), or a string corresponding to the content
  #                                         of a FASTA file (including newline characters).
  # :param (bool, optional) is_datafile: Either True (default) if passing the filepath to the data,
  #                                      or False if passing a string of raw data.
  # :rtype (list[dict]): A list of protein objects, each containing a
  #                      description (based on the header line) and the amino
  #                      acid sequence with, optionally, all non-amino-acid
  #                      letters removed.
  
  
  #ensure required argument is a string
  if (!typeof(datapath_or_datastring) == "character") {
    paste0("Please pass either the filepath to the data, or the data as a string")
  }
  
  raw_data <- list()
  
  #open file if given a path
  if (is_datafile) {
    file <- read_file(datapath_or_datastring)
    lines <- strsplit(file, '\n')
    if (!grepl(">", substring(lines[[1]][1], 1, 1))) {
      raw_data <- append(lines[[1]], ">", 0)
    } else {
      raw_data <- append(raw_data, lines)
    }
    
  } else{
    lines <- strsplit(datapath_or_datastring, '\n')
    if (!grepl(">", substring(lines[[1]][1], 1, 1))) {
      raw_data <- append(lines[[1]], ">", 0)
    } else {
      raw_data <- append(raw_data, lines)
    }
  }
  
  raw_data <- gsub("\r", "", raw_data[[1]])
  temp_file <- tempfile(pattern = "", fileext = ".fasta")
  lapply(lines, write, temp_file, append = TRUE)
  records <- readAAStringSet(temp_file)
  fasta_data = list()
  
  for (i in 1:length(records)) {
    blocks <-
      list("description" = DecodeDescription(names(records[i])),
           "sequence" = gsub("\r", "", records[[i]]))
    fasta_data[[length(fasta_data) + 1]] <- blocks
  }
  
  return(fasta_data)
}

DecodeDescription <- function(description) {
  # Parse the first line of a FASTA file using the specifications of
  # several different database headers (in _DATABASES).
  # :param (string) description: The header line with the initial '>'
  # removed.
  # :rtype (dict): A dictionary for which each key-value pair comprises
  # a property specified by the database used and the
  # value of that property given by the header. If the
  # database is not recognized, the keys are given as
  # 'desc-n' where n is the position of the property.
  
  if (is.null(description)) {
    return(list('-1' = "no description"))
  }
  
  decoded = list()
  desc <- strsplit(description, "|", fixed = TRUE)
  if (desc[[1]][1] %in% names(DATABASES)) {
    db_info <- DATABASES[[desc[[1]][1]]]
    if (grepl('sp|tr', desc[[1]][1])) {
      decoded[['Accession']] = desc[[1]][2]
      rs <-
        str_match_all(desc[[1]][3],
                      '([^\\s]+)(.*)\\ OS=(.*)\\ OX=(.*)\\ GN=(.*)\\ PE=(.*)\\ SV=(.*)$')
      for (i in 3:length(rs[[1]])) {
        el <- setNames(rs[[1]][i], db_info[[i]])
        decoded <- append(decoded, el)
      }
    } else{
      # shift by one, since first section in header describes the database
      for (i in 1:(length(desc[[1]]) - 1)) {
        el <- setNames(desc[[1]][i + 1], db_info[[i]])
        decoded <- append(decoded, el)
      }
    }
  } else {
    if (length(desc[[1]]) > 1) {
      for (i in 1:(length(desc[[1]]) - 1)) {
        el <- setNames(desc[[1]][i + 1], as.character(i))
        decoded <- append(decoded, el)
      }
    } else{
      decoded['Header'] = desc[[1]]
    }
  }
  return(decoded)
}