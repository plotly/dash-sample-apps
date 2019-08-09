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



read_fasta <- function(datapath_or_datastring, is_datafile=TRUE){
  
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
  if(!typeof(datapath_or_datastring)=="character"){
    paste0("Please pass either the filepath to the data, or the data as a string")
  }
  
  raw_data <- list()
  
  #open file if given a path
  if(is_datafile){
    file <- read_file(datapath_or_datastring)
    lines <- strsplit(file,'\n')
    if(!grepl(">",substring(lines[[1]][1],1,1))){
      raw_data <- append(lines[[1]],">",0)
    } else {
      raw_data <- append(raw_data,lines)
    }
    
  } else{
    lines <- strsplit(datapath_or_datastring,'\n')
    if(!grepl(">",substring(lines[[1]][1],1,1))){
      raw_data <- append(lines[[1]],">",0)
    } else {
      raw_data <- append(raw_data,lines)
    }
  }

  raw_data <- gsub("\r","",raw_data[[1]])
  temp_file <- tempfile(pattern = "",fileext=".fasta")
  lapply(lines, write, temp_file, append=TRUE)
  records <- readAAStringSet(temp_file)
  fasta_data = list()
  
  for (i in 1:length(records)){
    blocks <-list("description"=decode_description(names(records[i])),"sequence"=gsub("\r","",records[[i]]))
    fasta_data[[length(fasta_data)+1]] <- blocks
  }
  
  return(fasta_data)
}

decode_description <- function(description){
  # Parse the first line of a FASTA file using the specifications of
  # several different database headers (in _DATABASES).
  # :param (string) description: The header line with the initial '>'
  # removed.
  # :rtype (dict): A dictionary for which each key-value pair comprises
  # a property specified by the database used and the
  # value of that property given by the header. If the
  # database is not recognized, the keys are given as
  # 'desc-n' where n is the position of the property.
  
  if(is.null(description)){
    return(list('-1'="no description"))
  }
    
  decoded = list()
  #description <- names(op[1])
  desc <- strsplit(description,"|",fixed=TRUE)
  if(desc[[1]][1] %in% names(DATABASES)){
    db_info <- DATABASES[[desc[[1]][1]]]
    if(grepl('sp|tr',desc[[1]][1])){
      decoded[['Accession']] = desc[[1]][2]
      rs <- str_match_all(desc[[1]][3], '([^\\s]+)(.*)\\ OS=(.*)\\ OX=(.*)\\ GN=(.*)\\ PE=(.*)\\ SV=(.*)$')
      for(i in 3:length(rs[[1]])){
        el<-setNames(rs[[1]][i],db_info[[i]])
        decoded <- append(decoded,el)
      }
    } else{
      # shift by one, since first section in header describes the database
      for(i in 1:(length(desc[[1]])-1)){
        el<-setNames(desc[[1]][i+1],db_info[[i]])
        decoded <- append(decoded,el)
      }
    }
  } else {
    if(length(desc[[1]])>1){
      for(i in 1:(length(desc[[1]])-1)){
        el<-setNames(desc[[1]][i+1],as.character(i))
        decoded <- append(decoded,el)
      }
    } else{
      decoded['Header']=desc[[1]]
    }
  }
  return(decoded)
}



# dataset3
# kk<-readAAStringSet("test3.fasta")
# 
# k<-readAAStringSet("./assets/sample_data/sequence_viewer_NX_P02768.fasta")
# read_file("./assets/sample_data/sequence_viewer_NX_P02768.fasta")
# names(k)
# paste0(k[[3]])
# lines <- read_file("./assets/sample_data/sequence_viewer_NX_P02768.fasta")
# lines
# lines <- strsplit(lines,'\n')
# lines
# lines <-">nxp|NX_P02768-1|ALB|Serum albumin|Iso 1\r\nMKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\r\nEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\r\nERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\r\nFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\r\nARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\r\nECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\r\nRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\r\nQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\r\nLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\r\nSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\r\nAASQAALGL\r\n>nxp|NX_P02768-2|ALB|Serum albumin|Iso 2\r\nMKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKAWAVARLSQRFPKAEF\r\nAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSH\r\nCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLR\r\nLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALL\r\nVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTP\r\nVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTA\r\nLVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL\r\n>nxp|NX_P02768-3|ALB|Serum albumin|Iso 3\r\nMKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\r\nEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\r\nERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYETTLEKCCAAADPHEC\r\nYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRN\r\nLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFS\r\nALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMD\r\nDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL\r\n"
# lines <- strsplit(lines,'\n')
# lines[[1]]
# lines<-gsub("\r","",lines[[1]])
# lines
# lapply(lines, write, "test.fasta", append=TRUE)
# sink("test.fasta")
# #writeLines(unlist(lapply(lines, paste, collapse="")))
# sink()
# read_file("test.fasta")
# x<-paste(x, collapse = '\n')
# sink()
# read_file("test.fasta")
# 
# gsub("\r","",read_file("test.fasta"))
# raw_data <- list()
# raw_data <- append(raw_data,lines)
# raw_data
# lapply(lines, writeLines, "test2.txt")
# read_file("test2.txt")
# ll <- readAAStringSet("test2.txt")
# paste0(ll[[2]])
# lapply(lines, cat, file="test2.txt", append=TRUE)
# read_file("test2.txt")
# 
# raw_data[[2]]
# 1:(length(desc[[1]])-1)
# desc[[1]]
# length(desc[[1]])-1
# kk <- str_match_all(desc[[1]][3], '([^\\s]+)(.*)\\ OS=(.*)\\ OX=(.*)\\ GN=(.*)\\ PE=(.*)\\ SV=(.*)$')
# kk[[1]][4]
# grep('([^\\s]+)(.*)\\ OS=(.*)\\ OX=(.*)\\ GN=(.*)\\ PE=(.*)\\ SV=(.*)$',desc[[1]][3])
# 
# 
# for(i in 3:length(db_info)){
#    db_info[[i]] = kk[[1]][i]
#    decoded <- append(decoded,db_info[[i]])
# }
# kk[[1]][3]
# db_info
# for(i in 1:(length(desc[[1]])-1)){
#   decoded <- append(decoded, as.character(desc[[1]][i+1]))
# }
# decoded
# desc[[1]][2]
# 
# for(i in 1:(length(desc[[1]])-1)){
#   print(i)
# }
# length(op)
# paste0(records[1])
# paste0(gsub("\r","",records[[3]]))
# ok <- list()
# for (i in 1:length(op)){
#   #ll<-setNames("description",decode_description(names(op[i])))
#   ll <-list("description"=decode_description(names(op[i])),"sequence"=paste0(op[i]))
#   ok[[length(ok)+1]] <- ll
# }
# 
# decode_description(names())
# 
# ok[[1]]
# 
# paste0(op[[1]])
# k<-list(description=decode_description(names(op[1])))
# k
# names(op[i])
