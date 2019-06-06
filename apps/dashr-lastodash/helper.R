library(readr)
library(stringr)
library(gdata)
library(tidyr)
library(dplyr)

las_object <- setRefClass("las_object", fields = list(
  V="data.frame", 
  W="data.frame", 
  C="data.frame", 
  P="data.frame", 
  O="data.frame", 
  A="data.frame")
)

convert_las <- function(las_path) {
  las_data <- las_object()
  las_string <- read_file(las_path)
  las_sections <- str_split(las_string, "~")
  for (section in las_sections[[1]]) {
    # section contains version and wrap mode information
    if (startsWith(section, "V")) {
      print("~V")
      linesV <- str_split(section, "\r\n")
      contentV <- linesV[[1]][-c(1, length(linesV[[1]]))]
      las_data$V <- data.frame(x = contentV) %>% 
        separate(col=x, into=c("label", "x"), sep="\\.", extra="merge") %>% 
        separate(col=x, into=c("value", "description"), sep=":\\s+", extra="merge")
    }
    # section contains well identification
    if (startsWith(section, "W")) {
      print("~W")
      linesW <- str_split(section, "\r\n")
      contentW <- linesW[[1]][-c(1, length(linesW[[1]]))]
      if ((length(contentW) >= 2) & startsWith(contentW, "#")) {
        contentW <- contentW[-c(1, 2)]
      }
      #las_data$W <- read.table(text=contentW, header=FALSE, sep="", strip.white=TRUE) %>%  
      las_data$W <- data.frame(x = contentW) %>%
        separate(col=x, into=c("mnemonic", "x"), sep="\\.", extra="merge") %>% 
        separate(col=x, into=c("unit", "x"), sep="\\s+", extra="merge") %>% 
        separate(col=x, into=c("value", "description"), sep=":\\s+", extra="merge")
    }
    # section contains curve information
    if (startsWith(section, "C")) {
      print("~C")
      linesC <- str_split(section, "\r\n")
      contentC <- linesC[[1]][-c(1, length(linesC[[1]]))]
      if ((length(contentC) >= 2) & startsWith(contentC, "#")) {
        contentC <- contentC[-c(1, 2)]
      }
      las_data$C <- data.frame(x = contentC) %>%
        separate(col=x, into=c("mnemonic", "x"), sep="\\.", extra="merge") %>%
        separate(col=x, into=c("unit", "x"), sep="\\s+", extra="merge") %>%
        separate(col=x, into=c("api code", "curve description"), sep=":\\s+", extra="merge")
    }
    # section contains parameters or constants
    if (startsWith(section, "P")) {
      print("~P")
    }
    # section contains other information such as comments
    if (startsWith(section, "O")) {
      print("~O")
    }
    # section contains ASCII log data
    if (startsWith(section, "A")) {
      print("~A")
      #linesA <- str_split(section, "\r\n")
      #contentA <- linesA[[1]][-c(1, length(linesC[[1]]))]
      #tableA <- read.table(text=contentA, header=FALSE, sep="", col.names=as.vector(tableC['mnemonic']))
      contentA <- str_sub(section, 2, -2)
      las_data$A <- read.table(text=contentA, header=TRUE, sep="")
      #las_data$A <- as.data.frame(apply(dataA,2,function(x)str_trim(x,"left")))
      #las_data$A <- sub("[[:space:]]+$","",read.table(text=contentA, header=TRUE, sep=""))
    }
  }
  return (las_data)
}

get_item <- function(df, name1, name2, key1) {
  key2 <- ""
  for (i in 1:nrow(df)) {
    if (str_trim(df[i,name1],"left")==key1) {
      key2 <- df[i,name2]
    }
  }
  return (key2)
}

get_plot_list <- function(nplots) {
  lapply(seq_len(nplots), function(x) plot_ly())
}
