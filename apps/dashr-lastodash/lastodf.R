library(readr)
library(stringr)
library(gdata)
library(tidyr)
library(dplyr)

LASobject <- setRefClass("LASobject", fields = list(
  V="data.frame", 
  W="data.frame", 
  C="data.frame", 
  P="data.frame", 
  O="data.frame", 
  A="data.frame")
)

convertLAS <- function(LASpath) {
  LASdata <- LASobject()
  LASstring <- read_file(LASpath)
  LASsections <- str_split(LASstring, "~")
  for (section in LASsections[[1]]) {
    # section contains version and wrap mode information
    if (startsWith(section, "V")) {
      print("~V")
      linesV <- str_split(section, "\r\n")
      contentV <- linesV[[1]][-c(1, length(linesV[[1]]))]
      LASdata$V <- data.frame(x = contentV) %>% 
        separate(col=x, into=c("label", "x"), sep="\\.", extra="merge") %>% 
        separate(col=x, into=c("value", "description"), sep=":", extra="merge")
    }
    # section contains well identification
    if (startsWith(section, "W")) {
      print("~W")
      linesW <- str_split(section, "\r\n")
      contentW <- linesW[[1]][-c(1, length(linesW[[1]]))]
      if ((length(contentW) >= 2) & startsWith(contentW, "#")) {
        contentW <- contentW[-c(1, 2)]
      }
      LASdata$W <- data.frame(x = contentW) %>% 
        separate(col=x, into=c("mnemonic", "x"), sep="\\.", extra="merge") %>% 
        separate(col=x, into=c("unit", "x"), sep=" ", extra="merge") %>% 
        separate(col=x, into=c("value", "description"), sep=":", extra="merge")
    }
    # section contains curve information
    if (startsWith(section, "C")) {
      print("~C")
      linesC <- str_split(section, "\r\n")
      contentC <- linesC[[1]][-c(1, length(linesC[[1]]))]
      if ((length(contentC) >= 2) & startsWith(contentC, "#")) {
        contentC <- contentC[-c(1, 2)]
      }
      LASdata$C <- data.frame(x = contentC) %>%
      separate(col=x, into=c("mnemonic", "x"), sep="\\.", extra="merge") %>%
      separate(col=x, into=c("unit", "x"), sep=" ", extra="merge") %>%
      separate(col=x, into=c("api code", "curve description"), sep=":", extra="merge")
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
      LASdata$A <- read.table(text=contentA, header=TRUE, sep="")
    }
  }
  return (LASdata)
}