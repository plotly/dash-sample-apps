# R script to run author supplied code, typically used to install additional R packages
# contains placeholders which are inserted by the compile script
# NOTE: this script is executed in the chroot context; check paths!

r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org"
options(repos=r)



# ======================================================================
# packages go here
remotes::install_github("plotly/dashR", dependencies=TRUE, upgrade=TRUE)
remotes::install_github("plotly/rasterly", ref = "64f215b", upgrade=TRUE)

install.packages("stringr")
install.packages("rentrez")
install.packages("ape")
install.packages("readr")
install.packages("Biostrings")

