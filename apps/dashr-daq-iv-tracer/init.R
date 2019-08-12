# R script to run author supplied code, typically used to install additional R packages
# contains placeholders which are inserted by the compile script
# NOTE: this script is executed in the chroot context; check paths!

r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org"
options(repos=r)

# ======================================================================
# packages go here
install.packages("remotes")

# installs Rcpp, rlang, BH
install.packages("later")
install.packages("jsonlite")

install.packages("listenv")

install.packages("magrittr")

install.packages("manhattanly")

install.packages("RCurl")
install.packages("readr")


# fiery and friends
install.packages("https://cloud.r-project.org/src/contrib/routr_0.3.0.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/fiery_1.1.1.tar.gz", type="source", repos=NULL)




remotes::install_github("plotly/dashR", dependencies=TRUE)
remotes::install_github("plotly/dash-html-components")
remotes::install_github("plotly/dash-core-components")
remotes::install_github("plotly/dash-daq", dependencies=TRUE)

