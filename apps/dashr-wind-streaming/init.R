# R script to run author supplied code, typically used to install additional R packages
# contains placeholders which are inserted by the compile script
# NOTE: this script is executed in the chroot context; check paths!

r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org"
options(repos=r)

# ======================================================================

# packages go here
install.packages("remotes")

remotes::install_github("plotly/dash-html-components", dependencies=TRUE)
remotes::install_github("plotly/dash-core-components", dependencies=TRUE)
remotes::install_github("plotly/dashR", dependencies=TRUE)

install.packages(c("DBI", "RSQLite", "stringr", "glue", "data.table", "plotly", "VGAM"))
