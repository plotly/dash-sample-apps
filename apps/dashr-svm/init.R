# R script to run author supplied code, typically used to install additional R packages
# contains placeholders which are inserted by the compile script
# NOTE: this script is executed in the chroot context; check paths!

r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org"
options(repos=r)

remotes::install_github("plotly/dash-table", ref="4.0.0-cran")
remotes::install_github("plotly/dash-html-components", ref="1.0.0-cran")
remotes::install_github("plotly/dash-core-components", ref="1.0.0-cran")
remotes::install_github("plotly/dashR", ref="0.1.0-cran")

install.packages("ipred", dependencies = TRUE)
install.packages("caret", dependencies = TRUE)
install.packages("kernlab", dependenceis = TRUE)
install.packages("ROCR", dependencies = TRUE)
install.packages("e1071", dependencies = TRUE)

