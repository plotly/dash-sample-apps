r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org"
options(repos=r)

remotes::install_github("plotly/dash-html-components")
remotes::install_github("plotly/dash-core-components")
remotes::install_github("plotly/dash-table")
remotes::install_github("plotly/dashR")
remotes::install_github("plotly/dash-bio")
remotes::install_github("plotly/dash-daq")
