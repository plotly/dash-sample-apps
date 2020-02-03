r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org"
options(repos=r)

install.packages("remotes")

install.packages("jsonlite")
install.packages("listenv")
install.packages("tm")
install.packages("topicmodels")
install.packages("tidytext")
install.packages("Rtsne")
install.packages("rlist")
install.packages("RColorBrewer")

# installs magrittr, promises, R6
remotes::install_version("httpuv", version = "1.4.5.1", repos = "http://cloud.r-project.org", upgrade="never")

# installs crayon, digest, htmltools, mime, sourcetools, xtable
remotes::install_version("shiny", version = "1.2.0", repos = "http://cloud.r-project.org", upgrade="never")

# installs askpass, assertthat, base64enc, cli, colorspace, crosstalk, curl, data.table, dplyr, fansi, ggplot2, glue, gtable, hexbin, htmlwidgets, httr, labeling, lattice, lazyeval, mgcv, munsell, nlme, openssl, pillar, pkgconfig, plogr, plyr, purrr, RColorBrewer, reshape2, scales, stringi, stringr, sys, tibble, tidyr, tidyselect, utf8, viridisLite, withr, yaml
remotes::install_version("plotly", version = "4.9.0", repos = "http://cloud.r-project.org", upgrade="never")

# fiery and friends
install.packages("routr", type="source")
install.packages("fiery", type="source")


#remotes::install_github("plotly/dashR", dependencies=FALSE)
#remotes::install_github("plotly/dash-html-components")
#remotes::install_github("plotly/dash-core-components")
#remotes::install_github("plotly/dash-table")

remotes::install_github("plotly/dashR", upgrade = TRUE)

