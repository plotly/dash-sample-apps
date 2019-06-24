# require(tidyverse)
# require(dashR)
# require(dashCoreComponents)
# require(dashHtmlComponents)
# require(dashPlayer)
# require(plotly)
# require(rlist)
# require(stringr)

# require(data.table)
# require(glue)
# require(compiler)

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

# installs magrittr, promises, R6
remotes::install_version("httpuv", repos = "http://cloud.r-project.org")

# installs crayon, digest, htmltools, mime, sourcetools, xtable
remotes::install_version("shiny", repos = "http://cloud.r-project.org", upgrade="never")

# installs askpass, assertthat, base64enc, cli, colorspace, crosstalk, curl, data.table, dplyr, fansi, ggplot2, glue, gtable, hexbin, htmlwidgets, httr, labeling, lattice, lazyeval, mgcv, munsell, nlme, openssl, pillar, pkgconfig, plogr, plyr, purrr, RColorBrewer, reshape2, scales, stringi, stringr, sys, tibble, tidyr, tidyselect, utf8, viridisLite, withr, yaml
remotes::install_version("plotly", repos = "http://cloud.r-project.org", upgrade="never")

install.packages("https://cloud.r-project.org/src/contrib/assertthat_0.2.1.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/xml2_1.2.0.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/triebeard_0.3.0.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/Archive/urltools/urltools_1.7.2.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/jsonlite_1.6.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/webutils_0.6.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/brotli_1.2.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/reqres_0.2.2.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/uuid_0.1-2.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/base64enc_0.1-3.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/codetools_0.2-16.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/globals_0.12.4.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/Archive/future/future_1.11.1.1.tar.gz", type="source", repos=NULL)

# fiery and friends
install.packages("https://cloud.r-project.org/src/contrib/routr_0.3.0.tar.gz", type="source", repos=NULL)
install.packages("https://cloud.r-project.org/src/contrib/fiery_1.1.1.tar.gz", type="source", repos=NULL)

remotes::install_github("plotly/dash-table", ref="4.0.0-cran")
remotes::install_github("plotly/dash-html-components", ref="1.0.0-cran")
remotes::install_github("plotly/dash-core-components", ref="1.0.0-cran")
remotes::install_github("plotly/dashR", ref="0.1.0-cran")
install.packages("dashPlayer_0.0.1.tar.gz", type="source", repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/data.table_1.12.2.tar.gz', repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/glue_1.3.1.tar.gz', repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/stringr_1.4.0.tar.gz', repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/magrittr_1.5.tar.gz', repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/stringi_1.4.3.tar.gz', repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/jsonlite_1.6.tar.gz', repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/iterators_1.0.10.tar.gz', repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/foreach_1.4.4.tar.gz', repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/XML_3.98-1.20.tar.gz', repos=NULL)
install.packages('https://cran.r-project.org/src/contrib/rlist_0.4.6.1.tar.gz', repos=NULL)

##DashPlayer- install the local version inside the repo
