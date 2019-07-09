# R script to run author supplied code, typically used to install additional R packages
# contains placeholders which are inserted by the compile script
# NOTE: this script is executed in the chroot context; check paths!

r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org"
options(repos=r)

# ======================================================================

# installs askpass, assertthat, base64enc, cli, colorspace, crosstalk, curl, data.table, dplyr, fansi, ggplot2, glue, gtable, hexbin, htmlwidgets, httr, labeling, lattice, lazyeval, mgcv, munsell, nlme, openssl, pillar, pkgconfig, plogr, plyr, purrr, RColorBrewer, reshape2, scales, stringi, stringr, sys, tibble, tidyr, tidyselect, utf8, viridisLite, withr, yaml
remotes::install_version("plotly", version = "4.9.0", repos = "http://cloud.r-project.org", upgrade="never")

remotes::install_github("plotly/dashR", dependencies=FALSE)
remotes::install_github("plotly/dash-html-components")
remotes::install_github("plotly/dash-core-components")