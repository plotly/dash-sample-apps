# The original python dash-svm-explorer pulls datasets from sklearn
# The package clusterSim in CRAN can be used to generate moons and circles
# the data is generated here and saved locally.

# required dependency installation:
# BiocManager::install("genefilter")
# install.packages("clusterSim")

library(clusterSim)
library(data.table)


# Functions to Generate Datasets:
generate_moons <- function(n){
  moons <- clusterSim::shapes.two.moon(n / 2)
  moons_dt <- data.table(cbind(moons$data, moons$clusters))
  colnames(moons_dt) <- c("X1", "X2", "Y")
  moons_dt[, Y := Y - 1]
  moons_dt
}
generate_circles <- function(n){
  circles <- clusterSim::shapes.bulls.eye(n / 2)
  circles_dt <- data.table(cbind(circles$data, circles$clusters))
  colnames(circles_dt) <- c("X1", "X2", "Y")
  circles_dt[, Y := Y - 1]
  circles_dt
}
generate_lin_sep <- function(n, noise){
  linsep <- list()
  x <- sample(c(rep(0, n / 2), rep(1, n / 2)))
  linsep <- data.table(cbind(cbind(x, x), x))
  colnames(linsep) <- c("X1", "X2", "Y")
  linsep
}

DTm <- generate_moons(500)
DTc <- generate_circles(500)
DTl <- generate_lin_sep(500)

allDT <- list(moons = DTm, circles = DTc, linear = DTl)

saveRDS(object = allDT, file = "data.rds", version = 2)
