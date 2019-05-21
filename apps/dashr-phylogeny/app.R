library(dashR)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(data.table)
library(ape)

# TODO: add flu virus?
species <- list(
  "Avian", "Ebola", "Dengue", "Lassa", "Measles", "Mumps", "Zika"
)

# Get necessary data for world map (Exclude antarctica)
world <- fread(
  "data/country_iso_list_long_lat_iso3letters.csv",
  key = "Country"
)
world <- world[!Country == "Antarctica"]

# Get X coordinates for nodes and clades
getXCoordinates <- function(tree){
  if (is.rooted(tree)){
    xcoords <- node.depth.edgelength(tree)
  } else {
    xcoords <- node.depth(tree)
  }
  xcoords
}

# Get Y coordinates for nodes and clades
getYCoordinates <- function(tree, dist = 1.3){
  ycoords <- node.height(tree)
  ycoords
}

# Draw the tree branches
getCladeLines <- function(tree, fig){
  # Adapted from plot.phylo function in package "ape"
  xx <- node.depth.edgelength(tree)
  yy <- node.height(tree)
  Ntip <- Ntip(tree)
  Nnode <- Nnode(tree)
  nodes <- (Ntip + 1):(Ntip + Nnode)
  x0v <- xx[nodes]
  y0v <- y1v <- numeric(Nnode)
  edge <- tree$edge
  NodeInEdge1 <- vector("list", Nnode)
  e1 <- edge[, 1]
  for (i in seq_along(e1)){
    j <- e1[i] - Ntip
    NodeInEdge1[[j]] <- c(NodeInEdge1[[j]], i)
  }
  for (i in 1:Nnode){
    j <- NodeInEdge1[[i]]
    tmp <- range(yy[edge[j, 2]])
    y0v[i] <- tmp[1]
    y1v[i] <- tmp[2]
  }
  x0h <- xx[edge[, 1]]
  x1h <- xx[edge[, 2]]
  y0h <- yy[edge[, 2]]
  fig %>%
    add_segments(
      inherit = FALSE, x0h, y0h, x1h, y0h,
      line = list(color = "black", width = 1),
      showlegend = FALSE
      ) %>%
    add_segments(
      inherit = FALSE, x0v, y0v, x0v, y1v,
      line = list(color = "black", width = 1),
      showlegend = FALSE
    )
}

# Create the tree graph
createTree <- function(virus_name, tree_file, metadata_file){
  tree <- read.tree(file = tree_file)
  x_coords <- getXCoordinates(tree)
  y_coords <- getYCoordinates(tree)
  DT <- data.table(
    x = x_coords,
    y = y_coords,
    Strain = c(tree$tip.label, tree$node.label)
  )
  DT[, Clade := !Strain %like% "NODE_.*"]
  df <- fread(metadata_file)
  regions <- fread("utils/regions.csv")
  DT2 <- merge(DT, df, by = "Strain", all.x = TRUE)
  DT2 <- merge(DT2, regions, by = "Country", all.x = TRUE)
  nb_genome <- DT2[Clade == TRUE, .N]
  graph_title <- sprintf(
    paste0(
      "Phylogeny of %s Virus<br>",
      " %s genomes colored according to region and country"
    ),
    paste0(
      toupper(substring(virus_name, 1, 1)),
      substring(virus_name, 2)
    ),
    nb_genome
  )
  p <- plot_ly(
    data = DT2[Clade == "FALSE"],
    x = ~x, y = ~y,
    color = I("#646464"),
    type = "scatter", mode = "markers",
    name = "</br> Nodes",
    hoverinfo = "text",
    text = ~paste(
      sprintf("Node: %s </br>", Strain),
      sprintf("X: %s", round(x, 3)),
      sprintf("Y: %s", round(y, 3)),
      sep = "</br> "
    )
  ) %>%
    add_markers(
      inherit = FALSE,
      data = DT2[Clade == "TRUE"], x = ~x, y = ~y,
      color = ~Country,
      type = "scatter", mode = "markers",
      hoverinfo = "text",
      text = ~paste(
        sprintf("%s </br>", Strain),
        sprintf("Country: %s", Country),
        sprintf("Region: %s", region),
        sprintf("Collection Date: %s", Date),
        sprintf("Journal: %s", Journal),
        sprintf("Authors: ", Authors),
        sep = "</br>"
      )
    ) %>%
    layout(
      title = graph_title,
      font = list(family = "Roboto", size = "14"),
      hovermode = "closest",
      margin = list(
        t = 70
      ),
      autosize = TRUE,
      xaxis = list(
        showline = FALSE,
        zeroline = FALSE,
        showgrid = FALSE,
        showticklabels = FALSE,
        title = ""
      ),
      yaxis = list(
        showline = FALSE,
        zeroline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        title = ""
      )
    )
  getCladeLines(tree, p)
}

# Create the map graph
createMapBubbleYear <- function(virus_name, metadata_file_stat,
                                 min_date, max_date){
  df <- fread(metadata_file_stat)
  # replace USA with United States
  df[df$Country == "USA", "Country"] <- "United States"
  df <- df[(Year >= min_date) & (Year <= max_date)]
  df <- df[, .(count = .N, code = ISO3[1]), by = Country]
  full_df <- merge(
    df, world,
    by.x = "Country", by.y = "Country", all.y = TRUE
  )
  full_df[is.na(full_df)] <- 0
  g <- list(
    showframe = FALSE,
    showcoastlines = FALSE,
    projection = list(type = "Mercator")
  )
  plot_geo(full_df) %>%
    add_trace(
      z = ~count, color = ~count,
      colorscale = "Blues", locations = ~ISO3,
      reversescale = TRUE,
      hoverinfo = "text",
      text = ~sprintf(
        "%s<br>Number of cases: %s",
        Country, count
      )
    ) %>%
    colorbar(
      title = sprintf(
        "Number of %s Virus Cases",
        paste0(
          toupper(substring(virus_name, 1, 1)),
          substring(virus_name, 2)
        )
      ),
      titleside = "right",
      autotick = FALSE,
      len = 1
    ) %>%
    layout(
      title = sprintf(
        "%s Virus Cases<br>Between %s and %s",
        paste0(
          toupper(substring(virus_name, 1, 1)),
          substring(virus_name, 2)
        ),
        min_date, max_date
      ),
      font = list(family = "Roboto", size = "14"),
      margin = list(
        t = 35
      ),
      geo = g,
      autosize = TRUE
    )
}

# Create the curve graph
createCurveLine <- function(metadata_file_stat, virus_name,
                            min_date, max_date){
  df <- fread(metadata_file_stat)
  df <- df[(Year >= min_date) & (Year <= max_date)]
  casesByYear <- CJ(min_date:max_date, unique(df[, Country]))
  colnames(casesByYear) <- c("Year", "Country")
  casesByYear <- merge(
    casesByYear, df[, .N, c("Year", "Country")],
    all.x = TRUE
  )
  casesByYear[is.na(casesByYear)] <- 0
  plot_ly(
    data = casesByYear, x = ~Year, y = ~N,
    type = "scatter", mode = "markers+lines",
    color = ~Country,
    marker = list(color = "black"),
    line = list(width = 1)
  ) %>%
    layout(
      title = sprintf(
        paste0(
          "Number of %s Virus Cases<br>",
          "Between %s and %s by Country"
        ),
        paste0(
          toupper(substring(virus_name, 1, 1)),
          substring(virus_name, 2)
        ),
        min_date, max_date
      ),
      xaxis = list(
        autotick = FALSE,
        tickmode = "array",
        tickvals = seq(min_date, max_date, slicer(min_date, max_date)$step),
        dtick = 1
      ),
      yaxis = list(
        title = sprintf(
          "Number of %s Virus Cases",
          paste0(
            toupper(substring(virus_name, 1, 1)),
            substring(virus_name, 2)
          )
        )
      )
    )
}

# Create histogram
createHistogram <- function(metadata_file_stat, virus_name,
                            min_date, max_date){
  df <- fread(metadata_file_stat)
  df <- df[(Year >= min_date) & (Year <= max_date)]
  plot_ly(
    data = df, x = ~Country,
    histfunc = "count", type = "histogram",
    hoverinfo = "y"
  ) %>%
    layout(
      title = sprintf(
        paste0(
          "Distribution of %s Virus Cases<br>",
          "Between %s and %s"
        ),
        paste0(
          toupper(substring(virus_name, 1, 1)),
          substring(virus_name, 2)
        ),
        min_date, max_date
      ),
      xaxis = list(
        title = "",
        tickangle = 45
      )
    )
}

createPathsFile <- function(virus_name, ...){
  dir <- paste0("data/", virus_name, "/")
  lvls_dir <- paste(c(list(...), ""), collapse = "/")
  lvls_file <- paste(c("", list(...)), collapse = "_")
  tree_file <- paste0(
    dir, lvls_dir, "nextstrain_", virus_name, lvls_file, "_tree.new"
  )
  metadata_file <- paste0(
    dir, lvls_dir, "nextstrain_", virus_name, lvls_file, "_metadata.csv"
  )
  stat_file <- paste0(
    dir, lvls_dir, "stat_year_nextstrain_",
    virus_name, lvls_file, "_metadata.csv"
  )
  list(
    tree_file = tree_file,
    metadata_file = metadata_file,
    stat_file = stat_file
  )
}

slicer <- function(min_date, max_date){
  n_years <- max_date - min_date
  step <- ifelse(
    n_years < 5,
    1, ifelse(
      n_years <= 10,
      2, ifelse(
        n_years <= 50,
        5, 10
      )
    )
  )
  marks <- as.list(
    setNames(
      seq(min_date, max_date, by = step),
      seq(min_date, max_date, by = step)
    )
  )
  list(marks = marks, step = step)
}

virus_name <- "measles"
files <- createPathsFile(virus_name)
tree_file <- files$tree_file
metadata_file_stat <- files$stat_file
metadata_file <- files$metadata_file
min_date <- 1977
max_date <- 2016

#######################################################################
app <- Dash$new(name = "DashR Phylogeny")

app$layout(
  htmlDiv(
    list(
      htmlDiv(
        list(
          htmlDiv(
            list(
              htmlH2(
                "Phylogeny trees and global spread of 6 viruses",
                id = "title"
              ),
              htmlImg(
                src = "https://user-images.githubusercontent.com/37411533/56820587-55d56e80-681a-11e9-91b6-026e1551b338.png"
              )
            )
          )
        ),
        className = "banner"
      ),
      # Body
      htmlDiv(
        className = "container",
        list(
          htmlDiv(
            list(
              htmlDiv(
                list(
                  htmlDiv(
                    list(
                      htmlH6("Dataset"),
                      dccDropdown(
                        id = "d_virus-name",
                        options = lapply(
                          species, function(x){
                            list(label = x, value = x)
                          }
                        ),
                        value = "Measles"
                      ),
                      htmlDiv(
                        id = "output-container",
                        children = "Selected Virus"
                      ),
                      htmlDiv(
                        id = "controls-container_mumps",
                        children = list(
                          htmlSpan(
                            "Region:", style = list(fontStyle = "italic")
                          ),
                          dccDropdown(
                            id = "d_mumps",
                            options = list(
                              list(label = "Global", value = "global"),
                              list(label = "North America", value = "na")
                            ),
                            value = "global"
                          )
                        ),
                        style = list(display = "none")
                      ),
                      htmlDiv(
                        id = "controls-container_dengue",
                        children = list(
                          htmlSpan(
                            "Serotype:", style = list(fontStyle = "italic")
                          ),
                          dccDropdown(
                            id = "d_dengue",
                            options = lapply(
                              list("All", "Denv1", "Denv2", "Denv3"),
                              function(x){
                                list(label = x, value = tolower(x))
                              }
                            ),
                            value = "all"
                          )
                        ),
                        style = list(display = "none")
                      ),
                      htmlDiv(
                        id = "controls-container_lassa",
                        children = list(
                          htmlSpan(
                            "RNA:", style = list(fontStyle = "italic")
                          ),
                          dccDropdown(
                            id = "d_lassa",
                            options = lapply(
                              list("s", "l"),
                              function(x){
                                list(label = x, value = x)
                              }
                            ),
                            value = "s"
                          )
                        ),
                        style = list(display = "none")
                      ),
                      htmlDiv(
                        id = "controls-container_avian",
                        children = list(
                          htmlSpan(
                            "Subtype:", style = list(fontStyle = "italic")
                          ),
                          dccDropdown(
                            id = "d_avian_opt1",
                            options = list(
                              list(
                                label = "h7n9",
                                value = "h7n9",
                                display = "block"
                              )
                            ),
                            value = "h7n9"
                          ),
                          htmlSpan(
                            "RNA segment:", style = list(fontStyle = "italic")
                          ),
                          dccDropdown(
                            id = "d_avian_opt2",
                            # RNA segments?
                            options = lapply(
                              list(
                                "ha", "mp", "na", "ns",
                                "np", "pa", "pb2", "pb1"
                              ),
                              function(x){
                                list(label = x, value = x)
                              }
                            ),
                            value = "ha"
                          )
                        ),
                        style = list(display = "none")
                      ),
                      htmlDiv(
                        id = "controls-container_flu",
                        children = list(
                          dccDropdown(
                            id = "d_flu_opt1",
                            options = lapply(
                              list("h3n2", "h1n1pdm", "vic", "yam"),
                              function(x){
                                list(label = x, value = x)
                              }
                            ),
                            value = "h3n2"
                          ),
                          dccDropdown(
                            id = "d_flu_opt2",
                            options = lapply(
                              list("ha", "na"),
                              function(x){
                                list(label = x, value = x)
                              }
                            ),
                            value = "ha"
                          ),
                          dccDropdown(
                            id = "d_flu_opt3",
                            options = lapply(
                              list("2y", "3y", "6y", "12y"),
                              function(x){
                                list(label = x, value = x)
                              }
                            ),
                            value = "3y"
                          )
                        ),
                        style = list(display = "none")
                      ),
                      htmlBr(),
                      htmlBr(),
                      htmlH6("Data Range"),
                      htmlDiv(
                        id = "id-slicer",
                        children = list(
                          dccRangeSlider(
                            id = "id-year",
                            min = min_date,
                            max = max_date,
                            step = 1,
                            marks = slicer(min_date, max_date)$marks,
                            value = list(min_date, max_date)
                          )
                        ),
                        style = list(margin = "0 1.5rem")
                      ),
                      htmlBr(),
                      htmlBr(),
                      htmlDiv(id = "output-container-range-slider"),
                      htmlBr(),
                      htmlBr(),
                      dccGraph(
                        id = "curve-line-graph",
                        figure = createCurveLine(
                          metadata_file_stat,
                          virus_name,
                          min_date,
                          max_date
                        ),
                        style = list(maxHeight = 700)
                      )
                    )
                  )
                ),
                className = "four columns",
                style = list(marginTop = "1rem")
              ),
              htmlBr(),
              htmlDiv(
                className = "eight columns",
                style = list(marginTop = "10"),
                children = htmlDiv(
                  list(
                    htmlDiv(
                      id = "right-top-graph",
                      children = dccGraph(
                        id = "tree-graph",
                        figure = createTree(
                          virus_name,
                          tree_file,
                          metadata_file
                        ),
                        style = list(height = "1000")
                      )
                    )
                  ),
                  className = "row"
                )
              )
            )
          ),
          htmlDiv(
            list(
              htmlDiv(
                list(
                  dccGraph(
                    id = "graph_map",
                    figure = createMapBubbleYear(
                      virus_name,
                      metadata_file_stat,
                      min_date,
                      max_date
                    )
                  )
                ),
                className = "six columns"
              ),
              htmlDiv(
                list(
                  dccGraph(
                    id = "id-histo",
                    figure = createHistogram(
                      metadata_file_stat,
                      virus_name,
                      min_date,
                      max_date
                    )
                  )
                ),
                className = "six columns"
              )
            ),
            style = list(marginTop = "5rem"),
            className = "twelve columns"
          )
        )
      )
    )
  )
)

app$callback(
  output("output-container", "children"),
  list(input("d_virus-name", "value")),
  function(virus_name){
    sprintf("You have selected the %s virus", virus_name)
  }
)

app$callback(
  output("controls-container_mumps", "style"),
  list(input("d_virus-name", "value")),
  function(virus_name){
    if (virus_name == "Mumps"){
      return(list(display = "block"))
    }
    return(list(display = "none"))
  }
)

app$callback(
  output("controls-container_dengue", "style"),
  list(input("d_virus-name", "value")),
  function(virus_name){
    if (virus_name == "Dengue"){
      return(list(display = "block"))
    }
    return(list(display = "none"))
  }
)

app$callback(
  output("controls-container_lassa", "style"),
  list(input("d_virus-name", "value")),
  function(virus_name){
    if (virus_name == "Lassa"){
      return(list(display = "block"))
    }
    return(list(display = "none"))
  }
)

app$callback(
  output("controls-container_avian", "style"),
  list(input("d_virus-name", "value")),
  function(virus_name){
    if (virus_name == "Avian"){
      return(list(display = "block"))
    }
    return(list(display = "none"))
  }
)

app$callback(
  output("controls-container_flu", "style"),
  list(input("d_virus-name", "value")),
  function(virus_name){
    if (virus_name == "Flu"){
      return(list(display = "block"))
    }
    return(list(display = "none"))
  }
)

app$callback(
  output("tree-graph", "figure"),
  list(
    input("d_virus-name", "value"),
    input("d_mumps", "value"),
    input("d_dengue", "value"),
    input("d_lassa", "value"),
    input("d_avian_opt1", "value"),
    input("d_avian_opt2", "value"),
    input("d_flu_opt1", "value"),
    input("d_flu_opt2", "value"),
    input("d_flu_opt3", "value")
  ),
  function(virus_name, mumps, dengue, lassa,
           avian_opt1, avian_opt2, flu_opt1,
           flu_opt2, flu_opt3){
    virus_name <- tolower(virus_name)
    if (virus_name %in% c("ebola", "zika", "measles")){
      files <- createPathsFile(virus_name)
    } else if (virus_name == "mumps"){
      files <- createPathsFile(virus_name, mumps)
    } else if (virus_name == "dengue"){
      files <- createPathsFile(virus_name, dengue)
    } else if (virus_name == "lassa"){
      files <- createPathsFile(virus_name, lassa)
    } else if (virus_name == "avian"){
      files <- createPathsFile(
        virus_name, avian_opt1, avian_opt2
      )
    } else if (virus_name == "flu"){
      files <- createPathsFile(
        virus_name, flu_opt1, flu_opt2, flu_opt3
      )
    }
    createTree(virus_name, files$tree_file, files$metadata_file)
  }
)

app$callback(
  output("graph_map", "figure"),
  list(
    input("d_virus-name", "value"),
    input("d_mumps", "value"),
    input("d_dengue", "value"),
    input("d_lassa", "value"),
    input("d_avian_opt1", "value"),
    input("d_avian_opt2", "value"),
    input("d_flu_opt1", "value"),
    input("d_flu_opt2", "value"),
    input("d_flu_opt3", "value"),
    input("id-year", "value")
  ),
  function(virus_name, mumps, dengue, lassa,
           avian_opt1, avian_opt2, flu_opt1,
           flu_opt2, flu_opt3, id_year){
    virus_name <- tolower(virus_name)
    if (virus_name %in% c("ebola", "zika", "measles")){
      files <- createPathsFile(virus_name)
    } else if (virus_name == "mumps"){
      files <- createPathsFile(virus_name, mumps)
    } else if (virus_name == "dengue"){
      files <- createPathsFile(virus_name, dengue)
    } else if (virus_name == "lassa"){
      files <- createPathsFile(virus_name, lassa)
    } else if (virus_name == "avian"){
      files <- createPathsFile(
        virus_name, avian_opt1, avian_opt2
      )
    } else if (virus_name == "flu"){
      files <- createPathsFile(
        virus_name, flu_opt1, flu_opt2, flu_opt3
      )
    }
    min_date <- id_year[[1]]
    max_date <- id_year[[2]]
    createMapBubbleYear(
      virus_name, files$stat_file, min_date, max_date
    )
  }
)

app$callback(
  output("id-slicer", "children"),
  list(
    input("d_virus-name", "value"),
    input("d_mumps", "value"),
    input("d_dengue", "value"),
    input("d_lassa", "value"),
    input("d_avian_opt1", "value"),
    input("d_avian_opt2", "value"),
    input("d_flu_opt1", "value"),
    input("d_flu_opt2", "value"),
    input("d_flu_opt3", "value")
  ),
  function(virus_name, mumps, dengue, lassa,
           avian_opt1, avian_opt2, flu_opt1,
           flu_opt2, flu_opt3, id_year){
    virus_name <- tolower(virus_name)
    if (virus_name %in% c("ebola", "zika", "measles")){
      files <- createPathsFile(virus_name)
    } else if (virus_name == "mumps"){
      files <- createPathsFile(virus_name, mumps)
    } else if (virus_name == "dengue"){
      files <- createPathsFile(virus_name, dengue)
    } else if (virus_name == "lassa"){
      files <- createPathsFile(virus_name, lassa)
    } else if (virus_name == "avian"){
      files <- createPathsFile(
        virus_name, avian_opt1, avian_opt2
      )
    } else if (virus_name == "flu"){
      files <- createPathsFile(
        virus_name, flu_opt1, flu_opt2, flu_opt3
      )
    }
    df <- fread(files$stat_file)
    min_year <- min(df[, Year])
    max_year <- max(df[, Year])
    dccRangeSlider(
      id = "id-year",
      min = min_year,
      max = max_year,
      step = 1,
      marks = slicer(min_year, max_year)$marks,
      value = list(min_year, max_year)
    )
  }
)

app$callback(
  output("curve-line-graph", "figure"),
  list(
    input("d_virus-name", "value"),
    input("d_mumps", "value"),
    input("d_dengue", "value"),
    input("d_lassa", "value"),
    input("d_avian_opt1", "value"),
    input("d_avian_opt2", "value"),
    input("d_flu_opt1", "value"),
    input("d_flu_opt2", "value"),
    input("d_flu_opt3", "value"),
    input("id-year", "value")
  ),
  function(virus_name, mumps, dengue, lassa,
           avian_opt1, avian_opt2, flu_opt1,
           flu_opt2, flu_opt3, id_year){
    virus_name <- tolower(virus_name)
    if (virus_name %in% c("ebola", "zika", "measles")){
      files <- createPathsFile(virus_name)
    } else if (virus_name == "mumps"){
      files <- createPathsFile(virus_name, mumps)
    } else if (virus_name == "dengue"){
      files <- createPathsFile(virus_name, dengue)
    } else if (virus_name == "lassa"){
      files <- createPathsFile(virus_name, lassa)
    } else if (virus_name == "avian"){
      files <- createPathsFile(
        virus_name, avian_opt1, avian_opt2
      )
    } else if (virus_name == "flu"){
      files <- createPathsFile(
        virus_name, flu_opt1, flu_opt2, flu_opt3
      )
    }
    min_date <- id_year[[1]]
    max_date <- id_year[[2]]
    createCurveLine(files$stat_file, virus_name, min_date, max_date)
  }
)

app$callback(
  output("id-histo", "figure"),
  list(
    input("d_virus-name", "value"),
    input("d_mumps", "value"),
    input("d_dengue", "value"),
    input("d_lassa", "value"),
    input("d_avian_opt1", "value"),
    input("d_avian_opt2", "value"),
    input("d_flu_opt1", "value"),
    input("d_flu_opt2", "value"),
    input("d_flu_opt3", "value"),
    input("id-year", "value")
  ),
  function(virus_name, mumps, dengue, lassa,
           avian_opt1, avian_opt2, flu_opt1,
           flu_opt2, flu_opt3, id_year){
    virus_name <- tolower(virus_name)
    if (virus_name %in% c("ebola", "zika", "measles")){
      files <- createPathsFile(virus_name)
    } else if (virus_name == "mumps"){
      files <- createPathsFile(virus_name, mumps)
    } else if (virus_name == "dengue"){
      files <- createPathsFile(virus_name, dengue)
    } else if (virus_name == "lassa"){
      files <- createPathsFile(virus_name, lassa)
    } else if (virus_name == "avian"){
      files <- createPathsFile(
        virus_name, avian_opt1, avian_opt2
      )
    } else if (virus_name == "flu"){
      files <- createPathsFile(
        virus_name, flu_opt1, flu_opt2, flu_opt3
      )
    }
    min_date <- id_year[[1]]
    max_date <- id_year[[2]]
    createHistogram(files$stat_file, virus_name, min_date, max_date)
  }
)

app$run_server()
