source("helperfunctions.R")

# Data processing for surf_reg_model_both
pts <- read_mniobj("data/surf_reg_model_both.obj")[1] # ~ Takes around 22 sec
tri <- read_mniobj("data/surf_reg_model_both.obj")[2] # ~ Takes around 22 sec

intensities <- lapply(as.list(read.table('data/aal_atlas.txt')), as.list)[[1]]

default_colorscale <- list(
  list(0, '#0c3383'),
  list(0.25, '#0a88ba'),
  list(0.5, '#f2d338'),
  list(0.75, '#f28f38'),
  list(1, '#d91e1e')
)

initialGraphdata <- plotly_triangular_mesh(
  pts,
  tri,
  intensities,
  colorscale = default_colorscale,
  flatshading = FALSE,
  showscale = FALSE,
  reversescale = FALSE,
  plot_edges = FALSE
)

initialGraphdata <- initialGraphdata[[1]]
initialGraphdata['name'] <- 'human_atlas'
initialGraphdata[c(2, 3, 4, 8, 9, 10)] <-
  lapply(initialGraphdata[c(2, 3, 4, 8, 9, 10)], unlist)

# Data processing for realct
pts <- read_mniobj('data/realct.obj')[1]
tri <- read_mniobj('data/realct.obj')[2]
intensities <- lapply(as.list(read.table('data/realct.txt')), as.list)[[1]]
realct <- plotly_triangular_mesh(
  pts,
  tri,
  intensities,
  colorscale = DEFAULT_COLORSCALE,
  flatshading = FALSE,
  showscale = FALSE,
  reversescale = FALSE,
  plot_edges = FALSE
)

realct <- realct[[1]]
realct[c(2,3,4,8,9,10)] <- lapply(realct[c(2,3,4,8,9,10)], unlist)

# Data processing for surf_reg_model_both
pts <- read_mniobj('data/surf_reg_model_both.obj')[1]
tri <- read_mniobj('data/surf_reg_model_both.obj')[2]
intensities <- lapply(as.list(read.table('data/aal_atlas.txt')),as.list)[[1]]
surfreg <- plotly_triangular_mesh(
  pts,
  tri,
  intensities,
  colorscale = DEFAULT_COLORSCALE,
  flatshading = FALSE,
  showscale = FALSE,
  reversescale = FALSE,
  plot_edges = FALSE
)

surfreg <- surfreg[[1]]
surfreg[c(2,3,4,8,9,10)] <- lapply(surfreg[c(2,3,4,8,9,10)], unlist)

# Data processing for mouse_surf
pts <- read_mniobj('data/mouse_surf.obj')[1]
tri <- read_mniobj('data/mouse_surf.obj')[2]
intensities <- lapply(as.list(read.table('data/mouse_map.txt')),as.list)[[1]]
mouse_map <- plotly_triangular_mesh(
  pts,tri,intensities,
  colorscale=DEFAULT_COLORSCALE, flatshading=FALSE,
  showscale=FALSE, reversescale=FALSE, plot_edges=FALSE)

mouse_map[[1]][c(2,3,4,8,9,10)] = lapply(mouse_map[[1]][c(2,3,4,8,9,10)], unlist)

# Data processing for mouse_brain_outline
pts <- read_mniobj('data/mouse_brain_outline.obj')[1]
tri <- read_mniobj('data/mouse_brain_outline.obj')[2]
outer_mesh <- plotly_triangular_mesh(
  pts,
  tri,
  intensities,
  colorscale = DEFAULT_COLORSCALE,
  flatshading = FALSE,
  showscale = FALSE,
  reversescale = FALSE,
  plot_edges = FALSE
)[1]

outer_mesh[[1]]['opacity'] <- 0.5
outer_mesh[[1]]['colorscale'] <- 'Greys'
outer_mesh[[1]][c(2,3,4,8,9,10)] <- lapply(outer_mesh[[1]][c(2,3,4,8,9,10)], unlist)
mouse_map[[2]] <- outer_mesh[[1]]

# Saving to file
save(
  'initialGraphdata',
  'realct',
  'surfreg',
  'mouse_map',
  file = 'alldata.RData',
  compress = TRUE,
  version = 2
)
