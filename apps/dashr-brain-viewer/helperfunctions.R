# Helper Functions ------------------------------

read_mniobj <- function(file) {

  triangulate_polygons <- function(list_vertex_indices) {
    j = 1
    otherlist = list()
    for (k in seq(1, length(list_vertex_indices), 3)) {
      otherlist[[j]] = list_vertex_indices[k:(k + 2)]
      j = j + 1
    }
    return(otherlist)
  }

  fp = readLines(con = file)
  n_vert <- list()
  n_poly <- list()
  k = 0
  j = 1
  list_indices <- list(list())
  fp <- as.list(fp)

  for (i in as.numeric(labels(fp))) {
    if (i == 1) {
      n_vert = as.numeric(as.list(strsplit(trimws(fp[i]), " ")[[1]])[7])
      vertices <- matrix(0, nrow = n_vert, ncol = 3)
    } else if (i <= n_vert) {
      vertices[i - 1,] <-
        unlist(lapply(strsplit(trimws(fp[i]), " ")[[1]], as.numeric))
    } else if (i > 2 * (n_vert) + 6) {
      if (lapply(fp[i], is.empty) == TRUE) {
        k = 1
      } else if (k == 1) {
        list_indices[[j]] <- strsplit(trimws(fp[i]), " ")
        j = j + 1
      }
    }
  }

  list_indices <-  lapply(as.list(unlist(list_indices)),
                          function(x) { as.integer(as.list(x))})

  faces <- as.array(triangulate_polygons(list_indices))

  return(list(vertices, faces))
}


plotly_triangular_mesh <- function(
  vertices,
  faces,
  intensities = NULL,
  colorscale = NULL,
  flatshading = FALSE,
  showscale = FALSE,
  reversescale = FALSE,
  plot_edges = FALSE) {

  I = list()
  J = list()
  K = list()

  for (i in 1:length(faces[[1]])) {
    I[i] = faces[[1]][[i]][1]
    J[i] = faces[[1]][[i]][2]
    K[i] = faces[[1]][[i]][3]
  }

  X = list()
  Y = list()
  Z = list()
  X = as.list(vertices[[1]][,1])
  Y = as.list(vertices[[1]][,2])
  Z = as.list(vertices[[1]][,3])

  if (length(intensities) == 0) {
    intensities <- standard_intensity(x, y, z)
  }

  if (!is.null(intensities[['__call__']])) {
    intensity <- intensities(x, y, z)
  } else if (class(intensities) == 'list' |
             class(intensities) == 'array') {
    intensity <- intensities
  } else {
    sprintf("intensities can be either a function or a list, np.array")
  }

  mesh <- list(
    type = 'isosurface',
    x = X,
    y = Y,
    z = Z,
    colorscale = colorscale,
    intensity = as.list(intensities),
    flatshading = flatshading,
    i = I,
    j = J,
    k = K,
    name = '',
    showscale = showscale,
    lighting = list(
      ambient = 0.18,
      diffuse = 1,
      fresnel =  0.1,
      specular = 1,
      roughness = 0.1,
      facenormalsepsilon = 1e-6,
      vertexnormalsepsilon = 1e-12
    ),
    lightposition = list(x = 100, y = 200, z = 0)
  )

  if (isTRUE(showscale)) {
    mesh4 <- list(colorbar = list(
      thickness = 20,
      ticklen = 4,
      len = 0.75
    ))
    mesh <- modifyList(mesh, mesh4)
  }

  if (length(plot_edges) == 0) {
    return(list(mesh))
  } else {
    tri_vertices = vertices[['faces']]
    Xe = list()
    Ye = list()
    Ze = list()

    for (T in tri_vertices) {
      for (k in 1:4) {
        Xe[[k]] <- Xe + list(T[[k %% 3]][1])
        Ye[[k]] <- Xe + list(T[[k %% 3]][2])
        Ze[[k]] <- Xe + list(T[[k %% 3]][3])

        lines <- list(
          type = 'scatter3d',
          x = Xe,
          y = Ye,
          z = Ze,
          mode = 'lines',
          name = '',
          line = list(color = '#464646', width = 1)
        )
      }
    }
  }
  return(list(mesh,lines))
}
