pts <- read_mniobj("surf_reg_model_both.obj")[1]
tri <- read_mniobj("surf_reg_model_both.obj")[2]
intensities <- lapply(as.list(read.table('aal_atlas.txt')),as.list)[[1]]

data=plotly_triangular_mesh(pts, tri, intensities,
                            colorscale=DEFAULT_COLORSCALE, flatshading=FALSE,
                            showscale=FALSE, reversescale=FALSE, plot_edges=FALSE)
data = data[[1]]
data['name'] = 'human_atlas'
data[c(2,3,4,8,9,10)] = lapply(data[c(2,3,4,8,9,10)], unlist)


saveRDS(data,file = 'initialGraphdata.RData')

pts = read_mniobj('realct.obj')[1]
tri = read_mniobj('realct.obj')[2]
intensities <- lapply(as.list(read.table('realct.txt')),as.list)[[1]]
data = plotly_triangular_mesh(pts,tri,intensities, 
                              colorscale=DEFAULT_COLORSCALE, flatshading=FALSE,
                              showscale=FALSE, reversescale=FALSE, plot_edges=FALSE)

data = data[[1]]
data[c(2,3,4,8,9,10)] = lapply(data[c(2,3,4,8,9,10)], unlist)

saveRDS(data,file = 'realct.RData')

pts = read_mniobj('surf_reg_model_both.obj')[1]
tri = read_mniobj('surf_reg_model_both.obj')[2]
intensities <- lapply(as.list(read.table('aal_atlas.txt')),as.list)[[1]]
data = plotly_triangular_mesh(pts,tri,intensities, 
                              colorscale=DEFAULT_COLORSCALE, flatshading=FALSE,
                              showscale=FALSE, reversescale=FALSE, plot_edges=FALSE)

data = data[[1]]
data[c(2,3,4,8,9,10)] = lapply(data[c(2,3,4,8,9,10)], unlist)

saveRDS(data,file = 'surfreg.RData')

pts = read_mniobj('mouse_surf.obj')[1]
tri = read_mniobj('mouse_surf.obj')[2]
intensities <- lapply(as.list(read.table('mouse_map.txt')),as.list)[[1]]
data = plotly_triangular_mesh(pts,tri,intensities, 
                              colorscale=DEFAULT_COLORSCALE, flatshading=FALSE,
                              showscale=FALSE, reversescale=FALSE, plot_edges=FALSE)

#data = data[[1]]
data[[1]][c(2,3,4,8,9,10)] = lapply(data[[1]][c(2,3,4,8,9,10)], unlist)

pts = read_mniobj('mouse_brain_outline.obj')[1]
tri = read_mniobj('mouse_brain_outline.obj')[2]
outer_mesh = plotly_triangular_mesh(pts,tri,intensities, 
                                    colorscale=DEFAULT_COLORSCALE, flatshading=FALSE,
                                    showscale=FALSE, reversescale=FALSE, plot_edges=FALSE)[1]
outer_mesh[[1]]['opacity'] = 0.5
outer_mesh[[1]]['colorscale'] = 'Greys'
outer_mesh[[1]][c(2,3,4,8,9,10)] = lapply(outer_mesh[[1]][c(2,3,4,8,9,10)], unlist)
data[[2]] = outer_mesh[[1]]

saveRDS(data,file = 'mouse.RData')
