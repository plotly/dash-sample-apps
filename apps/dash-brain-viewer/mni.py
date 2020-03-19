import pathlib
import numpy as np


DATA_PATH = pathlib.Path(__file__).parent.joinpath("data").resolve()

default_colorscale = [
    [0, "rgb(12,51,131)"],
    [0.25, "rgb(10,136,186)"],
    [0.5, "rgb(242,211,56)"],
    [0.75, "rgb(242,143,56)"],
    [1, "rgb(217,30,30)"],
]


def read_mniobj(file):
    """
    Parses an obj file.
    
    :params file: file name in data folder
    :returns: a tuple
    """

    def triangulate_polygons(list_vertex_indices):
        for k in range(0, len(list_vertex_indices), 3):
            yield list_vertex_indices[k : k + 3]

    with open(DATA_PATH.joinpath(file)) as fp:
        num_vertices = 0
        matrix_vertices = []
        k = 0
        list_indices = []

        for i, line in enumerate(fp):
            if i == 0:
                num_vertices = int(line.split()[6])
                matrix_vertices = np.zeros([num_vertices, 3])
            elif i <= num_vertices:
                matrix_vertices[i - 1] = list(map(float, line.split()))
            elif i > 2 * num_vertices + 5:
                if not line.strip():
                    k = 1
                elif k == 1:
                    list_indices.extend(line.split())

    list_indices = [int(i) for i in list_indices]
    faces = np.array(list(triangulate_polygons(list_indices)))
    return matrix_vertices, faces


def plotly_triangular_mesh(
    vertices,
    faces,
    intensities=None,
    colorscale="Viridis",
    flatshading=False,
    showscale=False,
    plot_edges=False,
):

    x, y, z = vertices.T
    I, J, K = faces.T

    if intensities is None:
        intensities = z

    mesh = {
        "type": "mesh3d",
        "x": x,
        "y": y,
        "z": z,
        "colorscale": colorscale,
        "intensity": intensities,
        "flatshading": flatshading,
        "i": I,
        "j": J,
        "k": K,
        "name": "",
        "showscale": showscale,
        "lighting": {
            "ambient": 0.18,
            "diffuse": 1,
            "fresnel": 0.1,
            "specular": 1,
            "roughness": 0.1,
            "facenormalsepsilon": 1e-6,
            "vertexnormalsepsilon": 1e-12,
        },
        "lightposition": {"x": 100, "y": 200, "z": 0},
    }

    if showscale:
        mesh["colorbar"] = {"thickness": 20, "ticklen": 4, "len": 0.75}

    if plot_edges is False:
        return [mesh]

    lines = create_plot_edges_lines(vertices, faces)
    return [mesh, lines]


def create_plot_edges_lines(vertices, faces):
    tri_vertices = vertices[faces]
    Xe = []
    Ye = []
    Ze = []
    for T in tri_vertices:
        Xe += [T[k % 3][0] for k in range(4)] + [None]
        Ye += [T[k % 3][1] for k in range(4)] + [None]
        Ze += [T[k % 3][2] for k in range(4)] + [None]

    # define the lines to be plotted
    lines = {
        "type": "scatter3d",
        "x": Xe,
        "y": Ye,
        "z": Ze,
        "mode": "lines",
        "name": "",
        "line": {"color": "rgb(70,70,70)", "width": 1},
    }
    return lines


def create_mesh_data(option):

    data = []
    if option == "human":
        vertices, faces = read_mniobj("realct.obj")
        intensities = np.loadtxt(DATA_PATH.joinpath("realct.txt"))
    elif option == "human_atlas":
        vertices, faces = read_mniobj("surf_reg_model_both.obj")
        intensities = np.loadtxt(DATA_PATH.joinpath("aal_atlas.txt"))
    elif option == "mouse":
        vertices, faces = read_mniobj("mouse_surf.obj")
        intensities = np.loadtxt(DATA_PATH.joinpath("mouse_map.txt"))
    else:
        raise ValueError

    data = plotly_triangular_mesh(
        vertices, faces, intensities, colorscale=default_colorscale
    )

    if option == "mouse":
        vertices, faces = read_mniobj("mouse_brain_outline.obj")
        outer_mesh = plotly_triangular_mesh(vertices, faces)[0]
        outer_mesh["opacity"] = 0.5
        outer_mesh["colorscale"] = "Greys"
        data.append(outer_mesh)

    data[0]["name"] = option
    return data
