import pathlib 

axis_template = {
    "showbackground": True,
    "backgroundcolor": "#141414",
    "gridcolor": "rgb(255, 255, 255)",
    "zerolinecolor": "rgb(255, 255, 255)",
}

plot_layout = {
    "title": "",
    "margin": {"t": 0, "b": 0, "l": 0, "r": 0},
    "font": {"size": 12, "color": "white"},
    "showlegend": False,
    "plot_bgcolor": "#141414",
    "paper_bgcolor": "#141414",
    "scene": {
        "xaxis": axis_template,
        "yaxis": axis_template,
        "zaxis": axis_template,
        "aspectratio": {"x": 1, "y": 1.2, "z": 1},
        "camera": {"eye": {"x": 1.25, "y": 1.25, "z": 1.25}},
        "annotations": [],
    },
}

DATA_PATH = pathlib.Path(__file__).parent.joinpath("data").resolve()

default_colorscale = [
    [0, "rgb(12,51,131)"],
    [0.25, "rgb(10,136,186)"],
    [0.5, "rgb(242,211,56)"],
    [0.75, "rgb(242,143,56)"],
    [1, "rgb(217,30,30)"],
]

default_colorscale_index = [ea[1] for ea in default_colorscale]