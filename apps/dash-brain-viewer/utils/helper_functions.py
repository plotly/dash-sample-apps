def add_marker(x, y, z):
    """Create a plotly marker dict."""

    return {
        "x": [x],
        "y": [y],
        "z": [z],
        "mode": "markers",
        "marker": {"size": 25, "line": {"width": 3}},
        "name": "Marker",
        "type": "scatter3d",
        "text": ["Click point to remove annotation"],
    }


def add_annotation(x, y, z):
    """Create plotly annotation dict."""

    return {
        "x": x,
        "y": y,
        "z": z,
        "font": {"color": "black"},
        "bgcolor": "white",
        "borderpad": 5,
        "bordercolor": "black",
        "borderwidth": 1,
        "captureevents": True,
        "ay": -100,
        "arrowcolor": "white",
        "arrowwidth": 2,
        "arrowhead": 0,
        "text": "Click here to annotate<br>(Click point to remove)",
    }


def marker_in_points(points, marker):
    """
    Checks if the marker is in the list of points.

    :params points: a list of dict that contains x, y, z
    :params marker: a dict that contains x, y, z
    :returns: index of the matching marker in list
    """

    for index, point in enumerate(points):
        if (
            point["x"] == marker["x"]
            and point["y"] == marker["y"]
            and point["z"] == marker["z"]
        ):
            return index
    return None
