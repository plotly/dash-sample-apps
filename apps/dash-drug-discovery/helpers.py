import dash_html_components as html


def make_dash_table(selection, df):
    """ Return a dash defintion of an HTML table from a Pandas dataframe. """

    df_subset = df.loc[df["NAME"].isin(selection)]
    table = []

    for index, row in df_subset.iterrows():
        rows = []
        rows.append(html.Td([row["NAME"]]))
        rows.append(html.Td([html.Img(src=row["IMG_URL"])]))
        rows.append(html.Td([row["FORM"]]))
        rows.append(
            html.Td([html.A(href=row["PAGE"], children="Datasheet", target="_blank")])
        )
        table.append(html.Tr(rows))

    return table


def _add_markers(figure_data, molecules, plot_type="scatter3d"):
    """
    Add markers on the plot graph.

    :params figure_data: the graph data
    :params molecules: list of selected molecules
    :params plot_type: plot type (scatter3d, histogram2d, scatter)
    :returns: plotly graph trace list
    """

    drug_data = figure_data[0]
    list_of_drugs = drug_data["text"]

    # get the axis index for each drug
    indices = [index for index, value in enumerate(list_of_drugs) if value in molecules]

    if plot_type == "histogram2d":
        plot_type = "scatter"

    traces = []
    for point_number in indices:
        trace = {
            "x": [drug_data["x"][point_number]],
            "y": [drug_data["y"][point_number]],
            "marker": {"color": "red", "size": 16, "opacity": 0.6, "symbol": "cross"},
            "type": plot_type,
        }
        if plot_type == "scatter3d":
            trace["z"] = [drug_data["z"][point_number]]
        traces.append(trace)
    return traces


def _create_axis(axis_type, variation="Linear", title=None):
    """
    Creates a 2d or 3d axis.

    :params axis_type: 2d or 3d axis
    :params variation: axis type (log, line, linear, etc)
    :parmas title: axis title
    :returns: plotly axis dictionnary
    """

    if axis_type not in ["3d", "2d"]:
        return None

    default_style = {
        "background": "rgb(230, 230, 230)",
        "gridcolor": "rgb(255, 255, 255)",
        "zerolinecolor": "rgb(255, 255, 255)",
    }

    if axis_type == "3d":
        return {
            "showbackground": True,
            "backgroundcolor": default_style["background"],
            "gridcolor": default_style["gridcolor"],
            "title": title,
            "type": variation,
            "zerolinecolor": default_style["zerolinecolor"],
        }

    if axis_type == "2d":
        return {
            "xgap": 10,
            "ygap": 10,
            "backgroundcolor": default_style["background"],
            "gridcolor": default_style["gridcolor"],
            "title": title,
            "zerolinecolor": default_style["zerolinecolor"],
            "color": "#444",
        }


def _black_out_axis(axis):
    axis["showgrid"] = False
    axis["zeroline"] = False
    axis["color"] = "white"
    return axis


def _create_layout(layout_type, xlabel, ylabel):
    """ Return dash plot layout. """

    base_layout = {
        "font": {"family": "Raleway"},
        "hovermode": "closest",
        "margin": {"r": 20, "t": 0, "l": 0, "b": 0},
        "showlegend": False,
    }

    if layout_type == "scatter3d":
        base_layout["scene"] = {
            "xaxis": _create_axis(axis_type="3d", title=xlabel),
            "yaxis": _create_axis(axis_type="3d", title=ylabel),
            "zaxis": _create_axis(axis_type="3d", title=xlabel, variation="log"),
            "camera": {
                "up": {"x": 0, "y": 0, "z": 1},
                "center": {"x": 0, "y": 0, "z": 0},
                "eye": {"x": 0.08, "y": 2.2, "z": 0.08},
            },
        }

    elif layout_type == "histogram2d":
        base_layout["xaxis"] = _black_out_axis(
            _create_axis(axis_type="2d", title=xlabel)
        )
        base_layout["yaxis"] = _black_out_axis(
            _create_axis(axis_type="2d", title=ylabel)
        )
        base_layout["plot_bgcolor"] = "black"
        base_layout["paper_bgcolor"] = "black"
        base_layout["font"]["color"] = "white"

    elif layout_type == "scatter":
        base_layout["xaxis"] = _create_axis(axis_type="2d", title=xlabel)
        base_layout["yaxis"] = _create_axis(axis_type="2d", title=ylabel)
        base_layout["plot_bgcolor"] = "rgb(230, 230, 230)"
        base_layout["paper_bgcolor"] = "rgb(230, 230, 230)"

    return base_layout


def create_plot(
    x,
    y,
    z,
    size,
    color,
    name,
    xlabel="LogP",
    ylabel="pkA",
    zlabel="Solubility (mg/ml)",
    plot_type="scatter3d",
    markers=[],
):

    colorscale = [
        [0, "rgb(244,236,21)"],
        [0.3, "rgb(249,210,41)"],
        [0.4, "rgb(134,191,118)"],
        [0.5, "rgb(37,180,167)"],
        [0.65, "rgb(17,123,215)"],
        [1, "rgb(54,50,153)"],
    ]

    data = [
        {
            "x": x,
            "y": y,
            "z": z,
            "mode": "markers",
            "marker": {
                "colorscale": colorscale,
                "colorbar": {"title": "Molecular<br>Weight"},
                "line": {"color": "#444"},
                "reversescale": True,
                "sizeref": 45,
                "sizemode": "diameter",
                "opacity": 0.7,
                "size": size,
                "color": color,
            },
            "text": name,
            "type": plot_type,
        }
    ]

    if plot_type in ["histogram2d", "scatter"]:
        del data[0]["z"]

    if plot_type == "histogram2d":
        # Scatter plot overlay on 2d Histogram
        data[0]["type"] = "scatter"
        data.append(
            {
                "x": x,
                "y": y,
                "type": "histogram2d",
                "colorscale": "Greys",
                "showscale": False,
            }
        )

    layout = _create_layout(plot_type, xlabel, ylabel)

    if len(markers) > 0:
        data = data + _add_markers(data, markers, plot_type=plot_type)

    return {"data": data, "layout": layout}
