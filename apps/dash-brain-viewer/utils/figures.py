from utils.helper_functions import add_marker, add_annotation, marker_in_points
from dash.exceptions import PreventUpdate


def brain_graph_handler(click_data, val, colorscale, figure, current_anno):
    """Listener on colorscale, option picker, and graph on click to update the graph."""

    # new option select
    if figure["data"][0]["name"] != val:
        figure["data"] = create_mesh_data(val)
        figure["layout"] = plot_layout
        cs = [[i / (len(colorscale) - 1), rgb] for i, rgb in enumerate(colorscale)]
        figure["data"][0]["colorscale"] = cs
        return figure

    # modify graph markers
    if click_data is not None and "points" in click_data:

        y_value = click_data["points"][0]["y"]
        x_value = click_data["points"][0]["x"]
        z_value = click_data["points"][0]["z"]

        marker = add_marker(x_value, y_value, z_value)
        point_index = marker_in_points(figure["data"], marker)

        # delete graph markers
        if len(figure["data"]) > 1 and point_index is not None:

            figure["data"].pop(point_index)
            anno_index_offset = 2 if val == "mouse" else 1
            try:
                figure["layout"]["scene"]["annotations"].pop(
                    point_index - anno_index_offset
                )
            except Exception as error:
                print(error)
                pass

        # append graph markers
        else:

            # iterate through the store annotations and save it into figure data
            if current_anno is not None:
                for index, annotations in enumerate(
                    figure["layout"]["scene"]["annotations"]
                ):
                    for key in current_anno.keys():
                        if str(index) in key:
                            figure["layout"]["scene"]["annotations"][index][
                                "text"
                            ] = current_anno[key]

            figure["data"].append(marker)
            figure["layout"]["scene"]["annotations"].append(
                add_annotation(x_value, y_value, z_value)
            )

    cs = [[i / (len(colorscale) - 1), rgb] for i, rgb in enumerate(colorscale)]
    figure["data"][0]["colorscale"] = cs

    return figure


def save_annotations(relayout_data, current_data):
    """Update the annotations in the dcc store."""

    if relayout_data is None:
        raise PreventUpdate

    if current_data is None:
        return {}

    for key in relayout_data.keys():

        # to determine if the relayout has to do with annotations
        if "scene.annotations" in key:
            current_data[key] = relayout_data[key]

    return current_data
