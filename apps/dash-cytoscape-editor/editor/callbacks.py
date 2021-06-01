# pylint: disable=W0612,R1705
import json
from colour import Color

from dash.dependencies import Input, Output, State

from .constants import (
    ARROW_POSITIONS,
    LABEL_ELEMENT_TYPES_ALL,
    LABEL_ELEMENT_TYPES,
    ELEMENTS,
)


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def get_ids(elements):
    ids = []
    for n in range(len(elements)):
        curr_id = elements[n].get("data").get("id")
        if curr_id:
            ids.append(curr_id)
    return ids


def validate_positive(value):
    return min(0, value)


def validate_color(color, default="#999999"):
    """
    Check if a color is valid, if so returns the color, else return a default color
    :param color: The color to validate
    :param default: The default color
    :return: A string representing a color
    """
    if not color:
        return default

    try:
        # Converting 'deep sky blue' to 'deepskyblue'
        color = color.replace(" ", "")

        if color.startswith("rgb"):
            values = color.replace("rgb(", "").replace(")", "").split(",")

            if len(values) == 3 and all(0 <= int(v) <= 255 for v in values):
                return color

            return default

        Color(color)
        # if everything goes fine then return True
        return color
    except:  # noqa
        return default


def validate_px_percentage(value, default="0px"):
    if not value:
        return default
    elif "px" in value and is_float(value.replace("px", "")):
        return value
    elif "%" in value and is_float(value.replace("%", "")):
        return value
    else:
        return default


def assign_callbacks(app):
    # ############################## HIDING ###################################
    for n in range(1, 17):

        @app.callback(
            Output(f"div-pie-slice-{n}", "style"),
            [Input("dropdown-pie-slice-selected", "value")],
            [State(f"div-pie-slice-{n}", "id")],
        )
        def hide_div_pie_slice(current_slice_selected, div_id):
            if current_slice_selected != div_id:
                return {"display": "none"}
            else:
                return {"display": "block"}

    for pos in ARROW_POSITIONS:

        @app.callback(
            Output(f"div-arrow-position-{pos}", "style"),
            [Input("dropdown-arrow-position", "value")],
            [State(f"div-arrow-position-{pos}", "id")],
        )
        def hide_div_arrow_position(current_pos_selected, div_id):
            if current_pos_selected != div_id:
                return {"display": "none"}
            else:
                return {"display": "block"}

    for element in LABEL_ELEMENT_TYPES:

        @app.callback(
            Output(f"div-label-{element}", "style"),
            [Input("dropdown-select-element-label-styling", "value")],
            [State(f"div-label-{element}", "id")],
        )
        def hide_div_label_element(current_element_selected, div_id):
            if current_element_selected != div_id:
                return {"display": "none"}
            else:
                return {"display": "block"}

        @app.callback(
            Output(f"div-text-wrapping-{element}", "style"),
            [Input("dropdown-select-element-text-wrapping", "value")],
            [State(f"div-text-wrapping-{element}", "id")],
        )
        def hide_div_text_wrapping(current_element_selected, div_id):
            if current_element_selected != div_id:
                return {"display": "none"}
            else:
                return {"display": "block"}

    for element in LABEL_ELEMENT_TYPES_ALL:

        @app.callback(
            Output(f"div-text-margins-{element}", "style"),
            [Input("dropdown-select-element-text-margins", "value")],
            [State(f"div-text-margins-{element}", "id")],
        )
        def hide_div_text_margins(current_element_selected, div_id):
            if current_element_selected != div_id:
                return {"display": "none"}
            else:
                return {"display": "block"}

    @app.callback(
        Output("div-display-stylesheet-json", "children"),
        [Input("cytoscape", "stylesheet")],
    )
    def update_json_stylesheet_output(stylesheet):
        return json.dumps(stylesheet, indent=2)

    @app.callback(
        Output("div-display-elements-json", "children"),
        [Input("cytoscape", "elements")],
    )
    def update_json_elements_output(stylesheet):
        return json.dumps(stylesheet, indent=2)

    # ############################## STORING ##################################
    @app.callback(
        Output("div-storage-pie-background-color", "children"),
        [Input(f"input-pie-{n}-background-color", "value") for n in range(1, 17)],
    )
    def update_pie_color_storage(*args):
        args = [validate_color(color) for color in args]
        return json.dumps(
            dict(zip([f"pie-{i}-background-color" for i in range(1, 17)], args))
        )

    @app.callback(
        Output("div-storage-pie-background-size", "children"),
        [Input(f"input-pie-{n}-background-size", "value") for n in range(1, 17)],
    )
    def update_pie_size_storage(*args):
        return json.dumps(
            dict(zip([f"pie-{i}-background-size" for i in range(1, 17)], args))
        )

    @app.callback(
        Output("div-storage-pie-background-opacity", "children"),
        [Input(f"slider-pie-{n}-background-opacity", "value") for n in range(1, 17)],
    )
    def update_pie_opacity_storage(*args):
        return json.dumps(
            dict(zip([f"pie-{i}-background-opacity" for i in range(1, 17)], args))
        )

    @app.callback(
        Output("div-storage-arrow-color", "children"),
        [Input(f"input-{pos}-arrow-color", "value") for pos in ARROW_POSITIONS],
    )
    def update_arrow_color_storage(*args):
        args = [validate_color(color) for color in args]
        return json.dumps(
            dict(zip([f"{pos}-arrow-color" for pos in ARROW_POSITIONS], args))
        )

    @app.callback(
        Output("div-storage-arrow-shape", "children"),
        [Input(f"dropdown-{pos}-arrow-shape", "value") for pos in ARROW_POSITIONS],
    )
    def update_arrow_shape_storage(*args):
        return json.dumps(
            dict(zip([f"{pos}-arrow-shape" for pos in ARROW_POSITIONS], args))
        )

    @app.callback(
        Output("div-storage-arrow-fill", "children"),
        [Input(f"radio-{pos}-arrow-fill", "value") for pos in ARROW_POSITIONS],
    )
    def update_arrow_fill_storage(*args):
        return json.dumps(
            dict(zip([f"{pos}-arrow-fill" for pos in ARROW_POSITIONS], args))
        )

    # ############################## DISABLING ################################
    @app.callback(
        Output("input-background-image-height", "disabled"),
        [Input("radio-background-image-fit", "value")],
    )
    def disable_background_image_height(value):
        return value != "none"

    @app.callback(
        Output("input-background-image-width", "disabled"),
        [Input("radio-background-image-fit", "value")],
    )
    def disable_background_image_width(value):
        return value != "none"

    for side in ["source", "target"]:

        @app.callback(
            Output(f"input-{side}-endpoint-width", "disabled"),
            [Input(f"dropdown-{side}-endpoint-type", "value")],
        )
        def disable_side_endpoint_width(value):
            return value != "other"

        @app.callback(
            Output(f"input-{side}-endpoint-height", "disabled"),
            [Input(f"dropdown-{side}-endpoint-type", "value")],
        )
        def disable_side_endpoint_height(value):
            return value != "other"

    # ############################## CYTOSCAPE ################################
    @app.callback(
        Output("cytoscape", "elements"),
        [Input("dropdown-select-element-list", "value")],
    )
    def update_elements(dataset):
        return ELEMENTS[dataset]

    @app.callback(Output("cytoscape", "layout"), [Input("dropdown-layout", "value")])
    def update_layout(name):
        return {"name": name}

    @app.callback(
        Output("cytoscape", "stylesheet"),
        [
            Input(component, "value")
            for component in [
                # Node Body
                "input-node-content",
                "input-node-width",
                "input-node-height",
                "dropdown-node-shape",
                "input-node-color",
                "slider-node-opacity",
                "slider-node-blacken",
                "input-node-border-width",
                "dropdown-node-border-style",
                "input-node-border-color",
                "slider-node-border-opacity",
                "input-node-padding",
                "dropdown-node-padding-relative-to",
                "radio-node-compound-sizing",
                "input-node-compound-min-width",
                "input-node-compound-min-width-bias-left",
                "input-node-compound-min-width-bias-right",
                "input-node-compound-min-height",
                "input-node-compound-min-height-bias-top",
                "input-node-compound-min-height-bias-bottom",
                # Background Image
                "radio-use-background-image",
                "input-background-image-url",
                "radio-background-image-crossorigin",
                "slider-background-image-opacity",
                "input-background-image-width",
                "input-background-image-height",
                "radio-background-image-fit",
                "input-background-position-x",
                "input-background-position-y",
                "radio-background-width-relative-to",
                "radio-background-height-relative-to",
                "radio-use-pie-chart",
                "input-pie-size",
            ]
        ]
        + [
            Input(div, "children")
            for div in [
                "div-storage-pie-background-color",
                "div-storage-pie-background-size",
                "div-storage-pie-background-opacity",
            ]
        ]
        + [
            Input(component, "value")
            for component in [
                "input-edge-line-width",
                "dropdown-edge-curve-style",
                "input-edge-line-color",
                "radio-edge-line-style",
                "input-edge-loop-direction",
                "input-edge-loop-sweep",
                "radio-use-edge-arrow",
            ]
        ]
        + [
            Input(div, "children")
            for div in [
                "div-storage-arrow-color",
                "div-storage-arrow-shape",
                "div-storage-arrow-fill",
            ]
        ]
        + [
            Input(component, "value")
            for component in [
                "input-arrow-scale",
                "radio-use-edge-endpoints",
                "dropdown-source-endpoint-type",
                "input-source-endpoint-width",
                "input-source-endpoint-height",
                "dropdown-target-endpoint-type",
                "input-target-endpoint-width",
                "input-target-endpoint-height",
                "input-source-distance-from-node",
                "input-target-distance-from-node",
                # Components for Labels
                "radio-use-labels",
                "input-node-label",
                "input-edge-label",
                "input-edge-source-label",
                "input-edge-target-label",
                # Label Font Styling
                "input-node-label-color",
                "slider-node-label-text-opacity",
                "input-node-label-font-family",
                "input-node-label-font-size",
                "dropdown-node-label-font-style",
                "dropdown-node-label-font-weight",
                "dropdown-node-label-text-transform",
                "input-edge-label-color",
                "slider-edge-label-text-opacity",
                "input-edge-label-font-family",
                "input-edge-label-font-size",
                "dropdown-edge-label-font-style",
                "dropdown-edge-label-font-weight",
                "dropdown-edge-label-text-transform",
                # Label Text Wrapping
                "radio-node-label-text-wrap",
                "input-node-label-text-max-width",
                "radio-edge-label-text-wrap",
                "input-edge-label-text-max-width",
                # Label Alignment
                "radio-label-text-halign",
                "radio-label-text-valign",
                "input-label-source-text-offset",
                "input-label-target-text-offset",
                # Text Margins
                "input-node-text-margin-x",
                "input-node-text-margin-y",
                "input-edge-text-margin-x",
                "input-edge-text-margin-y",
                "input-source-text-margin-x",
                "input-source-text-margin-y",
                "input-target-text-margin-x",
                "input-target-text-margin-y",
            ]
        ],
    )
    def update_stylesheet(
        node_content,
        node_width,
        node_height,
        node_shape,
        node_color,
        node_opacity,
        node_blacken,
        node_border_width,
        node_border_style,
        node_border_color,
        node_border_opacity,
        node_padding,
        node_padding_relative_to,
        node_compound_sizing,
        node_compound_min_width,
        node_compound_min_width_bias_left,
        node_compound_min_width_bias_right,
        node_compound_min_height,
        node_compound_min_height_bias_top,
        node_compound_min_height_bias_bottom,
        use_background_image,
        background_image_url,
        background_image_crossorigin,
        background_image_opacity,
        background_image_width,
        background_image_height,
        background_image_fit,
        background_position_x,
        background_position_y,
        background_width_relative_to,
        background_height_relative_to,
        use_pie_chart,
        pie_size,
        storage_pie_background_color,
        storage_pie_background_size,
        storage_pie_background_opacity,
        edge_line_width,
        edge_curve_style,
        edge_line_color,
        edge_line_style,
        edge_loop_direction,
        edge_loop_sweep,
        use_edge_arrow,
        storage_arrow_color,
        storage_arrow_shape,
        storage_arrow_fill,
        arrow_scale,
        use_edge_endpoints,
        source_endpoint_type,
        source_endpoint_width,
        source_endpoint_height,
        target_endpoint_type,
        target_endpoint_width,
        target_endpoint_height,
        source_distance_from_node,
        target_distance_from_node,
        use_labels,
        node_label,
        edge_label,
        edge_source_label,
        edge_target_label,
        node_label_color,
        node_label_text_opacity,
        node_label_font_family,
        node_label_font_size,
        node_label_font_style,
        node_label_font_weight,
        node_label_text_transform,
        edge_label_color,
        edge_label_text_opacity,
        edge_label_font_family,
        edge_label_font_size,
        edge_label_font_style,
        edge_label_font_weight,
        edge_label_text_transform,
        node_label_text_wrap,
        node_label_text_max_width,
        edge_label_text_wrap,
        edge_label_text_max_width,
        label_text_halign,
        label_text_valign,
        label_source_text_offset,
        label_target_text_offset,
        node_text_margin_x,
        node_text_margin_y,
        edge_text_margin_x,
        edge_text_margin_y,
        source_text_margin_x,
        source_text_margin_y,
        target_text_margin_x,
        target_text_margin_y,
    ):
        def update_style(stylesheet, selector, addition):
            for style in stylesheet:
                if style["selector"] == selector:
                    style["style"].update(addition)

        # Validating Input
        node_color = validate_color(node_color)
        node_border_color = validate_color(node_border_color)
        node_padding = validate_px_percentage(node_padding)
        background_position_x = validate_px_percentage(background_position_x)
        background_position_y = validate_px_percentage(background_position_y)
        pie_size = validate_px_percentage(pie_size, default="100%")
        edge_line_color = validate_color(edge_line_color)

        stylesheet = [
            {
                "selector": "node",
                "style": {
                    "content": node_content,
                    "width": node_width,
                    "height": node_height,
                    "background-color": node_color,
                    "background-blacken": node_blacken,
                    "background-opacity": node_opacity,
                    "shape": node_shape,
                    "border-width": node_border_width,
                    "border-style": node_border_style,
                    "border-color": node_border_color,
                    "border-opacity": node_border_opacity,
                    "padding": node_padding,
                    "padding-relative-to": node_padding_relative_to,
                    "compound-sizing-wrt-labels": node_compound_sizing,
                    "min-width": node_compound_min_width,
                    "min-width-bias-left": node_compound_min_width_bias_left,
                    "min-width-bias-right": node_compound_min_width_bias_right,
                    "min-height": node_compound_min_height,
                    "min-height-bias-top": node_compound_min_height_bias_top,
                    "min-height-bias-bottom": node_compound_min_height_bias_bottom,
                },
            },
            {
                "selector": "edge",
                "style": {
                    "width": edge_line_width,
                    "curve-style": edge_curve_style,
                    "line-color": edge_line_color,
                    "line-style": edge_line_style,
                    "loop-direction": f"{edge_loop_direction}deg",
                    "loop-sweep": f"{edge_loop_sweep}deg",
                },
            },
        ]

        # Adds specified parameters if use background image is set to yes
        if use_background_image == "yes":
            if not background_image_url:
                background_image_url = "none"

            update_style(
                stylesheet=stylesheet,
                selector="node",
                addition={
                    "background-image": background_image_url,
                    "background-image-crossorigin": background_image_crossorigin,
                    "background-image-opacity": background_image_opacity,
                    "background-fit": background_image_fit,
                    "background-position-x": background_position_x,
                    "background-position-y": background_position_y,
                    "background-width-relative-to": background_width_relative_to,
                    "background-height-relative-to": background_height_relative_to,
                },
            )

        # If Background image fit is not set, we switch to using image width
        if background_image_fit == "none":
            if background_image_width is None:
                background_image_width = "auto"

            if background_image_height is None:
                background_image_height = "auto"

            update_style(
                stylesheet=stylesheet,
                selector="node",
                addition={
                    "background-width": background_image_width,
                    "background-height": background_image_height,
                },
            )

        if use_pie_chart == "yes":
            # Load json data from string format
            pie_background_color = json.loads(storage_pie_background_color)
            pie_background_size = json.loads(storage_pie_background_size)
            pie_background_opacity = json.loads(storage_pie_background_opacity)

            update_style(
                stylesheet=stylesheet,
                selector="node",
                addition={
                    "pie-size": pie_size,
                    **pie_background_color,
                    **pie_background_size,
                    **pie_background_opacity,
                },
            )

        if use_edge_arrow == "yes":
            arrow_color = json.loads(storage_arrow_color)
            arrow_shape = json.loads(storage_arrow_shape)
            arrow_fill = json.loads(storage_arrow_fill)

            update_style(
                stylesheet=stylesheet,
                selector="edge",
                addition={
                    "arrow-scale": arrow_scale,
                    **arrow_color,
                    **arrow_shape,
                    **arrow_fill,
                },
            )

        if use_edge_endpoints == "yes":
            if source_endpoint_type == "other":
                source_endpoint_width = validate_px_percentage(source_endpoint_width)
                source_endpoint_height = validate_px_percentage(source_endpoint_height)
                source_endpoint = f"{source_endpoint_width} {source_endpoint_height}"
            else:
                source_endpoint = source_endpoint_type

            if target_endpoint_type == "other":
                target_endpoint_width = validate_px_percentage(target_endpoint_width)
                target_endpoint_height = validate_px_percentage(target_endpoint_height)
                target_endpoint = f"{target_endpoint_width} {target_endpoint_height}"
            else:
                target_endpoint = target_endpoint_type

            update_style(
                stylesheet=stylesheet,
                selector="edge",
                addition={
                    "source-endpoint": source_endpoint,
                    "target-endpoint": target_endpoint,
                    "source-distance-from-node": source_distance_from_node,
                    "target-distance-from-node": target_distance_from_node,
                },
            )

        if use_labels == "yes":
            node_label_color = validate_color(node_label_color, default="black")
            edge_label_color = validate_color(edge_label_color, default="black")

            update_style(
                stylesheet=stylesheet,
                selector="node",
                addition={
                    "label": node_label,
                    # Font Styling
                    "color": node_label_color,
                    "text-opacity": node_label_text_opacity,
                    "font-family": node_label_font_family,
                    "font-size": node_label_font_size,
                    "font-style": node_label_font_style,
                    "font-weight": node_label_font_weight,
                    "text-transform": node_label_text_transform,
                    # Text Wrapping
                    "text-wrap": node_label_text_wrap,
                    "text-max-width": node_label_text_max_width,
                    # Label Alignment
                    "text-halign": label_text_halign,
                    "text-valign": label_text_valign,
                    # Text Margin
                    "text-margin-x": node_text_margin_x,
                    "text-margin-y": node_text_margin_y,
                },
            )

            update_style(
                stylesheet=stylesheet,
                selector="edge",
                addition={
                    "label": edge_label,
                    "source-label": edge_source_label,
                    "target-label": edge_target_label,
                    # Font Styling
                    "color": edge_label_color,
                    "text-opacity": edge_label_text_opacity,
                    "font-family": edge_label_font_family,
                    "font-size": edge_label_font_size,
                    "font-style": edge_label_font_style,
                    "font-weight": edge_label_font_weight,
                    "text-transform": edge_label_text_transform,
                    # Text Wrapping
                    "text-wrap": edge_label_text_wrap,
                    "text-max-width": edge_label_text_max_width,
                    # Label Alignment
                    "source-text-offset": label_source_text_offset,
                    "target-text-offset": label_target_text_offset,
                    # Text Margin
                    "text-margin-x": edge_text_margin_x,
                    "text-margin-y": edge_text_margin_y,
                    "source-text-margin-x": source_text_margin_x,
                    "source-text-margin-y": source_text_margin_y,
                    "target-text-margin-x": target_text_margin_x,
                    "target-text-margin-y": target_text_margin_y,
                },
            )

        return stylesheet
