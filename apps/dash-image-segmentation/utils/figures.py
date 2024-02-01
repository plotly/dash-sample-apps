import utils.plot_common as plot_common
from constants import (
    DEFAULT_IMAGE_PATH,
    DEFAULT_LABEL_CLASS,
    DEFAULT_STROKE_WIDTH,
)
import dash
import PIL.Image
from utils.helper_functions import class_to_color, show_segmentation
from time import time
from constants import compute_features, img, DEFAULT_LABEL_CLASS
from utils.shapes_to_segmentations import (
    compute_segmentations,
    blend_image_and_classified_regions_pil,
)

import io


def make_default_figure(
    images=[DEFAULT_IMAGE_PATH],
    stroke_color=class_to_color(DEFAULT_LABEL_CLASS),
    stroke_width=DEFAULT_STROKE_WIDTH,
    shapes=[],
):
    fig = plot_common.dummy_fig()
    plot_common.add_layout_images_to_fig(fig, images)
    fig.update_layout(
        {
            "dragmode": "drawopenpath",
            "shapes": shapes,
            "newshape.line.color": stroke_color,
            "newshape.line.width": stroke_width,
            "margin": dict(l=0, r=0, b=0, t=0, pad=4),
        }
    )
    return fig


def annotation_react(
    graph_relayoutData,
    any_label_class_button_value,
    stroke_width_value,
    show_segmentation_value,
    download_button_n_clicks,
    download_image_button_n_clicks,
    segmentation_features_value,
    sigma_range_slider_value,
    masks_data,
):
    classified_image_store_data = dash.no_update
    classifier_store_data = dash.no_update
    cbcontext = [p["prop_id"] for p in dash.callback_context.triggered][0]
    if cbcontext in ["segmentation-features.value", "sigma-range-slider.value"] or (
        ("Show segmentation" in show_segmentation_value)
        and (len(masks_data["shapes"]) > 0)
    ):
        segmentation_features_dict = {
            "intensity": False,
            "edges": False,
            "texture": False,
        }
        for feat in segmentation_features_value:
            segmentation_features_dict[feat] = True
        t1 = time()
        features = compute_features(
            img,
            **segmentation_features_dict,
            sigma_min=sigma_range_slider_value[0],
            sigma_max=sigma_range_slider_value[1],
        )
        t2 = time()
        print(t2 - t1)
    if cbcontext == "graph.relayoutData":
        if "shapes" in graph_relayoutData.keys():
            masks_data["shapes"] = graph_relayoutData["shapes"]
        else:
            return dash.no_update
    stroke_width = int(round(2 ** (stroke_width_value)))
    # find label class value by finding button with the most recent click
    if any_label_class_button_value is None:
        label_class_value = DEFAULT_LABEL_CLASS
    else:
        label_class_value = max(
            enumerate(any_label_class_button_value),
            key=lambda t: 0 if t[1] is None else t[1],
        )[0]

    fig = make_default_figure(
        stroke_color=class_to_color(label_class_value),
        stroke_width=stroke_width,
        shapes=masks_data["shapes"],
    )
    # We want the segmentation to be computed
    if ("Show segmentation" in show_segmentation_value) and (
        len(masks_data["shapes"]) > 0
    ):
        segimgpng = None
        try:
            feature_opts = dict(segmentation_features_dict=segmentation_features_dict)
            feature_opts["sigma_min"] = sigma_range_slider_value[0]
            feature_opts["sigma_max"] = sigma_range_slider_value[1]
            segimgpng, clf = show_segmentation(
                DEFAULT_IMAGE_PATH, masks_data["shapes"], features, feature_opts
            )
            if cbcontext == "download-button.n_clicks":
                classifier_store_data = clf
            if cbcontext == "download-image-button.n_clicks":
                classified_image_store_data = plot_common.pil_image_to_uri(
                    blend_image_and_classified_regions_pil(
                        PIL.Image.open(DEFAULT_IMAGE_PATH), segimgpng
                    )
                )
        except ValueError:
            # if segmentation fails, draw nothing
            pass
        images_to_draw = []
        if segimgpng is not None:
            images_to_draw = [segimgpng]
        fig = plot_common.add_layout_images_to_fig(fig, images_to_draw)
    fig.update_layout(uirevision="segmentation")
    return (
        fig,
        masks_data,
        "Current paintbrush width: %d" % (stroke_width,),
        classifier_store_data,
        classified_image_store_data,
    )
