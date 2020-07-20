import plotly.express as px
import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import plot_common
import json
from shapes_to_segmentations import (
    compute_segmentations,
    blend_image_and_classified_regions_pil,
)
from skimage import io as skio
from trainable_segmentation import multiscale_basic_features
import io
import base64
import PIL.Image
import pickle
from time import time
from joblib import Memory

memory = Memory("./joblib_cache", bytes_limit=3000000000, verbose=3)

compute_features = memory.cache(multiscale_basic_features)

DEFAULT_STROKE_WIDTH = 3  # gives line width of 2^3 = 8

DEFAULT_IMAGE_PATH = "assets/segmentation_img.jpg"

SEG_FEATURE_TYPES = ["intensity", "edges", "texture"]

# the number of different classes for labels
NUM_LABEL_CLASSES = 5
DEFAULT_LABEL_CLASS = 0
class_label_colormap = px.colors.qualitative.Light24
class_labels = list(range(NUM_LABEL_CLASSES))
# we can't have less colors than classes
assert NUM_LABEL_CLASSES <= len(class_label_colormap)

# Font and background colors associated with each theme
text_color = {"dark": "#95969A", "light": "#595959"}
card_color = {"dark": "#2D3038", "light": "#FFFFFF"}


def class_to_color(n):
    return class_label_colormap[n]


def color_to_class(c):
    return class_label_colormap.index(c)


img = skio.imread(DEFAULT_IMAGE_PATH)
features_dict = {}

app = dash.Dash(__name__)
server = app.server


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


def shapes_to_key(shapes):
    return json.dumps(shapes)


def store_shapes_seg_pair(d, key, seg, remove_old=True):
    """
    Stores shapes and segmentation pair in dict d
    seg is a PIL.Image object
    if remove_old True, deletes all the old keys and values.
    """
    bytes_to_encode = io.BytesIO()
    seg.save(bytes_to_encode, format="png")
    bytes_to_encode.seek(0)
    data = base64.b64encode(bytes_to_encode.read()).decode()
    if remove_old:
        return {key: data}
    d[key] = data
    return d


def look_up_seg(d, key):
    """ Returns a PIL.Image object """
    data = d[key]
    img_bytes = base64.b64decode(data)
    img = PIL.Image.open(io.BytesIO(img_bytes))
    return img


def generate_modal():
    with open("explanations.md", "r") as f:
        text_md = f.read()

    return html.Div(
        id="markdown",
        className="modal",
        style={"display": "none"},
        children=[
            html.Div(
                id="markdown-container",
                className="markdown-container",
                style={
                    "color": text_color["light"],
                    "backgroundColor": card_color["light"],
                },
                children=[
                    html.Div(
                        className="close-container",
                        children=html.Button(
                            "Close",
                            id="markdown_close",
                            n_clicks=0,
                            className="closeButton",
                            style={"color": "DarkBlue"},
                        ),
                    ),
                    html.Div(className="markdown-text", children=dcc.Markdown(text_md)),
                ],
            )
        ],
    )


app.layout = html.Div(
    id="app-container",
    children=[
        html.Div(
            id="banner",
            children=[
                html.H1(
                    "Interactive Machine Learning: Image Segmentation",
                    id="title",
                    className="seven columns",
                ),
                html.Img(id="logo", src=app.get_asset_url("dash-logo-new.png"),),
            ],
            className="twelve columns app-background",
        ),
        html.Div(
            id="description",
            children=[
                html.Div(
                    children=[
                        html.P(
                            'This is an example of interactive machine learning for image classification. To train the classifier, draw some marks on the picture using different colors for different parts, like in the example image. Then enable "Show segmentation" to see the classes a Random Forest Classifier gave to regions of the image, based on the marks you used as a guide. You may add more marks to clarify parts of the image where the classifier was not sucessful and the classification will update.',
                        ),
                        generate_modal(),  # modal window
                        html.Button(
                            id="learn-more-button",
                            children="Learn More",
                            n_clicks=0,
                            style={
                                "borderColor": "white",
                                "color": "#FFFFFF",
                                "backgroundColor": "DarkBlue",
                            },
                        ),
                    ],
                    className="ten columns",
                ),
                html.Img(
                    id="example-image",
                    src="assets/segmentation_img_example_marks.jpg",
                    className="two columns",
                ),
            ],
            className="twelve columns app-background",
        ),
        html.Div(
            id="main-content",
            children=[
                html.Div(
                    id="left-column",
                    children=[
                        dcc.Loading(
                            id="segmentations-loading",
                            type="circle",
                            children=[
                                # Graph
                                dcc.Graph(
                                    id="graph",
                                    figure=make_default_figure(),
                                    config={
                                        "modeBarButtonsToAdd": [
                                            "drawrect",
                                            "drawopenpath",
                                            "eraseshape",
                                        ]
                                    },
                                ),
                            ],
                        )
                    ],
                    className="six columns app-background",
                ),
                html.Div(
                    id="right-column",
                    children=[
                        html.H6("Label class"),
                        # Label class chosen with buttons
                        html.Div(
                            id="label-class-buttons",
                            children=[
                                html.Button(
                                    "%2d" % (n,),
                                    id={"type": "label-class-button", "index": n},
                                    style={"background-color": class_to_color(c)},
                                )
                                for n, c in enumerate(class_labels)
                            ],
                        ),
                        html.H6(id="stroke-width-display"),
                        # Slider for specifying stroke width
                        dcc.Slider(
                            id="stroke-width",
                            min=0,
                            max=6,
                            step=0.1,
                            value=DEFAULT_STROKE_WIDTH,
                        ),
                        # Indicate showing most recently computed segmentation
                        dcc.Checklist(
                            id="show-segmentation",
                            options=[
                                {
                                    "label": "Show segmentation",
                                    "value": "Show segmentation",
                                }
                            ],
                            value=[],
                        ),
                        html.H6("Features"),
                        dcc.Checklist(
                            id="segmentation-features",
                            options=[
                                {"label": l.capitalize(), "value": l}
                                for l in SEG_FEATURE_TYPES
                            ],
                            value=["intensity", "edges"],
                        ),
                        html.H6("Blurring parameter"),
                        dcc.RangeSlider(
                            id="sigma-range-slider",
                            min=0.01,
                            max=20,
                            step=0.01,
                            value=[0.5, 16],
                        ),
                        html.Div(
                            children=[
                                html.Button(
                                    "Download classifier", id="download-button"
                                ),
                                html.Span(
                                    "A script for using the classifier can be found in the source repository of this webapp: https://github.com/plotly/dash-sample-apps/dash-interactive-image-segmentation.",
                                    className="tooltiptext",
                                ),
                            ],
                            className="tooltip",
                        ),
                        html.A(id="download", download="classifier.json",),
                        html.Div(
                            children=[
                                html.Button(
                                    "Download classified image",
                                    id="download-image-button",
                                ),
                            ],
                        ),
                        html.A(id="download-image", download="classified-image.png",),
                    ],
                    className="six columns app-background",
                ),
            ],
            className="twelve columns",
        ),
        html.Div(
            id="no-display",
            children=[
                # Store for user created masks
                # data is a list of dicts describing shapes
                dcc.Store(id="masks", data={"shapes": []}),
                dcc.Store(id="classifier-store", data={}),
                dcc.Store(id="classified-image-store", data=""),
                dcc.Store(id="features_hash", data=""),
            ],
        ),
        html.Div(id="download-dummy"),
        html.Div(id="download-image-dummy"),
    ],
)


# Converts image classifier to a JSON compatible encoding and creates a
# dictionary that can be downloaded
# see use_ml_image_segmentation_classifier.py
def save_img_classifier(clf, label_to_colors_args, segmenter_args):
    clfbytes = io.BytesIO()
    pickle.dump(clf, clfbytes)
    clfb64 = base64.b64encode(clfbytes.getvalue()).decode()
    return {
        "classifier": clfb64,
        "segmenter_args": segmenter_args,
        "label_to_colors_args": label_to_colors_args,
    }


def show_segmentation(image_path, mask_shapes, features, segmenter_args):
    """ adds an image showing segmentations to a figure's layout """
    # add 1 because classifier takes 0 to mean no mask
    shape_layers = [color_to_class(shape["line"]["color"]) + 1 for shape in mask_shapes]
    label_to_colors_args = {
        "colormap": class_label_colormap,
        "color_class_offset": -1,
    }
    segimg, _, clf = compute_segmentations(
        mask_shapes,
        img_path=image_path,
        shape_layers=shape_layers,
        label_to_colors_args=label_to_colors_args,
        features=features,
    )
    # get the classifier that we can later store in the Store
    classifier = save_img_classifier(clf, label_to_colors_args, segmenter_args)
    segimgpng = plot_common.img_array_to_pil_image(segimg)
    return (segimgpng, classifier)


@app.callback(
    [
        Output("graph", "figure"),
        Output("masks", "data"),
        Output("stroke-width-display", "children"),
        Output("classifier-store", "data"),
        Output("classified-image-store", "data"),
    ],
    [
        Input("graph", "relayoutData"),
        Input(
            {"type": "label-class-button", "index": dash.dependencies.ALL},
            "n_clicks_timestamp",
        ),
        Input("stroke-width", "value"),
        Input("show-segmentation", "value"),
        Input("download-button", "n_clicks"),
        Input("download-image-button", "n_clicks"),
        Input("segmentation-features", "value"),
        Input("sigma-range-slider", "value"),
    ],
    [State("masks", "data"),],
)
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
        "Stroke width: %d" % (stroke_width,),
        classifier_store_data,
        classified_image_store_data,
    )


# set the download url to the contents of the classifier-store (so they can be
# downloaded from the browser's memory)
app.clientside_callback(
    """
function(the_store_data) {
    let s = JSON.stringify(the_store_data);
    let b = new Blob([s],{type: 'text/plain'});
    let url = URL.createObjectURL(b);
    return url;
}
""",
    Output("download", "href"),
    [Input("classifier-store", "data")],
)


# set the download url to the contents of the classified-image-store (so they can be
# downloaded from the browser's memory)
app.clientside_callback(
    """
function(the_image_store_data) {
    return the_image_store_data;
}
""",
    Output("download-image", "href"),
    [Input("classified-image-store", "data")],
)

# simulate a click on the <a> element when download.href is updated
app.clientside_callback(
    """
function (download_href) {
    let elem = document.querySelector('#download');
    elem.click()
    return "";
}
""",
    Output("download-dummy", "children"),
    [Input("download", "href")],
)

# simulate a click on the <a> element when download.href is updated
app.clientside_callback(
    """
function (download_image_href) {
    let elem = document.querySelector('#download-image');
    elem.click()
    return "";
}
""",
    Output("download-image-dummy", "children"),
    [Input("download-image", "href")],
)


# ======= Callback for modal popup =======
@app.callback(
    Output("markdown", "style"),
    [Input("learn-more-button", "n_clicks"), Input("markdown_close", "n_clicks")],
)
def update_click_output(button_click, close_click):
    if button_click > close_click:
        return {"display": "block"}
    else:
        return {"display": "none"}


if __name__ == "__main__":
    app.run_server(debug=True)
