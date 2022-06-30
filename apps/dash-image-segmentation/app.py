import plotly.express as px
import dash
from dash import Dash, html, dcc, Input, Output, State, callback, callback_context
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from time import time

from utils.helper_functions import (
    class_to_color,
    color_to_class,
    shapes_to_key,
    store_shapes_seg_pair,
    look_up_seg,
    save_img_classifier,
    show_segmentation,
)
from utils.figures import make_default_figure, annotation_react
from utils.components import (
    description,
    segmentation,
    sidebar,
    meta,
    header_items,
)

external_stylesheets = [dbc.themes.FLATLY, "assets/css/app.css"]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

server = app.server
app.title = "Interactive image segmentation based on machine learning"

app.layout = html.Div(
    [
        dmc.Header(height=70, padding="md", children=header_items),
        dbc.Container(
            [
                dbc.Row(description),
                dbc.Row(
                    id="app-content",
                    children=[dbc.Col(segmentation, md=8), dbc.Col(sidebar, md=4)],
                ),
                dbc.Row(dbc.Col(meta)),
            ],
            fluid=True,
        ),
    ]
)


@callback(
    Output("graph", "figure"),
    Output("masks", "data"),
    Output("stroke-width-display", "children"),
    Output("classifier-store", "data"),
    Output("classified-image-store", "data"),
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
    State("masks", "data"),
)
def return_annotation_react(
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
    return annotation_react(
        graph_relayoutData,
        any_label_class_button_value,
        stroke_width_value,
        show_segmentation_value,
        download_button_n_clicks,
        download_image_button_n_clicks,
        segmentation_features_value,
        sigma_range_slider_value,
        masks_data,
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
    Input("classified-image-store", "data"),
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
    Input("download-image", "href"),
)

# Callback for modal popup
@callback(
    Output("modal", "is_open"),
    Input("howto-open", "n_clicks"),
    Input("howto-close", "n_clicks"),
    State("modal", "is_open"),
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


# we use a callback to toggle the collapse on small screens
@callback(
    Output("navbar-collapse", "is_open"),
    Input("navbar-toggler", "n_clicks"),
    State("navbar-collapse", "is_open"),
)
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


if __name__ == "__main__":
    app.run_server()
