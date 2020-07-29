import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_core_components as dcc
import skimage.io
import skimage.transform
import io_utils


def wrap_div(e, className="control-wrapper"):
    return html.Span(children=[e], className=className)


DEFAULT_IMAGE = "assets/Euchondrus_septemdentatus_01.JPG"

app = dash.Dash(__name__)

server = app.server

app.layout = html.Div(
    id="root",
    children=[
        html.Div(
            id="main",
            children=[
                html.Div(
                    id="result-container",
                    children=[
                        html.Img(
                            id="logo",
                            src="assets/dash-logo-new.png",
                            alt="Plotly Inc. Logo",
                        ),
                        html.Div(
                            id="image-controls-container",
                            children=[
                                html.Div(
                                    id="image-container",
                                    children=[
                                        html.Img(
                                            id="image-display",
                                            alt="The transformed image",
                                        ),
                                    ],
                                ),
                                html.Div(
                                    id="controls",
                                    children=[
                                        wrap_div(
                                            dcc.Upload(
                                                id="uploader",
                                                children=html.Button(
                                                    title="Upload image",
                                                    children=html.Img(
                                                        src="assets/Upload_icon_-_Font_Awesome_-_Red.svg",
                                                        className="button-image",
                                                    ),
                                                ),
                                                multiple=False,
                                                className="inline_button",
                                            )
                                        ),
                                        wrap_div(
                                            html.A(
                                                id="downloader",
                                                download="image.png",
                                                children=[
                                                    html.Button(
                                                        title="Download image",
                                                        children=html.Img(
                                                            src="assets/Download_icon_-_Font_Awesome_-_Red.svg",
                                                            className="button-image",
                                                        ),
                                                        className="inline_button",
                                                    )
                                                ],
                                            )
                                        ),
                                        wrap_div(
                                            html.Div(
                                                "Rotation (degrees)",
                                                id="rotation-display",
                                                title="Adjust image rotation with slider",
                                            )
                                        ),
                                        dcc.Slider(
                                            id="rotation-slider",
                                            min=0,
                                            max=360,
                                            step=0.01,
                                            value=0,
                                        ),
                                    ],
                                ),
                            ],
                        ),
                    ],
                ),
                dcc.Store(
                    id="input-image", data=io_utils.mime_from_img_path(DEFAULT_IMAGE)
                ),
                html.Div(
                    id="sidebar",
                    children=[
                        html.H1(
                            id="title",
                            children=[
                                "Image processing template made with Dash and scikit-image"
                            ],
                        ),
                        html.P(
                            children=[
                                "This shows how an app for processing images can be made in ~100 lines of code with Dash and scikit-image. You can add your image processing algorithm by modifying this this app, which you can obtain here: ",
                                html.A(
                                    "https://github.com/plotly/dash-simple-image-processor",
                                    href="https://github.com/plotly/dash-simple-image-processor",
                                ),
                                ". Learn more about ",
                                html.A("Dash", href="https://plotly.com/dash/"),
                                " and ",
                                html.A(
                                    "scikit-image", href="https://scikit-image.org/"
                                ),
                                ".",
                            ]
                        ),
                        html.P(
                            children=[
                                "Credit: Upload and download icons based on ",
                                html.A(
                                    "these",
                                    href="https://commons.wikimedia.org/wiki/File:Upload_icon_-_Font_Awesome_-_Red.svg",
                                ),
                                ". Image of shells originally from ",
                                html.A(
                                    "here",
                                    href="https://commons.wikimedia.org/wiki/File:Euchondrus_septemdentatus_01.JPG",
                                ),
                                '. The font "Tex-Gyre-Heros" can be found ',
                                html.A(
                                    "here",
                                    href="https://ctan.org/pkg/tex-gyre-heros?lang=en",
                                ),
                                " and is used under the ",
                                html.A(
                                    "GUST Font License",
                                    href="https://ctan.org/license/gfl",
                                ),
                                ".",
                            ]
                        ),
                    ],
                ),
            ],
        ),
    ],
)


@app.callback(Output("input-image", "data"), [Input("uploader", "contents")])
def store_uploader_contents(uploader_contents):
    if uploader_contents is None:
        return dash.no_update
    return uploader_contents


@app.callback(
    [
        Output("image-display", "src"),
        Output("downloader", "href"),
        Output("rotation-display", "children"),
    ],
    [Input("input-image", "data"), Input("rotation-slider", "value")],
)
def process_image(input_image_data, rotation_slider_value):
    if input_image_data is None or rotation_slider_value is None:
        return dash.no_update
    im = io_utils.img_from_mime(input_image_data)

    # process the image here...
    im = skimage.transform.rotate(im, rotation_slider_value, resize=True)

    mimestr = io_utils.mime_from_img(im, format_str="jpeg")
    return (mimestr, mimestr, "Rotation: %.2f\u00B0" % rotation_slider_value)


if __name__ == "__main__":
    app.run_server()
