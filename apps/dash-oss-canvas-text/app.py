from PIL import Image
import base64
from io import BytesIO

import dash
import numpy as np
import dash_html_components as html
import dash_core_components as dcc
from dash_canvas import DashCanvas
from dash_canvas.utils import array_to_data_url, parse_jsonstring
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

# from iam_model.main import main

import pytesseract

app = dash.Dash(__name__)
server = app.server

canvas_width = 800
canvas_height = 200

app.layout = html.Div(
    [
        # Banner
        html.Div(
            [
                html.Img(src=app.get_asset_url("ocr-logo.png"), className="app__logo"),
                html.H4("Dash OCR", className="header__text"),
            ],
            className="app__header",
        ),
        # Canvas
        html.Div(
            [
                html.Div(
                    [
                        html.P(
                            "Write inside the canvas with your pencil and press Sign",
                            className="section_title",
                        ),
                        html.Div(
                            DashCanvas(
                                id="canvas",
                                lineWidth=8,
                                width=canvas_width,
                                height=canvas_height,
                                hide_buttons=[
                                    "zoom",
                                    "pan",
                                    "line",
                                    "pencil",
                                    "rectangle",
                                    "select",
                                ],
                                lineColor="black",
                                goButtonTitle="Sign",
                            ),
                            className="canvas-outer",
                            style={'margin-top': '1em'}
                        ),
                    ],
                    className="v-card-content",
                ),

                html.Div(
                    dcc.Markdown(id='reset-markdown', children="[CLEAR]({})".format(app.get_relative_path('/'))),
                    className='v-card-content-markdown-outer'
                ),
                # Annotation Geometry
                html.Div(
                    [
                        html.P(
                            "Handwriting annotation geometry", className="section_title"
                        ),
                        html.Img(id="my-image", width=canvas_width),
                    ],
                    style={'margin-top': '1em'},
                    className="v-card-content",
                ),
                # OCR output div
                html.Div(
                    [
                        html.B("Text Recognition Output", className="section_title"),
                        html.P("Pytesseract output: "),
                        dcc.Loading(dcc.Markdown(id="text-output", children=""))],
                    className="v-card-content",
                    style={'margin-top': '1em'}
                ),
                # html.Div(
                #     [
                #         html.P(
                #             "Handwriting annotation geometry", className="section_title"
                #         ),
                #         html.Img(id="my-image-iam", width=canvas_width,),
                #     ],
                #     className="v-card-content",
                # ),
                # OCR output div
                # html.Div(
                #     [
                #         html.P("IAM Trained Model Output: "),
                #         dcc.Loading(dcc.Markdown(id="text-output-iam", children=""))],
                #     className="v-card-content",
                # ),
            ],
            className="app__content",
        ),
    ]
)


@app.callback(
    [Output("text-output", "children"), Output("my-image", "src")],
    [Input("canvas", "json_data")],
)
def update_data(string):
    if string:
        emtpy_array = np.ones((150, 600), dtype=np.bool)  # Temp set canvas height with empty image
        empty_img = array_to_data_url(emtpy_array)

        try:
            mask = parse_jsonstring(string, shape=(canvas_height, canvas_width))
        except:
            return "Out of Bounding Box, click clear button and try again", empty_img
        # np.savetxt('data.csv', mask) use this to save the canvas annotations as a numpy array
        # print(mask)
        # Invert True and False
        mask = (~mask.astype(bool)).astype(int)
        # print(mask)

        image_string = array_to_data_url(
            (255 * mask).astype(np.uint8)
        )  # todo: include outputted image as well

        # this is from canvas.utils.image_string_to_PILImage(image_string)
        img = Image.open(
            BytesIO(base64.b64decode(image_string[22:]))
        )  # try save img to see what it looks like?

        # img.save("geeks2.png")
        # print('img', img)
        text = "{}".format(
            pytesseract.image_to_string(img, lang="eng", config="--psm 6")
        )
        # print('text', text)
        return (
            text,
            image_string,
        )  # todo : handle condition which ocr cannot recognize: return message: "enpty, try again"

    else:
        raise PreventUpdate


# @app.callback(Output("text-output-iam", "children"), [Input("canvas", "json_data")])
# def update_iam_output(string):
#     if string:
#         mask = parse_jsonstring(string, shape=(canvas_height, canvas_width))
#         # np.savetxt('data.csv', mask) use this to save the canvas annotations as a numpy array
#         # print(mask)
#         # Invert True and False
#         mask = (~mask.astype(bool)).astype(int)
#         # print(mask)
#
#         image_string = array_to_data_url(
#             (255 * mask).astype(np.uint8)
#         )  # todo: include outputted image as well
#
#         # this is from canvas.utils.image_string_to_PILImage(image_string)
#         img = Image.open(
#             BytesIO(base64.b64decode(image_string[22:]))
#         )  # try save img to see what it looks like?
#
#         img.save("my_writing.png")
#         # print('img', img)
#         # text = main(img)
#         text = "{}".format(str(main()))
#         return text  # todo : handle condition which ocr cannot recognize: return message: "enpty, try again"
#
#     else:
#         raise PreventUpdate


if __name__ == "__main__":
    app.run_server(debug=True)
