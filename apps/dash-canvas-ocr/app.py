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
                            "Write inside the canvas with your stylus and press Sign",
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
                                add_only=False,
                                lineColor="black",
                                goButtonTitle="Sign",
                            ),
                            className="canvas-outer",
                            style={"margin-top": "1em"},
                        ),
                    ],
                    className="v-card-content",
                ),
                html.Div(
                    html.Button(id="clear", children="clear"),
                    className="v-card-content-markdown-outer",
                ),
                html.Div(
                    [
                        html.B("Text Recognition Output", className="section_title"),
                        dcc.Loading(dcc.Markdown(id="text-output", children="")),
                    ],
                    className="v-card-content",
                    style={"margin-top": "1em"},
                ),
            ],
            className="app__content",
        ),
    ]
)


@app.callback(Output("canvas", "json_objects"), [Input("clear", "n_clicks")])
def clear_canvas(n):
    if n is None:
        return dash.no_update
    strings = ['{"objects":[ ]}', '{"objects":[]}']
    return strings[n % 2]


@app.callback(
    Output("text-output", "children"), [Input("canvas", "json_data")],
)
def update_data(string):
    if string:

        try:
            mask = parse_jsonstring(string, shape=(canvas_height, canvas_width))
        except:
            return "Out of Bounding Box, click clear button and try again"
        # np.savetxt('data.csv', mask) use this to save the canvas annotations as a numpy array
        # Invert True and False
        mask = (~mask.astype(bool)).astype(int)

        image_string = array_to_data_url((255 * mask).astype(np.uint8))

        # this is from canvas.utils.image_string_to_PILImage(image_string)
        img = Image.open(BytesIO(base64.b64decode(image_string[22:])))

        text = "{}".format(
            pytesseract.image_to_string(img, lang="eng", config="--psm 6")
        )
        return text
    else:
        raise PreventUpdate


if __name__ == "__main__":
    app.run_server(debug=True)
