from PIL import Image
import base64
from io import BytesIO

import dash
import numpy as np
import dash_html_components as html
from dash_canvas import DashCanvas
from dash_canvas.utils import array_to_data_url, parse_jsonstring
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

import pytesseract

app = dash.Dash(__name__)
server = app.server

canvas_width = 600
canvas_height = 450

app.layout = html.Div([
    # Banner
    html.Div(
        [
            html.Img(src=app.get_asset_url("logo.png"), className="app__logo"),
            html.H4("Dash OCR", className="header__text"),
        ],
        className="app__header",
    ),
    html.H6('Draw on image and press Save to show annotations geometry'),
    html.Div([
        DashCanvas(id='canvas',
                   lineWidth=2,
                   width=canvas_width,
                   hide_buttons=["zoom", "pan", "line", "pencil", "rectangle", "undo", "select"],
                   lineColor='black',
                   goButtonTitle='Sign'
                   ),
    ], className="five columns"),
    html.Div([
        html.Img(id='my-image', width=300),
        html.Div(id='text-output', children=''),
        html.P(id='my-output')
    ], className="five columns"),
])


@app.callback(Output('text-output', 'children'),
              [Input('canvas', 'json_data')])
def update_data(string):
    if string:
        mask = parse_jsonstring(string, shape=(canvas_height, canvas_width))
        # np.savetxt('data.csv', mask) use this to save the canvas annotations as a numpy array

        image_string = array_to_data_url((255 * mask).astype(np.uint8))  # todo: include outputted image as well

        # this is from canvas.utils.image_string_to_PILImage(image_string)
        img = Image.open(BytesIO(base64.b64decode(image_string[22:]))) #try save img to see what it looks like?
        print('img', img)
        text = pytesseract.image_to_string(img, lang='eng', config='--psm 7')
        print('text', text)
        return text #todo : handle condition which ocr cannot recognize: return message: "enpty, try again"

    else:
        raise PreventUpdate


if __name__ == '__main__':
    app.run_server(debug=True)
