import dash
import numpy as np
import dash_html_components as html
from dash_canvas import DashCanvas
from dash_canvas.utils import array_to_data_url, parse_jsonstring
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import pytesseract

#todo: add tesseract-ocr and libtesseract-dev to venv

app = dash.Dash(__name__)
server = app.server

canvas_width = 300
canvas_height = 200

app.layout = html.Div([
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
    # html.Img(id='my-image', width=300),
        html.P(id='my-output')
    ], className="five columns"),
    ])


@app.callback(Output('my-output', 'children'),
              [Input('canvas', 'json_data')])
def update_data(string):
    if string:
        mask = parse_jsonstring(string, (canvas_height, canvas_width))
        # np.savetxt('data.csv', mask) use this to save the canvas annotations as a numpy array
    else:
        raise PreventUpdate
    # array = array_to_data_url((255 * mask).astype(np.uint8)) #todo: include outputted image as well
    text = pytesseract.image_to_string(mask)
    return text


if __name__ == '__main__':
    app.run_server(debug=True)
