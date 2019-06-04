from glob import glob
import numpy as np
import pandas as pd
from skimage import io, data, transform
from time import sleep


import dash
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import dash_table

import dash_canvas
from dash_canvas.components import image_upload_zone
from dash_canvas.utils import (image_string_to_PILImage, array_to_data_url,
                               parse_jsonstring_line,
                               brightness_adjust, contrast_adjust)
from registration import register_tiles
from utils import StaticUrlPath

def tile_images(list_of_images, n_rows, n_cols):
    dtype = list_of_images[0].dtype
    if len(list_of_images) < n_rows * n_cols:
        white = np.zeros(list_of_images[0].shape, dtype=dtype)
        n_missing = n_rows * n_cols - len(list_of_images)
        list_of_images += [white, ] * n_missing
    return np.vstack([np.hstack(list_of_images[i_row*n_cols:
                                               i_row*n_cols + n_cols])
                      for i_row in range(n_rows)])


def untile_images(image_string, n_rows, n_cols):
    big_im = np.asarray(image_string_to_PILImage(image_string))
    tiles = [np.split(im, n_cols, axis=1) for im in np.split(big_im, n_rows)]
    return np.array(tiles)


def demo_data():
    im = data.immunohistochemistry()
    l_c = 128
    l_r = 180
    n_rows = 1
    n_cols = im.shape[1] // l_c
    init_i, init_j = 0, 0
    overlap_h = [5, 25]

    big_im = np.empty((n_rows * l_r, n_cols * l_c, 3), dtype=im.dtype)
    i = 0
    for j in range(n_cols):
        sub_im = im[init_i:init_i + l_r, init_j:init_j + l_c]
        big_im[i*l_r : (i+1)*l_r, j*l_c:(j+1)*l_c] = sub_im
        init_j += l_c - overlap_h[1]
        init_i += overlap_h[0]
    return big_im


def _sort_props_lines(props, height, width, ncols):
    props = pd.DataFrame(props)
    index_init = ncols * (((props['top']) - (props['height']) //2) // height) + \
                    (((props['left']) - (props['width']) //2)// width)
    index_end = ncols * (((props['top']) + (props['height']) //2) // height) + \
                    (((props['left']) + (props['width']) //2)// width)
    props['index_init'] = index_init
    props['index_end'] = index_end
    overlaps = {}
    for line in props.iterrows():
        overlaps[(line[1]['index_init'], line[1]['index_end'])] = (int(line[1]['height']),
                                                            int(line[1]['width']))
    return overlaps


def instructions():
    return html.Div(children=[
    dcc.Markdown("""
    - Choose the number of rows and columns of the mosaic
    - Upload images
    - Try automatic stitching by pressing
    the "Run stitching" button
    - If automatic stitching did not work,
    try adjusting the overlap parameter
    """)
    ])


height, width = 200, 500
canvas_width = 800
canvas_height = round(height * canvas_width / width)
scale = canvas_width / width

list_columns = ['length', 'width', 'height', 'left', 'top']
columns = [{"name": i, "id": i} for i in list_columns]

app = dash.Dash(__name__)
server = app.server
app.config.suppress_callback_exceptions = True


app.layout = html.Div([

    html.Div([
        html.H1(
            children='Stitching App'
        ), 
        instructions(),
        html.Button('LEARN MORE'),
        html.Label('Number of rows'),
        dcc.Input(
            id='nrows-stitch',
            type='number',
            value=1,
            name='number of rows',
        ),
        html.Label('Number of columns'),
        dcc.Input(
            id='ncolumns-stitch',
            type='number',
            value=4,
            name='number of columns',
        ),
    
        html.Label('Fraction of overlap (in [0-1] range)'),
        dcc.Input(
            id='overlap-stitch',
            type='number',
            value=0.15,
            min=0,
            max=1
        ),
        html.Label('Measured shifts between images'),
        dash_table.DataTable(
            id='table-stitch',
            columns=columns,
            editable=True,
        ),
        html.Br(),
        html.Button('Run stitching', id='button-stitch',
                    style={'color': 'red'}),
        html.Br()
    ], className="three columns"),
    html.Div([
        dcc.Tabs(
            id='stitching-tabs',
            value='canvas-tab',
            children=[
                dcc.Tab(
                    label='Image tiles',
                    value='canvas-tab',
                ),
                dcc.Tab(
                    label='Stitched Image',
                    value='result-tab',
                    ),
                dcc.Tab(
                    label='How to use this app',
                    value='help-tab',
                        )
            ]
            ),
        html.Div(id='tabs-content-example'),
        html.Button('Upload demo data', id='demo'),
        image_upload_zone('upload-stitch', multiple=True,
                          width=45),
	html.Div(id='sh_x', hidden=True),
        html.Div(id='stitched-res', hidden=True),
        dcc.Store(id='memory-stitch'),
    ], className="eight columns"),
   
    ])


@app.callback(Output('tabs-content-example', 'children'),
              [Input('stitching-tabs', 'value')])
def fill_tab(tab):
    if tab=='canvas-tab':
        return [dash_canvas.DashCanvas(
                    id='canvas-stitch',
                    width=canvas_width,
                    height=canvas_height,
                    scale=scale,
                    lineWidth=2,
                    lineColor='red',
                    tool="line",
                    hide_buttons=['pencil'],
                    image_content=array_to_data_url(
                        np.zeros((height, width), dtype=np.uint8)),
                    goButtonTitle='Estimate translation',
                    ),
                ]
    elif tab=='result-tab':
        return [
                dcc.Loading(id='loading-1', children=[
                html.Img(id='stitching-result',
                    src=array_to_data_url(
                        np.zeros((height, width), dtype=np.uint8)),
                    width=canvas_width)],
                type='circle'),
                html.Div([
                html.Label('Contrast'),
                dcc.Slider(id='contrast-stitch',
                            min=0,
                            max=1,
                            step=0.02,
                            value=0.5)],
                style={'width':'40%'}),
                html.Div([
                html.Label('Brightness'),
                dcc.Slider(id='brightness-stitch',
                            min=0,
                            max=1,
                            step=0.02,
                            value=0.5,)], 
                style={'width':'40%'}),
                ]
    else: 
        return [html.Img(id='bla', src='assets/stitching.gif',
                            width=canvas_width),
                ]


@app.callback(Output('stitching-tabs', 'value'),
            [Input('button-stitch', 'n_clicks')])
def change_focus(click):
    if click:
        return 'result-tab'
    return 'canvas-tab'


@app.callback(Output('table-stitch', 'data'),
            [Input('canvas-stitch', 'json_data')])
def estimate_translation(string):
    props = parse_jsonstring_line(string)
    if props is not None and len(props) > 0:
        df = pd.DataFrame(props, columns=list_columns)
        return df.to_dict("records")
    else:
        raise PreventUpdate




@app.callback(Output('memory-stitch', 'data'),
            [Input('button-stitch', 'n_clicks')])
def update_store(click):
    sleep(.1)
    return click



@app.callback(Output('stitching-result', 'src'),
            [Input('contrast-stitch', 'value'),
             Input('brightness-stitch', 'value'),
             Input('stitched-res', 'children')])
def modify_result(contrast, brightness, image_string):
    img = np.asarray(image_string_to_PILImage(image_string))
    img = contrast_adjust(img, contrast)
    img = brightness_adjust(img, brightness)
    return array_to_data_url(img)


@app.callback(Output('stitched-res', 'children'),
            [Input('button-stitch', 'n_clicks')],
            [State('nrows-stitch', 'value'),
            State('ncolumns-stitch', 'value'),
            State('overlap-stitch', 'value'),
            State('table-stitch', 'data'),
            State('sh_x', 'children')])
def modify_content(n_cl,
                   n_rows, n_cols, overlap, estimate, image_string, vals):
    blending = 1 in vals
    if image_string is None:
        raise PreventUpdate
    tiles = untile_images(image_string, n_rows, n_cols)
    if estimate is not None and len(estimate) > 0:
        overlap_dict = _sort_props_lines(estimate, tiles.shape[2],
                                            tiles.shape[3], n_cols)
    else:
        overlap_dict = None
    canvas = register_tiles(tiles, n_rows, n_cols,
                            overlap_global=overlap,
                            overlap_local=overlap_dict,
                            pad=np.max(tiles.shape[2:]),
                            blending=blending)
    return array_to_data_url(canvas)


@app.callback(Output('canvas-stitch', 'image_content'),
            [Input('sh_x', 'children')])
def update_canvas_image(im):
    return im


@app.callback(Output('upload-stitch', 'contents'),
             [Input('demo', 'n_clicks')])
def reset_contents(n_clicks):
    if n_clicks:
        return None

@app.callback(Output('sh_x', 'children'),
            [Input('upload-stitch', 'contents'),
            Input('upload-stitch', 'filename'),
            Input('demo', 'n_clicks'),
            Input('downsample', 'value'),],
            [State('nrows-stitch', 'value'),
            State('ncolumns-stitch', 'value')])
def upload_content(list_image_string, list_filenames, click, downsample,
                    n_rows, n_cols):
    if list_image_string is not None:
        order = np.argsort(list_filenames)
        image_list = [np.asarray(image_string_to_PILImage(
                            list_image_string[i])) for i in order]
        if downsample > 1:
            ratio = 1. / downsample
            multichannel = image_list[0].ndim > 2
            image_list = [transform.rescale(image, ratio,
                                    multichannel=multichannel,
                                    preserve_range=True).astype(np.uint8)
                                    for image in image_list]
        res = tile_images(image_list, n_rows, n_cols)
        return array_to_data_url(res)
    elif click:
        res = demo_data()
        tmp = array_to_data_url(res)
        return tmp
    else:
        raise PreventUpdate


if __name__ == '__main__':
    app.run_server(debug=True)

