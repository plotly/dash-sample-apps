import dash_core_components as dcc
import dash_html_components as html


def image_upload_zone(name, multiple=False, width=100):
    return dcc.Upload(
                id=name,
                children=[
                    'Drag and Drop or ',
                    html.A('Select an Image')],
                style={'width': str(width) + '%',
                       'height': '50px',
                       'lineHeight': '50px',
                       'borderWidth': '1px',
                       'borderStyle': 'dashed',
                       'borderRadius': '5px',
                       'textAlign': 'center'
                       },
                accept='image/*',
                multiple=multiple,
                )

