import base64
from io import BytesIO

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import tensorflow as tf
import tensorflow_hub as hub
from PIL import Image


def Header(name, app):
    title = html.H1(name, style={"margin-top": 5})
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 60}
    )
    link = html.A(logo, href="https://plotly.com/dash/")

    return dbc.Row([dbc.Col(title, md=8), dbc.Col(link, md=4)])


def preprocess_b64(image_enc):
    """Preprocess b64 string into TF tensor"""
    decoded = base64.b64decode(image_enc.split("base64,")[-1])
    hr_image = tf.image.decode_image(decoded)

    if hr_image.shape[-1] == 4:
        hr_image = hr_image[..., :-1]

    return tf.expand_dims(tf.cast(hr_image, tf.float32), 0)


def tf_to_b64(tensor, ext="jpeg"):
    buffer = BytesIO()

    image = tf.cast(tf.clip_by_value(tensor[0], 0, 255), tf.uint8).numpy()
    Image.fromarray(image).save(buffer, format=ext)

    encoded = base64.b64encode(buffer.getvalue()).decode("utf-8")

    return f"data:image/{ext};base64, {encoded}"


def image_card(src, header=None):
    return dbc.Card(
        [
            dbc.CardHeader(header),
            dbc.CardBody(html.Img(src=src, style={"width": "100%"})),
        ]
    )


# Load ML model
model = hub.load("https://tfhub.dev/captain-pool/esrgan-tf2/1")


app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

controls = [
    dcc.Upload(
        dbc.Card(
            "Drag and Drop or Click",
            body=True,
            style={
                "textAlign": "center",
                "borderStyle": "dashed",
                "borderColor": "black",
            },
        ),
        id="img-upload",
        multiple=False,
    )
]


app.layout = dbc.Container(
    [
        Header("Dash Image Enhancing with TensorFlow", app),
        html.Hr(),
        dbc.Row([dbc.Col(c) for c in controls]),
        html.Br(),
        dbc.Spinner(
            dbc.Row(
                [
                    dbc.Col(html.Div(id=img_id))
                    for img_id in ["original-img", "enhanced-img"]
                ]
            )
        ),
    ],
    fluid=False,
)


@app.callback(
    [Output("original-img", "children"), Output("enhanced-img", "children")],
    [Input("img-upload", "contents")],
    [State("img-upload", "filename")],
)
def enhance_image(img_str, filename):
    if img_str is None:
        return dash.no_update, dash.no_update

    # sr_str = img_str # PLACEHOLDER
    low_res = preprocess_b64(img_str)
    super_res = model(tf.cast(low_res, tf.float32))
    sr_str = tf_to_b64(super_res)

    lr = image_card(img_str, header="Original Image")
    sr = image_card(sr_str, header="Enhanced Image")

    return lr, sr


if __name__ == "__main__":
    app.run_server(debug=True)
