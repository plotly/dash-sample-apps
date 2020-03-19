import base64
from io import BytesIO

import numpy as np
from keras.models import load_model
from PIL import Image
import dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State
import pickle
import plotly.express as px

from helpers import load_mnist, parse_image, numpy_to_b64, create_img, label_mapping

train_images, train_labels = load_mnist("./fashion", subset="train")
test_images, test_labels = load_mnist("./fashion", subset="test")

all_labels = np.concatenate((train_labels, test_labels))

X_train = train_images.reshape(60000, 28, 28, 1)
X_test = test_images.reshape(10000, 28, 28, 1)

all_images = np.concatenate((X_train, X_test), axis=0)

train_X_hat = np.load("trained_data/train_tsne.npy")
test_X_hat = np.load("trained_data/test_tsne.npy")
all_X_hat = np.load("trained_data/all_images_tsne.npy")

model = load_model("trained_data/fashion_mnist_cnn.h5")

intro_text = """
This app applies T-SNE on the images from the [fashion mnist
dataset](https://github.com/zalandoresearch/fashion-mnist), which reduces each
image to a two dimensional embedding to visualize the similarity between
images. Clusters represent similar images, and the greater the distance between
two points, the less similar the images are. The app also allows to predict the class of each image using a simple
convolutional neural network (CNN). The T-SNE embeddings were generated using [RAPIDS
cuML](https://github.com/rapidsai/cuml), and the CNN was trained using Keras.

Hover over each point in the tsne graph to see the image it represents. You can
click an individual point to see the CNN's prediction for that point, as well
as the ground-truth label. You can also upload your own image to see how the CNN would classify it, as well as
to display the images' approximate location in the T-SNE space (new embeddings
are approximated using a linear model fit on the mnist images, using the
original mnist T-SNE embeddings as the dependent variable).
"""


def create_tsne_graph(data, uploaded_point=None):
    colors = px.colors.qualitative.Pastel
    traces = []
    for i, key in enumerate(label_mapping.keys()):
        # Training data
        idx = np.where(train_labels == key)
        x = all_X_hat[idx, 0].flatten()
        y = all_X_hat[idx, 1].flatten()
        if data in ["Train", "All"]:
            opacity = 0.9
            hoverinfo = "all"
            showlegend = True
            visible = True
        else:
            opacity = 0.5
            hoverinfo = "none"
            showlegend = False
            visible = "legendonly"
        trace = {
            "x": x,
            "y": y,
            "mode": "markers",
            "type": "scattergl",
            "marker": {"color": colors[i], "size": 3},
            "name": label_mapping[key],
            "text": label_mapping[key],
            "customdata": idx[0],
            "opacity": opacity,
            "hoverinfo": hoverinfo,
            "visible": visible,
            "showlegend": showlegend,
            "selected": {"marker": {"size": 10, "color": "black"}},
        }
        traces.append(trace)

    for i, key in enumerate(label_mapping.keys()):
        # Test data
        idx = np.where(test_labels == key)
        x = all_X_hat[(idx[0] + len(train_labels)), 0].flatten()
        y = all_X_hat[(idx[0] + len(train_labels)), 1].flatten()
        if data in ["Test", "All"]:
            opacity = 0.9
            hoverinfo = "all"
            showlegend = True if data == "Test" else False
            visible = True
        else:
            opacity = 0.5
            hoverinfo = "none"
            showlegend = False
            visible = "legendonly"
        trace = {
            "x": x,
            "y": y,
            "mode": "markers",
            "type": "scattergl",
            "marker": {"color": colors[i], "size": 3},
            "name": label_mapping[key],
            "text": label_mapping[key],
            "customdata": idx[0] + len(train_labels),
            "opacity": opacity,
            "hoverinfo": hoverinfo,
            "visible": visible,
            "showlegend": showlegend,
            "selected": {"marker": {"size": 10, "color": "black"}},
        }
        traces.append(trace)

    annotation = []

    if uploaded_point:
        annotation.append(
            {
                "x": uploaded_point[0][0],
                "y": uploaded_point[0][1],
                "xref": "x",
                "yref": "y",
                "text": "Predicted Embedding for Uploaded Image",
                "showarrow": True,
                "arrowhead": 1,
                "ax": 10,
                "ay": -40,
                "font": {"size": 20},
            }
        )

    layout = {
        "xaxis": {"visible": False},
        "yaxis": {"visible": False},
        "clickmode": "event+select",
        "annotations": annotation,
    }
    figure = {"data": traces, "layout": layout}
    return figure


app = dash.Dash(name=__name__)

server = app.server

app.css.config.serve_locally = False
app.config.suppress_callback_exceptions = True

header = html.Div(
    id="app-header",
    children=[
        html.Img(src=app.get_asset_url("dash-logo.png"), className="logo"),
        "Fashion MNIST Explorer: T-SNE and CNN",
    ],
)

app.layout = html.Div(
    children=[
        header,
        html.Br(),
        html.Details(
            id="intro-text",
            children=[html.Summary(html.B("About This App")), dcc.Markdown(intro_text)],
        ),
        # html.Div(html.Div(id="intro-text", children=dcc.Markdown(intro_text),),),
        html.Div(
            id="app-body",
            children=[
                html.Div(
                    id="control-card",
                    children=[
                        html.Span(
                            className="control-label",
                            children="Display Train or Test Data",
                        ),
                        dcc.Dropdown(
                            id="train-test-dropdown",
                            className="control-dropdown",
                            options=[
                                {"label": i, "value": i}
                                for i in ["Train", "Test", "All"]
                            ],
                            value="Train",
                        ),
                        html.Span(
                            className="control-label", children="Upload an Image"
                        ),
                        dcc.Upload(
                            id="img-upload",
                            className="upload-component",
                            children=html.Div(
                                ["Drag and Drop or ", html.A("Select Files")]
                            ),
                        ),
                        html.Div(id="output-img-upload"),
                    ],
                ),
                html.Div(
                    style={"width": "75vw"},
                    children=[
                        html.Div(
                            id="tsne-graph-div",
                            children=[
                                html.Div(
                                    id="tsne-graph-outer",
                                    children=[
                                        # html.Div(
                                        # id="intro-text",
                                        # children=dcc.Markdown(intro_text),
                                        # ),
                                        html.H3(
                                            className="graph-title",
                                            children="Fashion MNIST Images Reduced to 2 Dimensions with T-SNE",
                                        ),
                                        dcc.Graph(
                                            id="tsne-graph",
                                            figure=create_tsne_graph("Test"),
                                        ),
                                    ],
                                )
                            ],
                        ),
                        html.Div(
                            id="image-card-div",
                            children=[
                                html.Div(
                                    id="hover-point-outer",
                                    className="img-card",
                                    children=[
                                        html.Div(
                                            "Hover Point:", style={"height": "20%"}
                                        ),
                                        html.Br(),
                                        html.Br(),
                                        html.Br(),
                                        html.Img(
                                            id="hover-point-graph", className="image"
                                        ),
                                    ],
                                ),
                                html.Div(
                                    id="prediction-div",
                                    className="img-card",
                                    children=[
                                        html.Div(
                                            id="selected-data-graph-outer",
                                            children=[
                                                html.Div(
                                                    children=[
                                                        html.Div("Selected Point:"),
                                                        html.Div(
                                                            id="prediction",
                                                            children=[
                                                                "Click on a point to display the Network's prediction",
                                                                html.Br(),
                                                                html.Br(),
                                                            ],
                                                        ),
                                                    ],
                                                    style={"height": "20%"},
                                                ),
                                                html.Br(),
                                                html.Img(
                                                    id="selected-data-graph",
                                                    className="image",
                                                    src=create_img(np.zeros((28, 28))),
                                                ),
                                            ],
                                        )
                                    ],
                                ),
                            ],
                        ),
                    ],
                ),
            ],
        ),
    ]
)


@app.callback(
    Output("output-img-upload", "children"),
    [Input("img-upload", "contents")],
    [State("img-upload", "filename"), State("img-upload", "last_modified")],
)
def display_uploaded_img(contents, fname, date):
    if contents is not None:
        original_img, resized_img = parse_image(contents, fname, date)

        img = np.expand_dims(resized_img, axis=0)
        prediction_array = model.predict(img)
        prediction = np.argmax(prediction_array)

        children = [
            "Your uploaded image: ",
            html.Img(className="image", src=original_img),
            "Image fed the model: ",
            html.Img(className="image", src=create_img(resized_img)),
            f"The model thinks this is a {label_mapping[prediction]}",
            html.Button(
                id="clear-button", children="Remove Uploaded Image", n_clicks=0
            ),
        ]
        return children


@app.callback(Output("img-upload", "contents"), [Input("clear-button", "n_clicks")])
def clear_upload(n_clicks):
    if n_clicks >= 1:
        return None
    raise dash.exceptions.PreventUpdate


@app.callback(
    Output("tsne-graph", "figure"),
    [Input("train-test-dropdown", "value"), Input("img-upload", "contents")],
    [State("img-upload", "filename"), State("img-upload", "last_modified")],
)
def display_train_test(value, contents, fname, date):
    embedding_prediction = None
    if contents is not None:
        original_img, resized_img = parse_image(contents, fname, date)
        linear_model = pickle.load(
            open("trained_data/linear_model_embeddings.sav", "rb")
        )
        embedding_prediction = linear_model.predict(resized_img.reshape(1, -1)).tolist()
    return create_tsne_graph(value, embedding_prediction)


@app.callback(Output("hover-point-graph", "src"), [Input("tsne-graph", "hoverData")])
def display_selected_point(hoverData):
    if not hoverData:
        return create_img(train_images[0])
    idx = hoverData["points"][0]["customdata"]
    return create_img(all_images[idx])


@app.callback(
    [Output("selected-data-graph", "src"), Output("prediction", "children")],
    [Input("tsne-graph", "clickData")],
)
def display_selected_point(clickData):
    if not clickData:
        raise dash.exceptions.PreventUpdate

    idx = clickData["points"][0]["customdata"]
    img = np.expand_dims(all_images[idx], axis=0)
    prediction_array = model.predict(img)
    prediction = np.argmax(prediction_array)
    probability = np.round(prediction_array[0, prediction] * 100, 2)
    ground_truth = all_labels[idx]
    correct = prediction == ground_truth
    if correct:
        color = "green"
    else:
        color = "red"
    return [
        create_img(all_images[idx]),
        [
            f"prediction: {label_mapping[prediction]} ({probability}% certainty)",
            html.Br(),
            f"actual: {label_mapping[ground_truth]}",
        ],
    ]


if __name__ == "__main__":
    app.run_server(debug=False)
