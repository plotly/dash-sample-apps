import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from PIL import Image
import plotly.express as px
import requests
import numpy as np
import torch
import torchvision.transforms as T
from timm.data.constants import IMAGENET_DEFAULT_MEAN, IMAGENET_DEFAULT_STD


def format_label(label):
    label = " ".join(label.split()[1:])
    label = ",".join(label.split(",")[:3])
    return label


def Header(name, app):
    title = html.H3(name, style={"margin-top": 15})
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 60}
    )
    link = html.A(logo, href="https://plotly.com/dash/")
    return dbc.Row([dbc.Col(title, md=8), dbc.Col(link, md=4)])


# Load URLs and classes
RANDOM_URLS = open("urls/random_imagenet.txt").read().split("\n")[:-1]
CLASSES = np.array(
    [format_label(x) for x in open("imagenet_labels.txt").read().split("\n")]
)

# Load image transforms
transform = T.Compose(
    [
        T.Resize(224, interpolation=3),
        T.CenterCrop(224),
        T.ToTensor(),
        T.Normalize(IMAGENET_DEFAULT_MEAN, IMAGENET_DEFAULT_STD),
    ]
)

# Load model and send it to CPU or GPU
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = torch.hub.load(
    "facebookresearch/deit:main", "deit_base_patch16_224", pretrained=True
)
model.eval().to(DEVICE)


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

app.layout = dbc.Container(
    [
        Header("Classification with Data-efficient image Transformers (DeiT)", app),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    md=7,
                    children=[
                        dcc.Graph(id="stats-display"),
                        html.P("Input Image URL:"),
                        dbc.Input(
                            id="input-url",
                            placeholder='Insert URL and click on "Run"...',
                        ),
                    ],
                ),
                dbc.Col(
                    md=5,
                    children=[
                        dcc.Graph(id="image-display"),
                        html.P("Controls:"),
                        dbc.ButtonGroup(
                            [
                                dbc.Button(
                                    "Run", id="btn-run", color="primary", n_clicks=0
                                ),
                                dbc.Button(
                                    "Random",
                                    id="btn-random",
                                    color="primary",
                                    outline=True,
                                    n_clicks=0,
                                ),
                            ],
                            style={"width": "100%"},
                        ),
                    ],
                ),
            ]
        ),
    ],
    fluid=True,
)


@app.callback(
    [Output("btn-run", "n_clicks"), Output("input-url", "value")],
    [Input("btn-random", "n_clicks")],
    [State("btn-run", "n_clicks")],
)
def randomize(n_random, n_run):
    return n_run + 1, RANDOM_URLS[n_random % len(RANDOM_URLS)]


@app.callback(
    [Output("image-display", "figure"), Output("stats-display", "figure")],
    [Input("btn-run", "n_clicks"), Input("input-url", "n_submit")],
    [State("input-url", "value")],
)
def run_model(n_clicks, n_submit, url):
    try:
        im = Image.open(requests.get(url, stream=True).raw)
    except Exception as e:
        print(e)
        return px.scatter(title="Error: " + e)

    fig = px.imshow(im, title="Original Image")

    im_pt = transform(im).unsqueeze(0)
    with torch.no_grad():
        preds = torch.softmax(model(im_pt), dim=1)
        scores = preds.numpy().squeeze()

    topk_idx = scores.argsort()[::-1][:10]
    top_classes = CLASSES[topk_idx][::-1]
    top_scores = scores[topk_idx][::-1]

    scores_fig = px.bar(
        x=top_scores,
        y=top_classes,
        labels=dict(x="Confidence", y="Classes"),
        title="ImageNet predictions by DEiT-base",
        orientation="h",
    )

    return fig, scores_fig


if __name__ == "__main__":
    app.run_server(debug=True)
