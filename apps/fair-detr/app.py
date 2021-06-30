import base64
from io import BytesIO
import time

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
from PIL import Image
import requests

from model import detect, filter_boxes, detr, transform
from model import CLASSES, DEVICE


# Dash component wrappers
def Row(children=None, **kwargs):
    return html.Div(children, className="row", **kwargs)


def Column(children=None, width=1, **kwargs):
    nb_map = {
        1: "one",
        2: "two",
        3: "three",
        4: "four",
        5: "five",
        6: "six",
        7: "seven",
        8: "eight",
        9: "nine",
        10: "ten",
        11: "eleven",
        12: "twelve",
    }

    return html.Div(children, className=f"{nb_map[width]} columns", **kwargs)


# plotly.py helper functions
def pil_to_b64(im, enc="png"):
    io_buf = BytesIO()
    im.save(io_buf, format=enc)
    encoded = base64.b64encode(io_buf.getvalue()).decode("utf-8")
    return f"data:img/{enc};base64, " + encoded


def pil_to_fig(im, showlegend=False, title=None):
    img_width, img_height = im.size
    fig = go.Figure()
    # This trace is added to help the autoresize logic work.
    fig.add_trace(
        go.Scatter(
            x=[img_width * 0.05, img_width * 0.95],
            y=[img_height * 0.95, img_height * 0.05],
            showlegend=False,
            mode="markers",
            marker_opacity=0,
            hoverinfo="none",
            legendgroup="Image",
        )
    )

    fig.add_layout_image(
        dict(
            source=pil_to_b64(im),
            sizing="stretch",
            opacity=1,
            layer="below",
            x=0,
            y=0,
            xref="x",
            yref="y",
            sizex=img_width,
            sizey=img_height,
        )
    )

    # Adapt axes to the right width and height, lock aspect ratio
    fig.update_xaxes(
        showgrid=False, visible=False, constrain="domain", range=[0, img_width]
    )

    fig.update_yaxes(
        showgrid=False,
        visible=False,
        scaleanchor="x",
        scaleratio=1,
        range=[img_height, 0],
    )

    fig.update_layout(title=title, showlegend=showlegend)

    return fig


def add_bbox(
    fig,
    x0,
    y0,
    x1,
    y1,
    showlegend=True,
    name=None,
    color=None,
    opacity=0.5,
    group=None,
    text=None,
):
    fig.add_trace(
        go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y0, y0, y1, y1, y0],
            mode="lines",
            fill="toself",
            opacity=opacity,
            marker_color=color,
            hoveron="fills",
            name=name,
            hoverlabel_namelength=0,
            text=text,
            legendgroup=group,
            showlegend=showlegend,
        )
    )


# colors for visualization
COLORS = [
    "#fe938c",
    "#86e7b8",
    "#f9ebe0",
    "#208aae",
    "#fe4a49",
    "#291711",
    "#5f4b66",
    "#b98b82",
    "#87f5fb",
    "#63326e",
] * 50

RANDOM_URLS = open("random_urls.txt").read().split("\n")[:-1]
print("Running on:", DEVICE)

# Start Dash
app = dash.Dash(__name__)
server = app.server  # Expose the server variable for deployments

app.layout = html.Div(
    className="container",
    children=[
        Row(html.H1("Dash DETR Detection App")),
        Row(html.P("Input Image URL:")),
        Row(
            [
                Column(
                    width=8,
                    children=[
                        dcc.Input(
                            id="input-url",
                            style={"width": "100%"},
                            placeholder="Insert URL...",
                        ),
                    ],
                ),
                Column(html.Button("Run DETR", id="button-run", n_clicks=0), width=2),
                Column(
                    html.Button("Random Image", id="button-random", n_clicks=0), width=2
                ),
            ]
        ),
        Row(dcc.Graph(id="model-output", style={"height": "70vh"})),
        Row(
            [
                Column(
                    width=7,
                    children=[
                        html.P("Non-maximum suppression (IoU):"),
                        Row(
                            [
                                Column(
                                    width=3,
                                    children=dcc.Checklist(
                                        id="checklist-nms",
                                        options=[
                                            {"label": "Enabled", "value": "enabled"}
                                        ],
                                        value=[],
                                    ),
                                ),
                                Column(
                                    width=9,
                                    children=dcc.Slider(
                                        id="slider-iou",
                                        min=0,
                                        max=1,
                                        step=0.05,
                                        value=0.5,
                                        marks={0: "0", 1: "1"},
                                    ),
                                ),
                            ]
                        ),
                    ],
                ),
                Column(
                    width=5,
                    children=[
                        html.P("Confidence Threshold:"),
                        dcc.Slider(
                            id="slider-confidence",
                            min=0,
                            max=1,
                            step=0.05,
                            value=0.7,
                            marks={0: "0", 1: "1"},
                        ),
                    ],
                ),
            ]
        ),
    ],
)


@app.callback(
    [Output("button-run", "n_clicks"), Output("input-url", "value")],
    [Input("button-random", "n_clicks")],
    [State("button-run", "n_clicks")],
)
def randomize(random_n_clicks, run_n_clicks):
    return run_n_clicks + 1, RANDOM_URLS[random_n_clicks % len(RANDOM_URLS)]


@app.callback(
    [Output("model-output", "figure"), Output("slider-iou", "disabled")],
    [
        Input("button-run", "n_clicks"),
        Input("input-url", "n_submit"),
        Input("slider-iou", "value"),
        Input("slider-confidence", "value"),
        Input("checklist-nms", "value"),
    ],
    [State("input-url", "value")],
)
def run_model(n_clicks, n_submit, iou, confidence, checklist, url):
    apply_nms = "enabled" in checklist
    try:
        im = Image.open(requests.get(url, stream=True).raw)
    except:
        return go.Figure().update_layout(title="Incorrect URL")

    tstart = time.time()

    scores, boxes = detect(im, detr, transform, device=DEVICE)
    scores, boxes = filter_boxes(
        scores, boxes, confidence=confidence, iou=iou, apply_nms=apply_nms
    )

    scores = scores.data.numpy()
    boxes = boxes.data.numpy()

    tend = time.time()

    fig = pil_to_fig(
        im, showlegend=True, title=f"DETR Predictions ({tend-tstart:.2f}s)"
    )
    existing_classes = set()

    for i in range(boxes.shape[0]):
        class_id = scores[i].argmax()
        label = CLASSES[class_id]
        confidence = scores[i].max()
        x0, y0, x1, y1 = boxes[i]

        # only display legend when it's not in the existing classes
        showlegend = label not in existing_classes
        text = f"class={label}<br>confidence={confidence:.3f}"

        add_bbox(
            fig,
            x0,
            y0,
            x1,
            y1,
            opacity=0.7,
            group=label,
            name=label,
            color=COLORS[class_id],
            showlegend=showlegend,
            text=text,
        )

        existing_classes.add(label)

    return fig, not apply_nms


if __name__ == "__main__":
    app.run_server(debug=True)
