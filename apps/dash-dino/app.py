import time
from functools import partial

import dash
import dash_bootstrap_components as dbc
import dash_labs as dl
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import requests
import torch
import torch.nn as nn
from torchvision import transforms as pth_transforms
from PIL import Image
import numpy as np
from flask_caching import Cache


def download_img(url, size=(500, 500)):
    im = Image.open(requests.get(url, stream=True).raw).convert("RGB")
    im.thumbnail(size)
    return im


def compute_attentions(model, patch_size=16):
    def aux(img):
        """
        Source: https://github.com/facebookresearch/dino/blob/main/visualize_attention.py
        """
        # make the image divisible by the patch size
        w, h = (
            img.shape[1] - img.shape[1] % patch_size,
            img.shape[2] - img.shape[2] % patch_size,
        )
        img = img[:, :w, :h].unsqueeze(0)
        w_featmap = img.shape[-2] // patch_size
        h_featmap = img.shape[-1] // patch_size
        attentions = model.forward_selfattention(img)

        return attentions, w_featmap, h_featmap

    return aux


def apply_threshold(attentions, w_featmap, h_featmap, threshold, patch_size=16):
    nh = attentions.shape[1]  # number of head
    # we keep only the output patch attention
    attentions = attentions[0, :, 0, 1:].reshape(nh, -1)
    # we keep only a certain percentage of the mass
    val, idx = torch.sort(attentions)
    val /= torch.sum(val, dim=1, keepdim=True)
    cumval = torch.cumsum(val, dim=1)
    th_attn = cumval > (1 - threshold)
    idx2 = torch.argsort(idx)
    for head in range(nh):
        th_attn[head] = th_attn[head][idx2[head]]
    th_attn = th_attn.reshape(nh, w_featmap, h_featmap).float()
    # interpolate
    th_attn = (
        nn.functional.interpolate(
            th_attn.unsqueeze(0), scale_factor=patch_size, mode="nearest"
        )[0]
        .detach()
        .cpu()
        .numpy()
    )
    attentions = attentions.reshape(nh, w_featmap, h_featmap)
    attentions = nn.functional.interpolate(
        attentions.unsqueeze(0), scale_factor=patch_size, mode="nearest"
    )
    attentions = attentions[0].detach().cpu().numpy()

    return th_attn, attentions


# vars
default_url = "https://dl.fbaipublicfiles.com/dino/img.png"

# Load model
device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
print("Running on", device)
model = torch.hub.load("facebookresearch/dino:main", "dino_deits16").to(device)
transform = pth_transforms.Compose(
    [
        pth_transforms.ToTensor(),
        pth_transforms.Normalize((0.485, 0.456, 0.406), (0.229, 0.224, 0.225)),
    ]
)

# Initialize dash app and dash-labs template
app = dash.Dash(__name__, plugins=[dl.plugins.FlexibleCallbacks()])
cache = Cache(
    app.server,
    config={
        # try 'filesystem' if you don't want to setup redis
        "CACHE_TYPE": "filesystem",
        "CACHE_DIR": "flask_cache",
    },
)
tpl = dl.templates.DbcSidebar(
    title="DINO Demo (using Dash Labs)", theme=dbc.themes.DARKLY,
)

# memoize functions
predict = cache.memoize(timeout=300)(compute_attentions(model))
download_img = cache.memoize(timeout=300)(download_img)

# Define callback function
@app.callback(
    args=dict(
        url=tpl.textbox_input(default_url, label="Image URL", kind=dl.State),
        run=tpl.button_input("Run", label=""),
        head=tpl.dropdown_input(list(range(6)), value="0", label="Attention Head"),
        options=tpl.checklist_input(["use threshold", "overlay"], []),
        threshold=tpl.slider_input(0, 1, 0.6, 0.01),
    ),
    output=tpl.graph_output(),
    template=tpl,
)
def callback(url, run, threshold, head, options):
    try:
        im = download_img(url)
    except:
        return go.Figure().update_layout(title="Incorrect URL")

    ix = int(head)

    # Run model
    img = transform(im).to(device)
    attentions, w_featmap, h_featmap = predict(img)
    th_attn, scalar_attn = apply_threshold(attentions, w_featmap, h_featmap, threshold)

    if "use threshold" in options:
        attns = th_attn
    else:
        attns = scalar_attn

    if "overlay" in options:
        fig = px.imshow(im)
        fig.add_trace(go.Heatmap(z=attns[ix], opacity=0.55))
    else:
        fig = make_subplots(1, 2)
        fig.add_trace(go.Image(z=im), 1, 1)
        fig.add_trace(go.Heatmap(z=attns[ix], y=np.arange(attns.shape[1], 0, -1)), 1, 2)
        fig.update_xaxes(matches="x")

    return fig


app.layout = tpl.layout(app)


if __name__ == "__main__":
    app.run_server(debug=True)
