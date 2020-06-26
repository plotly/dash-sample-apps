import time

import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from transformers import BartTokenizer, BartForConditionalGeneration
import torch

device = 'cuda' if torch.cuda.is_available() else 'cpu'
# device = "cpu"
print(f"Device: {device}")

# Load Model
# pretrained = 'facebook/bart-large-xsum'
pretrained = "sshleifer/distilbart-xsum-12-3"
model = BartForConditionalGeneration.from_pretrained(pretrained)
tokenizer = BartTokenizer.from_pretrained(pretrained)

# Switch to cuda if available
model.to(device)
model.eval()
# Define app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

controls = dbc.Card(
    [
        dbc.FormGroup([
            dbc.Label("Output Length (# Tokens)"),
            dcc.Slider(
                id="max-length",
                min=10, max=50, value=30,
                marks={i: str(i) for i in range(10, 51, 10)},
            ),
        ]),
        dbc.FormGroup([
            dbc.Label("Beam Size"),
            dcc.Slider(
                id="num-beams",
                min=2, max=6, value=4,
                marks={i: str(i) for i in [2, 4, 6]},
            ),
        ]),
        dbc.FormGroup([
            dbc.Spinner([
                dbc.Button("Summarize", id="button-run"),
                html.Div(id="time-taken")
            ])
        ]),
    ],
    body=True, 
    style={'height': '275px'}
)


# Define Layout
app.layout = dbc.Container(fluid=True, children=[
    html.H1("Dash Automatic Summarization (with DistilBART)"),
    html.Hr(),
    dbc.Row([
        dbc.Col(width=5, children=[
            controls,
            dbc.Card(body=True, children=[
                dbc.FormGroup([
                    dbc.Label('Summarized Content'),
                    dcc.Textarea(
                        id="summarized-content",
                        style={"width": "100%", 'height': 'calc(75vh - 275px)'},
                    )
                ]),
            ])
        ]),

        dbc.Col(width=7, children=[
            dbc.Card(body=True, children=[
                dbc.FormGroup([
                    dbc.Label('Original Text (Paste here)'),
                    dcc.Textarea(
                        id="original-text", 
                        style={"width": "100%", 'height': '75vh'}
                    ),
                ])
            ])
        ]),
    ])
])


@app.callback(
    [Output("summarized-content", "value"), 
     Output("time-taken", "children")],
    [Input("button-run", "n_clicks"), 
     Input("max-length", "value"),
     Input("num-beams", "value")],
    [State("original-text", "value")],
)
def summarize(n_clicks, max_len, num_beams, original_text):
    if original_text is None or original_text == "":
        return "", "Did not run"

    t0 = time.time()

    inputs = tokenizer.batch_encode_plus([original_text], max_length=1024, return_tensors="pt")
    inputs = inputs.to(device)

    # Generate Summary
    summary_ids = model.generate(
        inputs["input_ids"], 
        num_beams=num_beams, 
        max_length=max_len, 
        early_stopping=True
    )
    out = [
        tokenizer.decode(g, skip_special_tokens=True, clean_up_tokenization_spaces=False)
        for g in summary_ids
    ]

    t1 = time.time()
    time_taken = f"Summarized on {device} in {t1-t0:.2f}s"

    return out[0], time_taken


if __name__ == "__main__":
    app.run_server(debug=True)