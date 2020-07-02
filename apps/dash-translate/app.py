import time

import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import torch
from transformers import MarianMTModel, MarianTokenizer
import nltk

nltk.download("punkt")
from nltk.tokenize import sent_tokenize

# Choose device and load model
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Device: {device}")
model_name = "Helsinki-NLP/opus-mt-en-ROMANCE"
tokenizer = MarianTokenizer.from_pretrained(model_name)
model = MarianMTModel.from_pretrained(model_name).to(device)

app = dash.Dash(__name__)
server = app.server

# Create app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server


# Define Layout
app.layout = dbc.Container(
    fluid=True,
    children=[
        html.H1("Dash Neural Machine Translation"),
        html.Hr(),
        dbc.Spinner(
            dbc.Row(
                [
                    dbc.Col(dbc.Button("Translate", id="button-translate"), width=2),
                    dbc.Col(
                        html.Div(id="time-output", style={"margin-top": "8px"}),
                        width=10,
                    ),
                ],
                style={"margin-bottom": "15px"},
            )
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.InputGroup(
                            [
                                dbc.InputGroupAddon(
                                    "Source Language", addon_type="prepend"
                                ),
                                dbc.Select(
                                    id="source-language",
                                    options=[{"label": "English", "value": "en"}],
                                    value="en",
                                ),
                            ]
                        ),
                        dbc.Textarea(
                            id="source-text",
                            style={"margin-top": "15px", "height": "65vh"},
                        ),
                    ]
                ),
                dbc.Col(
                    [
                        dbc.InputGroup(
                            [
                                dbc.InputGroupAddon(
                                    "Target Language", addon_type="prepend"
                                ),
                                dbc.Select(
                                    id="target-language",
                                    options=[
                                        {
                                            "label": v.replace(">>", "").replace(
                                                "<<", ""
                                            ),
                                            "value": v,
                                        }
                                        for v in tokenizer.supported_language_codes
                                    ],
                                    value=">>fr<<",
                                ),
                            ]
                        ),
                        dbc.Textarea(
                            id="target-text",
                            style={"margin-top": "15px", "height": "65vh"},
                        ),
                    ]
                ),
            ]
        ),
    ],
)


@app.callback(
    [Output("target-text", "value"), Output("time-output", "children")],
    [
        Input("button-translate", "n_clicks"),
        Input("source-language", "value"),
        Input("target-language", "value"),
    ],
    [State("source-text", "value")],
)
def translate(n_clicks, src_lang, tgt_lang, src_text):
    if src_text is None or src_text == "":
        return "", "Did not run."

    t0 = time.time()

    # Create the batched input
    text_in = [
        tgt_lang + " " + sent.replace("\n", " ") for sent in sent_tokenize(src_text)
    ]
    batch = tokenizer.prepare_translation_batch(text_in)
    for k in batch:
        batch[k] = batch[k].to(device)

    # Run the model and decode the output
    translated = model.generate(**batch)
    tgt_text = [tokenizer.decode(t, skip_special_tokens=True) for t in translated]

    t1 = time.time()
    time_output = f"Translated on {device} in {t1-t0:.2f}s"

    return " ".join(tgt_text), time_output


if __name__ == "__main__":
    app.run_server(debug=True)
