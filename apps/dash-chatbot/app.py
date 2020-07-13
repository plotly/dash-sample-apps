import time

import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from transformers import AutoModelWithLMHead, AutoTokenizer
import torch


# device = "cuda" if torch.cuda.is_available() else "cpu"
device = "cpu"
print(f"Device: {device}")

print("Start loading model...")
name = "microsoft/DialoGPT-medium"
tokenizer = AutoTokenizer.from_pretrained(name)
model = AutoModelWithLMHead.from_pretrained(name)

# Switch to cuda, eval mode, and FP16 for faster inference
if device == "cuda":
    model = model.half()
model.to(device)
model.eval()

print("Done.")


def textbox(text, box="other"):
    style = {
        "max-width": "55%",
        "width": "max-content",
        "padding": "10px 15px",
        "border-radius": "25px",
    }

    if box == "self":
        style["margin-left"] = "auto"
        style["margin-right"] = 0

        color = "primary"
        inverse = True

    elif box == "other":
        style["margin-left"] = 0
        style["margin-right"] = "auto"

        color = "light"
        inverse = False

    else:
        raise ValueError("Incorrect option for `box`.")

    return dbc.Card(text, style=style, body=True, color=color, inverse=inverse)


conversation = html.Div(
    style={
        "width": "80%",
        "max-width": "800px",
        "height": "70vh",
        "margin": "auto",
        "overflow-y": "auto",
    },
    id="display-conversation",
)

controls = dbc.InputGroup(
    style={"width": "80%", "max-width": "800px", "margin": "auto"},
    children=[
        dbc.Input(id="user-input", placeholder="Write to the chatbot...", type="text"),
        dbc.InputGroupAddon(dbc.Button("Submit", id="submit"), addon_type="append",),
    ],
)


# Define app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server


# Define Layout
app.layout = dbc.Container(
    fluid=True,
    children=[
        html.H1("Dash Chatbot (with DialoGPT)"),
        html.Hr(),
        dcc.Store(id="store-conversation", data=""),
        conversation,
        controls,
    ],
)


@app.callback(
    Output("display-conversation", "children"), [Input("store-conversation", "data")]
)
def update_display(chat_history):
    return [
        textbox(x, box="self") if i % 2 == 0 else textbox(x, box="other")
        for i, x in enumerate(chat_history.split(tokenizer.eos_token)[:-1])
    ]


@app.callback(
    [Output("store-conversation", "data"), Output("user-input", "value")],
    [Input("submit", "n_clicks"), Input("user-input", "n_submit")],
    [State("user-input", "value"), State("store-conversation", "data")],
)
def run_chatbot(n_clicks, n_submit, user_input, chat_history):
    if n_clicks == 0:
        return "", ""

    if user_input is None or user_input == "":
        return chat_history, ""

    # # temporary
    # return chat_history + user_input + "<|endoftext|>" + user_input + "<|endoftext|>", ""

    # encode the new user input, add the eos_token and return a tensor in Pytorch
    bot_input_ids = tokenizer.encode(
        chat_history + user_input + tokenizer.eos_token, return_tensors="pt"
    ).to(device)

    # generated a response while limiting the total chat history to 1000 tokens,
    chat_history_ids = model.generate(
        bot_input_ids, max_length=1024, pad_token_id=tokenizer.eos_token_id
    )
    chat_history = tokenizer.decode(chat_history_ids[0])

    return chat_history, ""


if __name__ == "__main__":
    app.run_server(debug=True)
