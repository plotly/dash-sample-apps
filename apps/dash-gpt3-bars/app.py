import black
import os
from textwrap import dedent

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dash import no_update
from dash.dependencies import Input, Output, State
import plotly.express as px
import openai


def Header(name, app):
    title = html.H1(name, style={"margin-top": 5})
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 60}
    )
    return dbc.Row([dbc.Col(title, md=8), dbc.Col(logo, md=4)])


# Authentication
openai.api_key = os.getenv("OPENAI_KEY")


# Define the prompt
desc = "Our zoo has three twenty giraffes, fourteen orangutans, 3 monkeys more than the number of giraffes we have."
code_exp = "px.bar(x=['giraffes', 'orangutans', 'monkeys'], y=[20, 14, 23], labels=dict(x='animals', y='count'), title='Our Zoo')"
formatted_exp = black.format_str(code_exp, mode=black.FileMode(line_length=50))

# Create
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

content_style = {"height": "475px"}

controls = [
    dbc.InputGroup(
        [
            dbc.Input(
                id="input-text",
                placeholder="Specify what you want GPT-3 to generate...",
            ),
            dbc.InputGroupAddon(
                dbc.Button("Submit", id="button-submit", color="primary"),
                addon_type="append",
            ),
        ]
    )
]
output_graph = [
    dbc.CardHeader("Plotly Express Graph"),
    dbc.CardBody(dcc.Graph(id="output-graph", style={"height": "400px"})),
]
output_code = [
    dbc.CardHeader("GPT-3 Generated Code"),
    dcc.Markdown(id="output-code", style={"margin": "50px 5px"}),
]

explanation = f"""
*GPT-3 can generate Plotly graphs from a simple description of what you want!
We only needed to give the following prompt to GPT-3:*

Description: **{desc}**

Code:
```
{code_exp}
```
"""
explanation_card = [
    dbc.CardHeader("What am I looking at?"),
    dbc.CardBody(dcc.Markdown(explanation)),
]

app.layout = dbc.Container(
    [
        Header("Dash GPT-3 Chart Generation", app),
        html.Hr(),
        html.Div(controls, style={"padding-bottom": "15px"}),
        dbc.Spinner(
            dbc.Row(
                [
                    dbc.Col(dbc.Card(comp, style=content_style))
                    for comp in [output_graph, output_code]
                ],
                style={"padding-bottom": "15px"},
            )
        ),
        dbc.Card(explanation_card),
    ],
    fluid=False,
)


@app.callback(
    [Output("output-graph", "figure"), Output("output-code", "children")],
    [Input("button-submit", "n_clicks"), Input("input-text", "n_submit")],
    [State("input-text", "value")],
)
def generate_graph(n_clicks, n_submit, text):
    if text is None:
        return dash.no_update, dash.no_update

    prompt = dedent(
        f"""
        description: {desc}
        code:
        {code_exp}

        description: {text}
        code:
        """
    ).strip("\n")

    response = openai.Completion.create(
        engine="davinci",
        prompt=prompt,
        max_tokens=100,
        stop=["description:", "code:"],
        temperature=0,
    )
    output = response.choices[0].text.strip()

    code = f"import plotly.express as px\nfig = {output}\nfig.show()"
    formatted = black.format_str(code, mode=black.FileMode(line_length=50))

    try:
        fig = eval(output).update_layout(margin=dict(l=35, r=35, t=35, b=35))
    except Exception as e:
        fig = px.line(title=f"Exception: {e}. Please try again!")

    return fig, f"```\n{formatted}\n```"


if __name__ == "__main__":
    app.run_server(debug=False)
