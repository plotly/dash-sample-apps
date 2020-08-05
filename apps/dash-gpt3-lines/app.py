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


# Load data
df = px.data.gapminder()

# Authentication
openai.api_key = os.getenv("OPENAI_KEY")


# Define the prompt
prompt = """
Our dataframe "df" only contains the following columns: country, continent, year, life expectancy (lifeExp), population (pop), GDP per capita (gdpPercap), the ISO alpha, the ISO numerical.

**Description**: The life expectancy in Oceania countries throughout the years.

**Code**: ```px.line(df.query("continent == 'Oceania'"), x='year', y='lifeExp', color='country', log_y=False, log_x=False)```
"""

# Create
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

content_style = {"height": "475px"}

chat_input = dbc.InputGroup(
    [
        dbc.Input(
            id="input-text", placeholder="Tell GPT-3 what you want to generate...",
        ),
        dbc.InputGroupAddon(
            dbc.Button("Submit", id="button-submit", color="primary"),
            addon_type="append",
        ),
    ]
)
output_graph = [
    dbc.CardHeader("Plotly Express Graph"),
    dbc.CardBody(dbc.Spinner(dcc.Graph(id="output-graph", style={"height": "425px"}))),
]
output_code = [
    dbc.CardHeader("GPT-3 Conversation Interface"),
    dbc.CardBody(
        dbc.Spinner(dcc.Markdown("", id="conversation-interface")),
        style={"height": "725px"},
    ),
]

explanation = f"""
*GPT-3 can generate Plotly graphs from a simple description of what you want!
We only needed to load the Gapminder dataset and give the following prompt to GPT-3:*

{prompt}
"""
explanation_card = [
    dbc.CardHeader("What am I looking at?"),
    dbc.CardBody(dcc.Markdown(explanation)),
]

left_col = [dbc.Card(output_graph), html.Br(), dbc.Card(explanation_card)]

right_col = [
    dbc.Card(output_code),
    html.Br(),
    chat_input,
]

app.layout = dbc.Container(
    [
        Header("Dash GPT-3 Chart Generation", app),
        html.Hr(),
        dbc.Row([dbc.Col(left_col, md=7), dbc.Col(right_col, md=5)]),
    ],
    fluid=True,
)


@app.callback(
    [
        Output("output-graph", "figure"),
        Output("conversation-interface", "children"),
        Output("input-text", "value"),
    ],
    [Input("button-submit", "n_clicks"), Input("input-text", "n_submit")],
    [State("input-text", "value"), State("conversation-interface", "children")],
)
def generate_graph(n_clicks, n_submit, text, conversation):
    if n_clicks is None and n_submit is None:
        return (dash.no_update,) * 3

    conversation += dedent(
        f"""
    **Description**: {text}

    **Code**: """
    )

    gpt_input = (prompt + conversation).replace("```", "").replace("**", "")

    response = openai.Completion.create(
        engine="davinci",
        prompt=gpt_input,
        max_tokens=200,
        stop=["Description:", "Code:"],
        temperature=0,
    )

    output = response.choices[0].text.strip()
    print(gpt_input)

    conversation += f"```{output}```\n"

    fig = eval(output)

    return fig, conversation, ""


if __name__ == "__main__":
    app.run_server(debug=False)
