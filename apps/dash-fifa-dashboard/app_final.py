import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import dash_table
from dash_table import DataTable
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc

# Dataset Processing

# importing data
data = pd.read_csv("archive/players_21.csv")

# data cleaning
nonusefulcolumns = ["sofifa_id", "player_url", "long_name", "league_rank"]
nonusefulattributes = data.loc[:, "player_traits":]

df = data.copy()
df = df[df["player_positions"] != "GK"]  # filtering out all the goalkeepers
df1 = df[df["age"] > 25]  # dataset for players over 25
df2 = df[df["age"] <= 25]  # dataset for players under 25

# variables for the analysis
skill_player = ["pace", "shooting", "passing", "dribbling", "defending", "physic"]
info_player = [
    "short_name",
    "nationality",
    "club_name",
    "age",
    "height_cm",
    "weight_kg",
]
labels_table = ["Player", "Nationality", "Club", "Age", "Height (cm)", "Weight (kg)"]
skills1 = [
    "skill_curve",
    "skill_dribbling",
    "skill_fk_accuracy",
    "skill_ball_control",
    "skill_long_passing",
]
player1 = "Lionel Andrés Messi Cuccittini"
player2 = "Kylian Mbappé Lottin"


###################################################   Interactive Components   #########################################
# choice of the players
players_options_over_25 = []
for i in df1.index:
    players_options_over_25.append(
        {"label": df1["long_name"][i], "value": df1["long_name"][i]}
    )

players_options_under_25 = []
for i in df2.index:
    players_options_under_25.append(
        {"label": df2["long_name"][i], "value": df2["long_name"][i]}
    )

dropdown_player_over_25 = dcc.Dropdown(
    id="player1",
    options=players_options_over_25,
    value="Lionel Andrés Messi Cuccittini",
)

dropdown_player_under_25 = dcc.Dropdown(
    id="player2", options=players_options_under_25, value="Kylian Mbappé Lottin"
)

dashtable_1 = dash_table.DataTable(
    id="table1",
    columns=[
        {"name": col, "id": info_player[idx]} for (idx, col) in enumerate(labels_table)
    ],
    data=df[df["long_name"] == player1].to_dict("records"),
    style_cell={"textAlign": "left", "font_size": "14px"},
    style_data_conditional=[
        {"if": {"row_index": "odd"}, "backgroundColor": "rgb(248, 248, 248)"}
    ],
    style_header={"backgroundColor": "rgb(230, 230, 230)", "fontWeight": "bold"},
)


dashtable_2 = dash_table.DataTable(
    id="table2",
    # columns=[{"name": i, "id": i} for i in info_player[::-1]],
    columns=[
        {"name": col, "id": info_player[::-1][idx]}
        for (idx, col) in enumerate(labels_table[::-1])
    ],
    data=df[df["long_name"] == player2].to_dict("records"),
    style_cell={"textAlign": "right", "font_size": "14px"},
    style_data_conditional=[
        {"if": {"row_index": "odd"}, "backgroundColor": "rgb(248, 248, 248)"}
    ],
    style_header={"backgroundColor": "rgb(230, 230, 230)", "fontWeight": "bold"},
)


data.drop(nonusefulcolumns, axis=1, inplace=True)
data.drop(nonusefulattributes, axis=1, inplace=True)
data["isOver25"] = data["age"] > 25
data["isOver25text"] = data["isOver25"].map(lambda x: "Over 25" if x else "Under 25")
################Components##############################################

options = [
    {"label": "Overall", "value": "overall"},
    {"label": "Potential", "value": "potential"},
    {"label": "Value", "value": "value_eur"},
    {"label": "Wage", "value": "wage_eur"},
    {"label": "Height", "value": "height_cm"},
    {"label": "Weight", "value": "weight_kg"},
    {"label": "Pace", "value": "pace"},
    {"label": "Shooting", "value": "shooting"},
    {"label": "Passing", "value": "passing"},
    {"label": "Dribbling", "value": "dribbling"},
    {"label": "Defending", "value": "defending"},
    {"label": "Physic", "value": "physic"},
]
top_10_leagues = [
    "Spain Primera Division",
    "Italian Serie A",
    "Spain Primera Division",
    "German 1. Bundesliga",
    "French Ligue 1",
    "English Premier League",
    "Portuguese Liga ZON SAGRES",
    "Belgian Jupiler Pro League",
    "Holland Eredivisie",
    "Russian Premier League",
]


leagues = [
    {"label": "Premier League", "value": "English Premier League"},
    {"label": "Bundesliga", "value": "German 1. Bundesliga"},
    {"label": "Serie A", "value": "Italian Serie A"},
    {"label": "Liga ZON SAGRES", "value": "Portuguese Liga ZON SAGRES"},
    {"label": "La Liga", "value": "Spain Primera Division"},
    {"label": "Ligue 1", "value": "French Ligue 1"},
    {"label": "Jupiler Pro League", "value": "Belgian Jupiler Pro League"},
    {"label": "Eredivisie", "value": "Holland Eredivisie"},
    {"label": "Russian Premier League", "value": "Russian Premier League"},
]

value_x = [
    {"label": "Overall", "value": "overall"},
    {"label": "Potential", "value": "potential"},
    {"label": "Value", "value": "value_eur"},
    {"label": "Wage", "value": "wage_eur"},
    {"label": "Height", "value": "height_cm"},
    {"label": "Weight", "value": "weight_kg"},
    {"label": "Pace", "value": "pace"},
    {"label": "Shooting", "value": "shooting"},
    {"label": "Passing", "value": "passing"},
    {"label": "Dribbling", "value": "dribbling"},
    {"label": "Defending", "value": "defending"},
    {"label": "Physic", "value": "physic"},
]


metric1_dropdown = dcc.Dropdown(id="drop1", options=options, value="overall")

metric2_dropdown = dcc.Dropdown(id="drop2", options=options, value="potential")

metric3_dropdown = dcc.Dropdown(id="drop3", options=options, value="value_eur")

metric_club_dropdown1 = dcc.Dropdown(
    id="club-drop1", options=leagues, value="English Premier League"
)
metric_scatter_dropdown1 = dcc.Dropdown(
    id="scatter-drop1", options=value_x, value="value_eur"
)
metric_scatter_dropdown2 = dcc.Dropdown(
    id="scatter-drop2", options=value_x, value="wage_eur"
)


age_slider = dcc.RangeSlider(
    id="age_slider",
    min=data["age"].min(),
    max=data["age"].max(),
    value=[data["age"].min(), data["age"].max()],
    step=1,
    marks={
        16: "16",
        20: "20",
        24: "24",
        28: "28",
        32: "32",
        36: "36",
        40: "40",
        44: "44",
        48: "48",
        52: "52",
    },
)

########Dash App Layout##########################

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

server = app.server

navbar = dbc.Navbar(
    [
        html.A(
            # Use row and col to control vertical alignment of logo / brand
            dbc.Row(
                [
                    dbc.Col(
                        html.Img(
                            src=app.get_asset_url("logo_white3.png"), height="130px"
                        ),
                        width=3,
                    ),
                    dbc.Col(
                        [
                            html.Label("FOOTBALL PLAYER GENERATIONS", id="label1"),
                            html.Label(
                                "Explore the differences between old-school and new talents",
                                className="label2",
                            ),
                            html.Br(),
                            html.Label(
                                "Dashboard created by: Catarina Pinheiro, Henrique Renda, Nguyen Phuc, Lorenzo Pigozzi",
                                className="label2",
                                style={"margin-bottom": ".34rem"},
                            ),
                        ],
                        width=8,
                    ),
                ],
                align="between",
                # no_gutters=True,
            ),
        ),
    ],
)

controls_player_1 = dbc.Card(
    [
        dbc.FormGroup(
            [
                html.Label("Choose an old-school Player:"),
                html.Br(),
                dropdown_player_over_25,
            ]
        ),
    ],
    body=True,
    className="controls_players",
)

controls_player_2 = dbc.Card(
    [
        dbc.FormGroup(
            [
                html.Label("Choose a new generation Player:"),
                html.Br(),
                dropdown_player_under_25,
            ]
        ),
    ],
    body=True,
    className="controls_players",
)

controls = dbc.Card(
    [
        dbc.FormGroup(
            [html.Label("Choose an Attribute:"), html.Br(), metric1_dropdown,]
        ),
        dbc.FormGroup(
            [html.Label("Choose an Attribute:"), html.Br(), metric2_dropdown,]
        ),
        dbc.FormGroup(
            [html.Label("Choose an Attribute:"), html.Br(), metric3_dropdown]
        ),
    ],
    body=True,
    className="controls",
)

controls_club = dbc.Card(
    [
        dbc.FormGroup(
            [html.Label("Choose a League:"), html.Br(), metric_club_dropdown1,]
        ),
        dbc.FormGroup(
            [
                html.Label("Choose an attribute for x:"),
                html.Br(),
                metric_scatter_dropdown1,
            ]
        ),
        dbc.FormGroup(
            [
                html.Label("Choose an attribute for y:"),
                html.Br(),
                metric_scatter_dropdown2,
            ]
        ),
    ],
    body=True,
    className="controls",
)

cards_1 = dbc.CardDeck(
    [
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Position", className="card-title1"),
                    html.Div(id="P_position1", className="card_info1"),
                ]
            ),
            className="attributes_card",
        ),
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Value", className="card-title1"),
                    html.Div(id="P_value1", className="card_info1"),
                ]
            ),
            className="attributes_card",
        ),
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Skill moves", className="card-title1"),
                    html.Div(id="P_skill1", className="card_info1"),
                ]
            ),
            className="attributes_card",
        ),
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Prefer foot", className="card-title1"),
                    html.Div(id="P_foot1", className="card_info1"),
                ]
            ),
            className="attributes_card",
        ),
    ]
)
cards_2 = dbc.CardDeck(
    [
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Prefer foot", className="card-title2"),
                    html.Div(id="P_foot2", className="card_info2"),
                ]
            ),
            className="attributes_card",
        ),
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Skill moves", className="card-title2"),
                    html.Div(id="P_skill2", className="card_info2"),
                ]
            ),
            className="attributes_card",
        ),
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Value", className="card-title2"),
                    html.Div(id="P_value2", className="card_info2"),
                ]
            ),
            className="attributes_card",
        ),
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Position", className="card-title2"),
                    html.Div(id="P_position2", className="card_info2"),
                ]
            ),
            className="attributes_card",
        ),
    ]
)
cards_3 = dbc.CardDeck(
    [
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Potential", className="card-title1"),
                    dcc.Graph(id="graph_example_1"),
                ]
            ),
            className="attributes_card",
        ),
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Skills", className="card-title1"),
                    dcc.Graph(id="graph_example_3"),
                ]
            ),
            className="attributes_card",
        ),
    ]
)
cards_4 = dbc.CardDeck(
    [
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Skills", className="card-title2"),
                    dcc.Graph(id="graph_example_4"),
                ]
            ),
            className="attributes_card",
        ),
        dbc.Card(
            dbc.CardBody(
                [
                    html.Div("Potential", className="card-title2"),
                    dcc.Graph(id="graph_example_2"),
                ]
            ),
            className="attributes_card",
        ),
    ]
)

tab1_content = (
    html.Div(
        [
            dbc.Card(
                dbc.CardBody(
                    [
                        html.H1("Player Comparison"),
                        html.Hr(),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Row(controls_player_1),
                                        dbc.Row(
                                            html.Img(
                                                src=app.get_asset_url("player_1.png"),
                                                className="playerImg",
                                            )
                                        ),
                                    ],
                                    sm=3,
                                ),
                                dbc.Col(
                                    dcc.Graph(id="graph_example"), sm=5, align="center"
                                ),
                                dbc.Col(
                                    [
                                        dbc.Row(controls_player_2),
                                        dbc.Row(
                                            html.Img(
                                                src=app.get_asset_url("player_3.png"),
                                                className="playerImg",
                                            )
                                        ),
                                    ],
                                    sm=3,
                                ),
                            ],
                            justify="between",
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dashtable_1,
                                        html.Br(),
                                        cards_1,
                                        html.Br(),
                                        cards_3,
                                    ],
                                    sm=6,
                                ),
                                dbc.Col(
                                    [
                                        dashtable_2,
                                        html.Br(),
                                        cards_2,
                                        html.Br(),
                                        cards_4,
                                    ],
                                    sm=6,
                                ),
                            ]
                        ),
                    ]
                )
            )
        ]
    ),
)


tab2_content = html.Div(
    [
        dbc.Card(
            dbc.CardBody(
                [
                    html.H1("League Analysis"),
                    html.Hr(),
                    dbc.Row(
                        [
                            dbc.Col(controls, sm=3),
                            dbc.Col(
                                dcc.Graph(
                                    id="league-graph1", className="LeagueBarPlot"
                                ),
                                sm=3,
                            ),
                            dbc.Col(
                                dcc.Graph(
                                    id="league-graph2", className="LeagueBarPlot"
                                ),
                                sm=3,
                            ),
                            dbc.Col(
                                dcc.Graph(
                                    id="league-graph3", className="LeagueBarPlot"
                                ),
                                sm=3,
                            ),
                        ],
                        align="center",
                    ),
                    dbc.Row(
                        [
                            dbc.Label(
                                "Select Age:",
                                style={"margin-left": "5%", "font-size": "20px"},
                            ),
                            dbc.Col(age_slider, align="center"),
                        ]
                    ),
                ]
            ),
        ),
        html.Br(),
        dbc.Card(
            dbc.CardBody(
                [
                    html.H1("Club Analysis"),
                    html.Hr(),
                    dbc.Row(
                        [
                            dbc.Col(controls_club, sm=2),
                            dbc.Col(
                                dcc.Graph(id="club-graph1", className="LeagueBarPlot"),
                                sm=5,
                            ),
                            dbc.Col(
                                dcc.Graph(id="club-graph2", className="LeagueBarPlot"),
                                sm=5,
                            ),
                        ],
                        align="center",
                    ),
                ]
            ),
            className=" mt-3",
        ),
    ]
)

app.layout = dbc.Container(
    [
        # html.H1("Fifa Players Analysis"),
        navbar,
        dbc.Tabs(
            [
                dbc.Tab(tab1_content, label="Players Comparison"),
                dbc.Tab(tab2_content, label="League & Club Analysis"),
            ],
        ),
    ],
    fluid=True,
)

#########Callbacks########################################


@app.callback(
    [
        Output(component_id="league-graph1", component_property="figure"),
        Output(component_id="league-graph2", component_property="figure"),
        Output(component_id="league-graph3", component_property="figure"),
    ],
    [
        Input(component_id="drop1", component_property="value"),
        Input(component_id="drop2", component_property="value"),
        Input(component_id="drop3", component_property="value"),
        Input(component_id="age_slider", component_property="value"),
    ],
)
###########Bar plot#######################################
def bar_plot(input_value1, input_value2, input_value3, age):

    filtered_by_age_data = data[(data["age"] >= age[0]) & (data["age"] <= age[1])]

    data_bar1 = dict(
        type="bar",
        y=filtered_by_age_data.groupby("league_name")
        .median()[input_value1]
        .sort_values(ascending=False)
        .head(5),
        x=filtered_by_age_data["league_name"].unique(),
    )

    data_bar2 = dict(
        type="bar",
        y=filtered_by_age_data.groupby("league_name")
        .median()[input_value2]
        .sort_values(ascending=False)
        .head(5),
        x=filtered_by_age_data["league_name"].unique(),
    )

    data_bar3 = (
        dict(
            type="bar",
            y=filtered_by_age_data.groupby("league_name")
            .median()[input_value3]
            .sort_values(ascending=False)
            .head(5),
            x=filtered_by_age_data["league_name"].unique(),
        ),
    )

    layout_bar1 = dict(xaxis=dict(title="League", tickangle=45),)

    layout_bar2 = dict(xaxis=dict(title="League", tickangle=45),)

    layout_bar3 = dict(xaxis=dict(title="League", tickangle=45),)
    fig1 = go.Figure(data=data_bar1, layout=layout_bar1)
    fig1.update_traces(
        marker_color="rgb(133,61,246)",
        marker_line_color="rgb(133,61,246)",
        marker_line_width=0.8,
        opacity=0.9,
    )
    fig1.update_layout(
        title_text=input_value1.capitalize(),
        title_x=0.5,
        margin=dict(l=70, r=40, t=60, b=40),
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        xaxis=dict(gridcolor="#e5e6dc", gridwidth=0.5),
        yaxis=dict(
            tickcolor="#e5e6dc", tickwidth=5, gridcolor="#e5e6dc", gridwidth=0.5
        ),
    )

    fig2 = go.Figure(data=data_bar2, layout=layout_bar2)
    fig2.update_traces(
        marker_color="rgb(158,50,249)",
        marker_line_color="rgb(133,61,246)",
        marker_line_width=0.8,
        opacity=0.9,
    )
    fig2.update_layout(
        title_text=input_value2.capitalize(),
        title_x=0.5,
        margin=dict(l=70, r=40, t=60, b=40),
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        xaxis=dict(gridcolor="#e5e6dc", gridwidth=0.5),
        yaxis=dict(
            tickcolor="#e5e6dc", tickwidth=5, gridcolor="#e5e6dc", gridwidth=0.5
        ),
    )
    fig3 = go.Figure(data=data_bar3, layout=layout_bar3)
    fig3.update_traces(
        marker_color="rgb(189,34,250)",
        marker_line_color="rgb(133,61,246)",
        marker_line_width=0.8,
        opacity=0.9,
    )
    fig3.update_layout(
        title_text=input_value3.capitalize(),
        title_x=0.5,
        margin=dict(l=70, r=40, t=60, b=40),
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        xaxis=dict(gridcolor="#e5e6dc", gridwidth=0.5),
        yaxis=dict(
            tickcolor="#e5e6dc", tickwidth=5, gridcolor="#e5e6dc", gridwidth=0.5
        ),
    )
    return fig1, fig2, fig3


# ----------------Callbacks for 2nd tab, clubs analysis----------------#
@app.callback(
    [
        Output(component_id="club-graph1", component_property="figure"),
        Output(component_id="club-graph2", component_property="figure"),
    ],
    [
        Input(component_id="club-drop1", component_property="value"),
        Input(component_id="scatter-drop1", component_property="value"),
        Input(component_id="scatter-drop2", component_property="value"),
    ],
)

###########Bar plot#######################################
def plots_clubs(league, x_val, y_val):
    # Scatter plot
    plot_df = (
        data[data["league_name"] == league]
        .sort_values("overall", ascending=False)
        .head(100)
    )
    fig1 = px.scatter(
        data_frame=plot_df,
        x=x_val,
        y=y_val,
        color="isOver25text",
        size="potential",
        color_discrete_sequence=["#5000bf", "rgb(255,171,0)"],
        hover_name="short_name",
        hover_data={
            "isOver25text": False,
            "age": True,
            "club_name": True,
            "overall": True,
            "potential": True,
            "player_positions": True,
        },
        title=("Top 100 players with highest overall rating in " + league),
    )
    fig1.update_layout(
        title=dict(font=dict(size=14)),
        legend=dict(
            title=dict(text="Age of players"),
            orientation="h",
            yanchor="bottom",
            y=1,
            xanchor="right",
            x=1,
        ),
        margin=dict(l=70, r=30, t=100, b=70),
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        xaxis=dict(gridcolor="#e5e6dc", gridwidth=0.5),
        yaxis=dict(
            tickcolor="#e5e6dc", tickwidth=5, gridcolor="#e5e6dc", gridwidth=0.5
        ),
    )
    fig1.update_traces(marker=dict(size=12, line=dict(width=2, color="DarkSlateGrey")))

    # Bar plot
    plot_df = (
        data[data["league_name"] == league]
        .groupby(["club_name", "isOver25"])
        .count()["short_name"]
        .reset_index()
    )
    plot_df = pd.pivot_table(
        plot_df,
        values="short_name",
        index=["club_name"],
        columns=["isOver25"],
        aggfunc=np.sum,
    ).reset_index()
    x = plot_df["club_name"]
    y = plot_df.iloc[:, -2:].div(plot_df.iloc[:, -2:].sum(axis=1), axis=0)
    fig2 = go.Figure()
    fig2.add_trace(
        go.Bar(
            y=x,
            x=y[True],
            name="Over 25",
            orientation="h",
            marker=dict(color="#5000bf", line=dict(color="DarkSlateGrey", width=1.1)),
        )
    )
    fig2.add_trace(
        go.Bar(
            y=x,
            x=y[False],
            name="Under 25",
            orientation="h",
            marker=dict(
                color="rgb(255,171,0)", line=dict(color="DarkSlateGrey", width=1.1)
            ),
        )
    )

    fig2.update_yaxes(tickfont=dict(size=10))
    fig2.update_layout(
        barmode="stack",
        title=dict(
            text="Players over 25 and under 25 years old by club", font=dict(size=14)
        ),
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        legend=dict(
            title=dict(text="Age of players"),
            orientation="h",
            yanchor="bottom",
            y=1,
            xanchor="right",
            x=0.95,
        ),
        margin=dict(r=30, t=100, b=70),
    )
    return fig1, fig2


# ----------------Callbacks for 1st tab, clubs analysis----------------#


@app.callback(
    [
        Output("graph_example", "figure"),
        Output("table1", "data"),
        Output("graph_example_1", "figure"),
        Output("graph_example_3", "figure"),
        Output("table2", "data"),
        Output("graph_example_2", "figure"),
        Output("graph_example_4", "figure"),
        Output("P_position1", "children"),
        Output("P_value1", "children"),
        Output("P_skill1", "children"),
        Output("P_foot1", "children"),
        Output("P_position2", "children"),
        Output("P_value2", "children"),
        Output("P_skill2", "children"),
        Output("P_foot2", "children"),
    ],
    [Input("player1", "value"), Input("player2", "value")],
)

###############################################   radar plot   #####################################################


def tab_1_function(player1, player2):

    # scatterpolar
    df1_for_plot = pd.DataFrame(df1[df1["long_name"] == player1][skill_player].iloc[0])
    df1_for_plot.columns = ["score"]
    df2_for_plot = pd.DataFrame(df2[df2["long_name"] == player2][skill_player].iloc[0])
    df2_for_plot.columns = ["score"]
    list_scores = [
        df1_for_plot.index[i].capitalize() + " = " + str(df1_for_plot["score"][i])
        for i in range(len(df1_for_plot))
    ]
    text_scores_1 = player1
    for i in list_scores:
        text_scores_1 += "<br>" + i

    list_scores = [
        df2_for_plot.index[i].capitalize() + " = " + str(df2_for_plot["score"][i])
        for i in range(len(df2_for_plot))
    ]
    text_scores_2 = player2
    for i in list_scores:
        text_scores_2 += "<br>" + i

    fig = go.Figure(
        data=go.Scatterpolar(
            r=df1_for_plot["score"],
            theta=df1_for_plot.index,
            fill="toself",
            marker_color="rgb(45,0,198)",
            opacity=1,
            hoverinfo="text",
            name=text_scores_1,
            text=[
                df1_for_plot.index[i] + " = " + str(df1_for_plot["score"][i])
                for i in range(len(df1_for_plot))
            ],
        )
    )
    fig.add_trace(
        go.Scatterpolar(
            r=df2_for_plot["score"],
            theta=df2_for_plot.index,
            fill="toself",
            marker_color="rgb(255,171,0)",
            hoverinfo="text",
            name=text_scores_2,
            text=[
                df2_for_plot.index[i] + " = " + str(df2_for_plot["score"][i])
                for i in range(len(df2_for_plot))
            ],
        )
    )

    fig.update_layout(
        polar=dict(
            hole=0.1,
            bgcolor="white",
            radialaxis=dict(
                visible=True,
                type="linear",
                autotypenumbers="strict",
                autorange=False,
                range=[30, 100],
                angle=90,
                showline=False,
                showticklabels=False,
                ticks="",
                gridcolor="black",
            ),
        ),
        width=550,
        height=550,
        margin=dict(l=80, r=80, t=20, b=20),
        showlegend=False,
        template="plotly_dark",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        font_color="black",
        font_size=15,
    )

    # table 1
    table_updated1 = df[df["long_name"] == player1].to_dict("records")

    # gauge plot 1
    df1_for_plot = pd.DataFrame(df1[df1["long_name"] == player1]["potential"])
    df1_for_plot["name"] = player2
    gauge1 = go.Figure(
        go.Indicator(
            domain={"x": [0, 1], "y": [0, 1]},
            value=df1_for_plot.potential.iloc[0],
            mode="gauge+number",
            gauge={"axis": {"range": [None, 100]}, "bar": {"color": "#5000bf"}},
        )
    )
    gauge1.update_layout(
        height=300,
        margin=dict(l=10, r=10, t=40, b=10),
        showlegend=False,
        template="plotly_dark",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        font_color="black",
        font_size=15,
    )
    # barplot 1
    df1_for_plot = pd.DataFrame(
        df1[df1["long_name"] == player1][skills1].iloc[0].reset_index()
    )
    df1_for_plot.rename(columns={df1_for_plot.columns[1]: "counts"}, inplace=True)
    df1_for_plot.rename(columns={df1_for_plot.columns[0]: "skills"}, inplace=True)
    barplot1 = px.bar(df1_for_plot, x="skills", y="counts")
    barplot1.update_traces(marker_color="#5000bf")
    barplot1.update_layout(
        height=300,
        margin=dict(l=10, r=10, t=20, b=0),
        showlegend=False,
        # yaxis={'visible': False, 'showticklabels': True},
        template="plotly_dark",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        font_color="black",
        font_size=10,
    )
    barplot1.update_yaxes(range=[1, 100])
    # table 2
    table_updated2 = df[df["long_name"] == player2].to_dict("records")

    # gauge plot 2
    df2_for_plot = pd.DataFrame(df2[df2["long_name"] == player2]["potential"])
    df2_for_plot["name"] = player2
    gauge2 = go.Figure(
        go.Indicator(
            domain={"x": [0, 1], "y": [0, 1]},
            value=df2_for_plot.potential.iloc[0],
            mode="gauge+number",
            gauge={"axis": {"range": [None, 100]}, "bar": {"color": "rgb(255,171,0)"}},
        )
    )
    gauge2.update_layout(
        height=300,
        margin=dict(l=10, r=10, t=40, b=10),
        showlegend=False,
        template="plotly_dark",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        font_color="black",
        font_size=15,
    )
    # bar plot 2
    df2_for_plot = pd.DataFrame(
        df2[df2["long_name"] == player2][skills1].iloc[0].reset_index()
    )
    df2_for_plot.rename(columns={df2_for_plot.columns[1]: "counts"}, inplace=True)
    df2_for_plot.rename(columns={df2_for_plot.columns[0]: "skills"}, inplace=True)
    barplot2 = px.bar(df2_for_plot, x="skills", y="counts")
    barplot2.update_traces(marker_color="rgb(255,171,0)")
    barplot2.update_layout(
        height=300,
        margin=dict(l=10, r=10, t=20, b=0),
        showlegend=False,
        template="plotly_dark",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        font_color="black",
        font_size=10,
    )
    barplot2.update_yaxes(range=[1, 100])
    # cards
    p_pos_1 = df1[df1["long_name"] == player1]["team_position"]
    p_value_1 = (
        str(df1[df1["long_name"] == player1]["value_eur"].values[0] / 1000000)
        + " M Euro"
    )
    p_skill_1 = df1[df1["long_name"] == player1]["skill_moves"]
    p_foot_1 = df1[df1["long_name"] == player1]["preferred_foot"]

    p_pos_2 = df2[df2["long_name"] == player2]["team_position"]
    p_value_2 = (
        str(df2[df2["long_name"] == player2]["value_eur"].values[0] / 1000000)
        + " M Euro"
    )
    p_skill_2 = df2[df2["long_name"] == player2]["skill_moves"]
    p_foot_2 = df2[df2["long_name"] == player2]["preferred_foot"]

    # outputs
    return (
        fig,
        table_updated1,
        gauge1,
        barplot1,
        table_updated2,
        gauge2,
        barplot2,
        p_pos_1,
        p_value_1,
        p_skill_1,
        p_foot_1,
        p_pos_2,
        p_value_2,
        p_skill_2,
        p_foot_2,
    )


if __name__ == "__main__":
    app.run_server(debug=True)
