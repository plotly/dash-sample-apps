# Imports
import dash
import ast

import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
from whitenoise import WhiteNoise

from chessboard import getChessboard, getHeatmap, getStackedBar, getBoard
from styles import *

# Read the .csv file with the preprocessed data.
url = "https://raw.githubusercontent.com/Exileus/DataVis2021_proj2/main/chess_app.csv"
df_original = pd.read_csv(
    url,
    sep=",",
    dtype={"pawns": int, "knights": int, "bishops": int, "rooks": int, "queens": int},
    converters={
        "wKing_sqr": ast.literal_eval,
        "bKing_sqr": ast.literal_eval,
        "wQueen_sqr": ast.literal_eval,
        "bQueen_sqr": ast.literal_eval,
        "wRook_sqr": ast.literal_eval,
        "bRook_sqr": ast.literal_eval,
        "wRook2_sqr": ast.literal_eval,
        "bRook2_sqr": ast.literal_eval,
        "wBishop_sqr": ast.literal_eval,
        "bBishop_sqr": ast.literal_eval,
        "wBishop2_sqr": ast.literal_eval,
        "bBishop2_sqr": ast.literal_eval,
        "wKnight_sqr": ast.literal_eval,
        "bKnight_sqr": ast.literal_eval,
        "wKnight2_sqr": ast.literal_eval,
        "bKnight2_sqr": ast.literal_eval,
    },
)

# Calculate min and max elo
min_elo, max_elo = df_original["avg_Elo"].min(), df_original["avg_Elo"].max()
max_moves = df_original["moves"].max()

# Define function to output an 8*8 dataframe based on a df and a list of column names to parse.


def board_output(df, col_list):
    brd = np.zeros((8, 8))
    for col_name in col_list:
        for tup in df[col_name]:
            if tup == (None, None):
                pass
            else:
                brd[tup] += 1

    return pd.DataFrame(brd)


# Define global variables for later.
g_color = "white_color"
g_piece = "King"
g_status, g_winner, g_time_control, g_game_type = ".*", ".*", ".*", ".*"
pieces_list = ["King", "Queen", "Rook", "Bishop", "Knight"]
# Define a dictionary to be used to update the board with the correct columns.
color_piece_dict = cp_dict = {
    ("white_color", "King"): ["wKing_sqr"],
    ("black_color", "King"): ["bKing_sqr"],
    ("white_color", "Queen"): ["wQueen_sqr"],
    ("black_color", "Queen"): ["bQueen_sqr"],
    ("white_color", "Rook"): ["wRook_sqr", "wRook2_sqr"],
    ("black_color", "Rook"): ["bRook_sqr", "bRook2_sqr"],
    ("white_color", "Bishop"): ["wBishop_sqr", "wBishop2_sqr"],
    ("black_color", "Bishop"): ["bBishop_sqr", "bBishop2_sqr"],
    ("white_color", "Knight"): ["wKnight_sqr", "wKnight2_sqr"],
    ("black_color", "Knight"): ["bKnight_sqr", "bKnight2_sqr"],
}

# Define an additional dict for dropdown status to use for callbacks.
dropdown_status_dict = st_dict = {
    "st_all": ".*",
    "st_draw": "draw",
    "st_mate": "mate",
    "st_resign": "resign",
    "st_outoftime": "outoftime",
}

dropdown_winner_dict = wn_dict = {
    "wn_all": ".*",
    "wn_white": "white",
    "wn_black": "black",
}


dropdown_time_control_dict = tc_dict = {
    "tc_all": ".*",
    "tc_bullet": "Bullet",
    "tc_blitz": "Blitz",
    "tc_classic": "Classical",
    "tc_none": "Correspondence",
}

dropdown_game_type_dict = gt_dict = {
    "gt_all": ".*",
    "gt_std": "game",
    "gt_tourney": "tournament",
}


popover_status = dbc.Popover(
    [
        dbc.PopoverHeader("Status of the Game"),
        dbc.PopoverBody(
            "Games can be over in a myriad of ways, either by checkmate, draw, player resignation, or when a player runs out of time. Filter the games by these conditions here."
        ),
    ],
    trigger="hover",
    target="dropdown_status",
    placement="left",
)

popover_time_control = dbc.Popover(
    [
        dbc.PopoverHeader("Time Control Filter"),
        dbc.PopoverBody(
            "Players have a specific time to make their moves. The games in the dataset follow this convention: Bullet Games (0-3 minutes), Blitz(3-10 minutes), Classical(10 minutes+). Note: Lichess uses a slight different system today."
        ),
    ],
    trigger="hover",
    target="dropdown_time_control",
    placement="left",
)

popover_game_type = dbc.Popover(
    [
        dbc.PopoverHeader("Type of Competitive Setting"),
        dbc.PopoverBody(
            "This dataset contains games played in specific tournaments, hosted by Lichess."
        ),
    ],
    trigger="hover",
    target="dropdown_game_type",
    placement="left",
)

about_this = html.Div(
    [
        dbc.Button("About this Visualization", id="abt_us"),
        dbc.Popover(
            [
                dbc.PopoverHeader("Powered by Lichess"),
                dbc.PopoverBody(
                    "This visualization is powered by a dataset of games played in April, 2017, sourced from the publically available lichess database.\nAuthors:Frederico Santos, Rupesh Baradi, Tiago Ramos.\nNova IMS,Data Visualization Course, 2021."
                ),
            ],
            trigger="click",
            target="abt_us",
        ),
        "stuff",
    ]
)


# Set stylesheets and app.
# ["https://codepen.io/chriddyp/pen/bWLwgP.css"]
FA = "https://use.fontawesome.com/releases/v5.12.1/css/all.css"
external_stylesheets = [dbc.themes.LUX, FA]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "CHESS KINGDOM"
server = app.server
server.wsgi_app = WhiteNoise(server.wsgi_app, root="static/")


# Defining app layout
margin_bottom = "30px"

# Banner

banner = dbc.Row(
    children=[
        dbc.Col(
            html.Img(
                src="/assets/apple-touch-icon.png",
                id="logo",
                style={"border-radius": "50%"},
            ),
            width=2,
            align="left",
        ),
        dbc.Col(
            html.H1("A Visualization of Endgame Chess Pieces"),
            align="center",
            width=10,
        ),
    ],
    style={"margin-bottom": "50px", "margin-top": "-30px"},
    align="center",
)

# Graph
graph = dbc.Row(
    style={"margin-left": "auto", "margin-right": "auto"},
    children=[
        dcc.Graph(
            id="chessboard",
            style={"margin-left": "auto", "margin-right": "auto"},
            config={
                "displayModeBar": False,
                "scrollZoom": False,
                "showAxisDragHandles": False,
            },
        )
    ],
)

# Stacked Bar
stacked_graph = dbc.Row(
    style={"margin-bottom": "30px"},
    justify="center",
    children=[
        dcc.Graph(
            id="stackedbar",
            animate=True,
            config={
                "displayModeBar": False,
                "scrollZoom": False,
                "showAxisDragHandles": False,
            },
        )
    ],
)


text_margin = "6px"

c_total_games = dbc.Row(
    style={"margin-bottom": "20px"},
    justify="center",
    children=[
        dbc.Col(
            children=[
                html.Div(id="game_count", style={"text-align": "center"}),
                html.Div(
                    "TOTAL GAMES",
                    style={"margin-left": text_margin, "text-align": "center"},
                ),
            ],
        ),
        dbc.Col(
            children=[
                html.Div(id="white_wins", style={"text-align": "center"}),
                html.Div(
                    "WINS BY WHITE",
                    style={"margin-left": text_margin, "text-align": "center"},
                ),
            ],
        ),
        dbc.Col(
            children=[
                html.Div(id="black_wins", style={"text-align": "center"}),
                html.Div(
                    "WINS BY BLACK",
                    style={"margin-left": text_margin, "text-align": "center"},
                ),
            ],
        ),
        dbc.Col(
            children=[
                html.Div(id="draw", style={"text-align": "center"}),
                html.Div(
                    "DRAWS", style={"margin-left": text_margin, "text-align": "center"}
                ),
            ],
        ),
    ],
)
# BLACK / WHITE
c_choose_side = dbc.Col(
    style={"margin-bottom": margin_bottom},
    children=[
        html.Div(
            str("Choose side").upper(),
            style={"text-align": "center", "margin-bottom": text_margin},
        ),
        dbc.Row(
            justify="center",
            children=[
                dbc.ButtonGroup(
                    style={"text-align": "center"},
                    children=[
                        dbc.Button(
                            "White",
                            color="secondary",
                            n_clicks=0,
                            id="white_color",
                            outline=True,
                            active=True,
                        ),
                        dbc.Button(
                            "Black",
                            color="dark",
                            n_clicks=0,
                            id="black_color",
                            outline=True,
                            active=False,
                        ),
                    ],
                ),
            ],
        ),
    ],
)

c_select_piece = dbc.Col(
    style={"margin-bottom": margin_bottom},
    width=9,
    children=[
        html.Div(
            str("Select Piece").upper(),
            style={"text-align": "center", "margin-bottom": text_margin},
        ),
        dbc.Row(
            justify="center",
            children=[
                dbc.ButtonGroup(
                    children=[
                        dbc.Button(
                            [
                                html.I(className=f"fas fa-chess-{name.lower()} mr-2"),
                                name,
                            ],
                            color="primary",
                            n_clicks=0,
                            outline=True,
                            id=name,
                            active=False,
                        )
                        for name in pieces_list
                    ],
                )
            ],
        ),
    ],
)
c_elo_slider = dbc.Col(
    style={
        "margin-bottom": margin_bottom,
        "margin-left": "auto",
        "margin-right": "auto",
    },
    width=12,
    children=[
        html.Div(
            str("Elo range").upper(),
            style={"text-align": "center", "margin-bottom": text_margin},
        ),
        dcc.RangeSlider(
            id="elo_slider",
            min=min_elo,
            max=max_elo,
            value=[min_elo, max_elo],
            step=10,
            pushable=1,
            allowCross=False,
            marks={
                i: str(i)
                for i in range(
                    int(min_elo) - 1,
                    int(max_elo) + 1,
                    int((max_elo - min_elo + 2) // 10),
                )
            },
        ),
    ],
)
c_moves_slider = dbc.Col(
    style={
        "margin-bottom": margin_bottom,
        "margin-left": "auto",
        "margin-right": "auto",
    },
    width=12,
    children=[
        html.Div(
            str("Number of Moves").upper(),
            style={"text-align": "center", "margin-bottom": text_margin},
        ),
        dcc.RangeSlider(
            id="moves_slider",
            min=1,
            max=max_moves,
            value=[0, max_moves],
            step=1,
            pushable=1,
            allowCross=False,
            marks={i: str(i) for i in range(0, max_moves, 5)},
        ),
    ],
)

dropdown_status = dbc.DropdownMenu(
    [
        dbc.DropdownMenuItem("Status", header=True),
        dbc.DropdownMenuItem("All", id="st_all", n_clicks=0),
        dbc.DropdownMenuItem("Draws", id="st_draw", n_clicks=0),
        dbc.DropdownMenuItem("Checkmate", id="st_mate", n_clicks=0),
        dbc.DropdownMenuItem("Resignation", id="st_resign", n_clicks=0),
        dbc.DropdownMenuItem("Time Forfeit", id="st_outoftime", n_clicks=0),
    ],
    label="Status",
    id="dropdown_status",
)

dropdown_winner = dbc.Collapse(
    dbc.DropdownMenu(
        [
            dbc.DropdownMenuItem("Winning Side", header=True),
            dbc.DropdownMenuItem("All", id="wn_all", n_clicks=0),
            dbc.DropdownMenuItem("White", id="wn_white", n_clicks=0),
            dbc.DropdownMenuItem("Black", id="wn_black", n_clicks=0),
        ],
        label="Winning Side",
        id="dropdown_winner",
    ),
    id="wn_menu",
)

dropdown_time_control = dbc.DropdownMenu(
    [
        dbc.DropdownMenuItem("Time Control", header=True),
        dbc.DropdownMenuItem("All", id="tc_all", n_clicks=0),
        dbc.DropdownMenuItem("Bullet", id="tc_bullet", n_clicks=0),
        dbc.DropdownMenuItem("Blitz", id="tc_blitz", n_clicks=0),
        # dbc.DropdownMenuItem("Rapid",id="tc_rpd",n_clicks=0), if this shows up later then include it.
        dbc.DropdownMenuItem("Classical", id="tc_classic", n_clicks=0),
        dbc.DropdownMenuItem("No Time Control", id="tc_none", n_clicks=0),
    ],
    label="Time Control",
    id="dropdown_time_control",
)

dropdown_game_type = dbc.DropdownMenu(
    [
        dbc.DropdownMenuItem("Game Type", header=True),
        dbc.DropdownMenuItem("All", id="gt_all", n_clicks=0),
        dbc.DropdownMenuItem("Standard", id="gt_std", n_clicks=0),
        dbc.DropdownMenuItem("Tournament", id="gt_tourney", n_clicks=0),
    ],
    label="Game Type",
    id="dropdown_game_type",
)

dropdown_menus = dbc.Row(
    style={"margin-bottom": margin_bottom},
    justify="center",
    children=[
        dropdown_status,
        popover_status,
        dropdown_winner,
        dropdown_time_control,
        popover_time_control,
        dropdown_game_type,
        popover_game_type,
    ],
)

app.layout = dbc.Jumbotron(
    style={"background-color": "#ebebeb"},  # ADD SETTINGS HERE
    children=[
        # Banner
        # Main Layout
        dbc.Row(  # ADD SETTINGS HERE
            children=[
                # PARAMETER SETTINGS COLUMN
                dbc.Col(
                    children=[
                        banner,
                        c_total_games,
                        stacked_graph,
                        dbc.Row(
                            style={"margin-bottom": margin_bottom},
                            children=[c_choose_side, c_select_piece],
                        ),
                        c_elo_slider,
                        c_moves_slider,
                        dropdown_menus,
                        about_this,
                    ]
                ),
                # CHESS BOARD COLUMN
                dbc.Col(width={"size": 6}, children=[graph]),
            ],
        ),
    ],
)


@app.callback(
    Output("chessboard", "figure"),
    Output("stackedbar", "figure"),
    Output("game_count", "children"),
    Output("white_wins", "children"),
    Output("black_wins", "children"),
    Output("draw", "children"),
    Output("wn_menu", "is_open"),
    Output("white_color", "active"),
    Output("black_color", "active"),
    Output("King", "active"),
    Output("Queen", "active"),
    Output("Rook", "active"),
    Output("Bishop", "active"),
    Output("Knight", "active"),
    Input("white_color", "n_clicks"),
    Input("black_color", "n_clicks"),
    Input("King", "n_clicks"),
    Input("Queen", "n_clicks"),
    Input("Rook", "n_clicks"),
    Input("Bishop", "n_clicks"),
    Input("Knight", "n_clicks"),
    Input("elo_slider", "value"),
    Input("moves_slider", "value"),
    Input("st_all", "n_clicks"),
    Input("st_draw", "n_clicks"),
    Input("st_mate", "n_clicks"),
    Input("st_resign", "n_clicks"),
    Input("st_outoftime", "n_clicks"),
    Input("wn_all", "n_clicks"),
    Input("wn_white", "n_clicks"),
    Input("wn_black", "n_clicks"),
    Input("tc_all", "n_clicks"),
    Input("tc_blitz", "n_clicks"),
    Input("tc_bullet", "n_clicks"),
    Input("tc_classic", "n_clicks"),
    Input("tc_none", "n_clicks"),
    Input("gt_all", "n_clicks"),
    Input("gt_std", "n_clicks"),
    Input("gt_tourney", "n_clicks"),
)
def update_chessboard(
    white_color,
    black_color,
    King,
    Queen,
    Rook,
    Bishop,
    Knight,
    elo_range,
    move_range,
    st_all,
    st_draw,
    st_mate,
    st_resign,
    st_outoftime,
    wn_all,
    wn_white,
    wn_black,
    tc_all,
    tc_blitz,
    tc_bullet,
    tc_classic,
    tc_none,
    gt_all,
    gt_std,
    gt_tourney,
):
    # Trigger button here, for when a button is pressed.
    trigger_button = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

    global g_status
    global g_winner
    global g_time_control
    global g_game_type

    if trigger_button in st_dict.keys():
        g_status = st_dict[trigger_button]

    elif trigger_button in wn_dict.keys():
        g_winner = wn_dict[trigger_button]

    elif trigger_button in tc_dict.keys():
        g_time_control = tc_dict[trigger_button]

    elif trigger_button in gt_dict.keys():
        g_game_type = gt_dict[trigger_button]

    # Filters go here.
    dff = df_original[
        (df_original["avg_Elo"] >= int(elo_range[0]))
        & (df_original["avg_Elo"] <= int(elo_range[1]))
        & (df_original["moves"] >= int(move_range[0]))
        & (df_original["moves"] <= int(move_range[-1]))
        & (df_original["victory_status"].str.contains(g_status))
        & (df_original["Winner"].str.contains(g_winner))
        & (df_original["Event"].str.contains(g_time_control))
        & (df_original["Event"].str.contains(g_game_type))
    ]

    # Before further manipulation, get the number of games from the filtered dataframe.
    game_count = dff.shape[0]
    game_results = dff.Winner.value_counts().to_dict()
    game_results_norm = np.round(
        dff.Winner.str.upper().value_counts(normalize=True), 2
    ).to_dict()

    if "white" in game_results.keys():
        white_wins = game_results["white"]
    else:
        white_wins = 0
    if "black" in game_results.keys():
        black_wins = game_results["black"]
    else:
        black_wins = 0
    if "draw" in game_results.keys():
        draw = game_results["draw"]
    else:
        draw = 0
    stackedbar = getStackedBar(game_results_norm)

    # Then retrieve the column of interest.
    global g_color
    global g_piece

    if trigger_button in ["white_color", "black_color"]:
        g_color = trigger_button
    if trigger_button in pieces_list:
        g_piece = trigger_button

    df = board_output(dff, cp_dict[g_color, g_piece])

    # Additionally:
    if g_status == "draw":
        is_open = False
    else:
        is_open = True
    # Additionaly pt.2:
    if g_color == "white_color":
        wc_act, bc_act = True, False
    else:
        wc_act, bc_act = False, True

    # Additionaly pt3:
    k_act, q_act, r_act, b_act, n_act = [x == g_piece for x in pieces_list]

    # Transform it for the heatmap.
    df = (
        df.stack()
        .reset_index()
        .rename(columns={"level_0": "rows", "level_1": "cols", 0: "freq"})
    )

    df["rows"] = df["rows"].replace({i: list(range(8))[::-1][i] for i in range(8)})
    chessboard = getChessboard(800)
    chessboard.add_trace(getHeatmap(dataframe=df))

    return (
        chessboard,
        stackedbar,
        game_count,
        white_wins,
        black_wins,
        draw,
        is_open,
        wc_act,
        bc_act,
        k_act,
        q_act,
        r_act,
        b_act,
        n_act,
    )


# Statring the dash app
if __name__ == "__main__":
    app.run_server(debug=True)
