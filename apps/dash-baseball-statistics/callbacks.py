# import dash IO and graph objects
from dash.dependencies import Input, Output

# Plotly graph objects to render graph plots
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Import dash html, bootstrap components, and tables for datatables
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_table

# Import app
from app import app

# Import custom data.py
import data

# Import data from data.py file
teams_df = data.teams
# Hardcoded list that contain era names
era_list = data.era_list

# Player Profiles data
player_df = data.players
team_players_df = data.team_players

# Statistical data
batter_df = data.batters
field_df = data.fielding
pitching_df = data.pitching


# This will update the team dropdown and the range of the slider
@app.callback(
    [
        Output("team-dropdown", "options"),
        Output("team-dropdown", "value"),
        Output("era-slider", "value"),
    ],
    [Input("era-dropdown", "value")],
)
def select_era(selected_era):
    # Check if selected era is equal to the value in the era list
    # Makes sure that teams and range are set to desired era
    if selected_era == era_list[1]["value"]:
        teams = data.dynamicteams(1)
        range = data.dynamicrange(1)
    elif selected_era == era_list[2]["value"]:
        teams = data.dynamicteams(2)
        range = data.dynamicrange(2)
    elif selected_era == era_list[3]["value"]:
        teams = data.dynamicteams(3)
        range = data.dynamicrange(3)
    elif selected_era == era_list[4]["value"]:
        teams = data.dynamicteams(4)
        range = data.dynamicrange(4)
    elif selected_era == era_list[5]["value"]:
        teams = data.dynamicteams(5)
        range = data.dynamicrange(5)
    elif selected_era == era_list[6]["value"]:
        teams = data.dynamicteams(6)
        range = data.dynamicrange(6)
    elif selected_era == era_list[7]["value"]:
        teams = data.dynamicteams(7)
        range = data.dynamicrange(7)
    else:
        teams = data.dynamicteams(0)
        range = data.dynamicrange(0)
    # Return team list, the initial value of the team list, and the range in the era
    return teams, teams[0]["value"], range


# Team championship datatable
@app.callback(
    [Output("team-data", "children")],
    [Input("team-dropdown", "value"), Input("era-slider", "value")],
)
def update_win_table(selected_team, year_range):
    # Create filter dataframe of requested team
    filter_team = teams_df[teams_df.team_id == selected_team]
    # I will revisit this again soon, it just doesnt seem efficient
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]

    # Filter unneccessary data
    Data = filter_year.drop(
        columns=[
            "team_id",
            "franchise_id",
            "div_id",
            "ghome",
            "g",
            "w",
            "l",
            "r",
            "ab",
            "h",
            "double",
            "triple",
            "hr",
            "bb",
            "so",
            "sb",
            "cs",
            "era",
            "cg",
            "sho",
            "sv",
            "ha",
            "hra",
            "bba",
            "soa",
            "e",
            "dp",
            "hbp",
            "sf",
            "ra",
            "er",
            "ipouts",
            "fp",
            "name",
            "park",
            "attendance",
            "bpf",
            "ppf",
            "team_id_br",
            "team_id_lahman45",
            "team_id_retro",
        ]
    )

    # Set data to variable if team won world series
    WIN = Data[Data.ws_win == "Y"]
    # Check if dataframe is empty (no world series won)
    if WIN.empty:
        # Set data to variable if team won a wild card
        WIN = Data[Data.wc_win == "Y"]
        # if the team did not win a wild card (dataframe is empty)
        if WIN.empty:
            # finally set data to variable if team won a division title
            WIN = Data[Data.div_win == "Y"]

    # Create empty data list for notification if needed
    data_note = []

    # If the team did not win any championships append alert and return data list
    if WIN.empty:
        data_note.append(
            html.Div(
                dbc.Alert(
                    "The selected team did not win any championships.", color="warning"
                )
            )
        )
        return data_note
    # else set and return datatable with team championship data
    else:
        data_note.append(
            html.Div(
                dash_table.DataTable(
                    data=WIN.to_dict("records"),
                    columns=[{"name": x, "id": x} for x in WIN],
                    style_as_list_view=True,
                    editable=False,
                    style_table={
                        "overflowY": "scroll",
                        "width": "100%",
                        "minWidth": "100%",
                    },
                    style_header={"backgroundColor": "#f8f5f0", "fontWeight": "bold"},
                    style_cell={"textAlign": "center", "padding": "8px"},
                )
            )
        )
        return data_note


# Callback to a W-L Bar Chart, takes data request from dropdown
@app.callback(
    Output("wl-bar", "figure"),
    [Input("team-dropdown", "value"), Input("era-slider", "value")],
)
def update_figure1(selected_team, year_range):
    filter_team = teams_df[teams_df.team_id == selected_team]
    # This feels like a hack
    # Checks if year_range is empty (NonType)
    if year_range:
        # Filter the years of the data to be within range
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        # I created this becuase i kept getting NonType errors with all of my graph call backs
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]
    # Create Bar Chart figure, Wins and Losses
    fig1 = go.Figure(
        data=[
            go.Bar(
                name="Wins",
                x=filter_year.year,
                y=filter_year.w,
                marker_color="#004687",
                opacity=0.8,
            ),
            go.Bar(
                name="Losses",
                x=filter_year.year,
                y=filter_year.l,
                marker_color="#AE8F6F",
                opacity=0.8,
            ),
        ]
    )
    # set x axes title and tick to only include year given no half year such as 1927.5
    fig1.update_xaxes(title="Year", tickformat="d")
    # set y axes to fixed selection range, user can only select data in the x axes
    fig1.update_yaxes(fixedrange=True)
    # Update figure, set hover to the X-Axis and establish title
    fig1.update_layout(
        hovermode="x",
        barmode="group",
        title="Win/Loss Performance",
        font={"color": "darkslategray"},
        paper_bgcolor="white",
        plot_bgcolor="#f8f5f0",
        legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1),
    )
    # return figure
    return fig1


# Call back to Line Graph, takes data request from dropdown
@app.callback(
    Output("batting-line", "figure"),
    [Input("team-dropdown", "value"), Input("era-slider", "value")],
)
def update_figure2(selected_team, year_range):
    # Create filter dataframe of requested team data
    filter_team = teams_df[teams_df.team_id == selected_team]
    # This still feels like a hack
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]

    # Calculate Slugging Average
    SLG = data.calculate_slg(filter_year)
    # Calculate BABIP
    BABIP = (filter_year.h - filter_year.hr) / (
        filter_year.ab - filter_year.so - filter_year.hr
    )
    # Calculete Batting Average
    # BAVG = filter_year.h / filter_year.ab

    # Create Line char figure using Slugging Average, BABIP, and Batting Average
    fig2 = go.Figure(
        data=[
            go.Scatter(
                name="Slugging Average",
                x=filter_year.year,
                y=SLG,
                mode="lines+markers",
                marker_color="Orange",
                opacity=0.9,
            ),
            go.Scatter(
                name="Batting Average Balls In Play",
                x=filter_year.year,
                y=BABIP,
                mode="lines+markers",
                marker_color="#005C5C",
                opacity=0.8,
            ),
            # go.Scatter(name='Batting Average', x=filter_year.year, y=BAVG, mode='lines+markers', marker_color='#0C2C56',opacity=0.8),
        ]
    )

    fig2.update_xaxes(title="Year", tickformat="d")
    fig2.update_yaxes(fixedrange=True)
    fig2.update_layout(
        hovermode="x",
        title="Batting Performance",
        font={"color": "darkslategray"},
        paper_bgcolor="white",
        plot_bgcolor="#f8f5f0",
        legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1),
    )

    return fig2


# Call back to Line Chart, Takes request data from dropdown menu
@app.callback(
    Output("feild-line", "figure"),
    [Input("team-dropdown", "value"), Input("era-slider", "value")],
)
def update_figure3(selected_team, year_range):
    # Create filter dataframe of requested team data
    filter_team = teams_df[teams_df.team_id == selected_team]
    # Im pretty sure this is a hack, there has to be a different approch
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]
    # Create line chart figure using Errors and Double Plays
    fig3 = go.Figure(
        data=[
            go.Bar(
                name="Errors",
                x=filter_year.year,
                y=filter_year.e,
                marker=dict(color="#5F259F"),
                opacity=0.7,
            ),
            go.Bar(
                name="Double Plays",
                x=filter_year.year,
                y=filter_year.dp,
                marker=dict(color="#005F61"),
                opacity=0.7,
            ),
        ]
    )

    fig3.update_xaxes(title="Year", tickformat="d")
    fig3.update_yaxes(fixedrange=True)
    fig3.update_layout(
        barmode="stack",
        hovermode="x",
        title="Feilding Performance",
        font={"color": "darkslategray"},
        paper_bgcolor="white",
        plot_bgcolor="#f8f5f0",
        legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1),
    )

    return fig3


# Call back to Pie Chart, takes data request from dropdown menu
@app.callback(
    Output("pitch-pie", "figure"),
    [Input("team-dropdown", "value"), Input("era-slider", "value")],
)
def update_figure4(selected_team, year_range):
    # Create filter dataframe of requested team data
    filter_team = teams_df[teams_df.team_id == selected_team]
    # IDK what to say, this is needed until i can investigate more.
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]
    # Create lists of team data: games, completed games (by single pitcher), shutouts, saves
    # Prevents a RuntimeError: invalid value encountered in longlong_scalars
    if filter_year.g.sum() == 0:
        G = 1
    else:
        G = filter_year.g.sum()

    CG = filter_year.cg.sum() / G
    SHO = filter_year.sho.sum() / G
    SV = filter_year.sv.sum() / G
    # Set an embedded list
    PCT = [CG, SHO, SV]

    # Create Pie chart figure of Complete games, Shutouts, and saves
    fig4 = go.Figure(
        go.Pie(values=PCT, labels=["Complete Games", "Shutouts", "Saves"], opacity=0.8)
    )
    fig4.update_layout(
        hovermode=False,
        title="% of Total Games Played",
        font={"color": "darkslategray"},
        paper_bgcolor="white",
        plot_bgcolor="#f8f5f0",
        legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1.15),
    )
    # Update figure trace marker colors
    fig4.update_traces(marker=dict(colors=["#462425", "#E35625", "#CEC6C0"]))

    # return figure
    return fig4


# Call back to Line Bubble Chart, take data request from dropdown menu
@app.callback(
    Output("pitch-bubble", "figure"),
    [Input("team-dropdown", "value"), Input("era-slider", "value")],
)
def update_figure5(selected_team, year_range):
    # Create filter dataframe of requested team
    filter_team = teams_df[teams_df.team_id == selected_team]
    # I will revisit this again soon, it just doesnt seem efficient
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]

    # Create lists of team data
    ERA = filter_year.era
    SOA = filter_year.soa
    BBA = filter_year.bba
    # Calculate K/BB ratio
    RATIO = SOA / BBA

    # Create line chart of K/BB ratio with ERA used for bubble size
    fig5 = go.Figure(
        data=go.Scatter(
            x=filter_year.year,
            y=RATIO,
            mode="markers",
            marker=dict(symbol="circle-open-dot", size=8.0 * ERA, color="#006BA6"),
            hovertemplate="K/BB: %{y:.2f}<extra></extra><br>" + "%{text}",
            text=["ERA: {}".format(i) for i in ERA],
        )
    )

    fig5.update_xaxes(title="Year", tickformat="d")
    fig5.update_yaxes(title="K/BB Ratio")
    fig5.update_layout(
        hovermode="x",
        title="K/BB Ratio with ERA Bubble",
        font={"color": "darkslategray"},
        paper_bgcolor="white",
        plot_bgcolor="#f8f5f0",
    )

    return fig5


# Callback to player dropdown menu
@app.callback(
    [Output("player-dropdown", "options"), Output("player-dropdown", "value")],
    [Input("team-dropdown", "value"), Input("era-slider", "value")],
)
def update_player_dropdown(selected_team, year_range):
    # Set filter dataframe of selected team
    filter_team = team_players_df[team_players_df.team_id == selected_team]

    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]

    # Use known names to set a key value pair for dropdown
    names = [
        {"label": k, "value": v}
        for k, v in zip(filter_year.known_name, filter_year.player_id)
    ]
    # create empty list
    player_list = []
    # append to list if player is not in list
    [player_list.append(x) for x in names if x not in player_list]

    # Return given name and id key value pair to player dropdown
    return player_list, player_list[0]["value"]


# Callback to Player profile datatable
@app.callback(
    [Output("playerProfile", "data"), Output("playerProfile", "columns")],
    [Input("player-dropdown", "value")],
)
def update_profile_table(player):
    # Create player filter with selected player
    filter_player = player_df[player_df.player_id == player]

    # drop unneccesary columns
    data_filter = filter_player.drop(
        columns=[
            "player_id",
            "name_first",
            "name_last",
            "name_given",
            "retro_id",
            "bbref_id",
            "birth_month",
            "birth_day",
            "birth_country",
            "birth_city",
            "birth_state",
            "death_month",
            "death_day",
            "death_country",
            "death_city",
            "death_state",
            "final_game",
        ]
    )

    # Return player profile to datatable
    return data_filter.to_dict("records"), [{"name": x, "id": x} for x in data_filter]


# Callback to players batting datatable
@app.callback(
    [Output("batterTable", "data"), Output("batterTable", "columns")],
    [
        Input("player-dropdown", "value"),
        Input("team-dropdown", "value"),
        Input("era-slider", "value"),
    ],
)
def update_batter_table(player, selected_team, year_range):
    # take in the selected team
    filter_team = batter_df[batter_df.team_id == selected_team]
    # Set and filter year range
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]

    # Apply filter to player ID in batters dataframe
    filter_batter = filter_year[filter_year.player_id == player]
    # drop unneccesary columns
    data_filter = filter_batter.drop(
        columns=["player_id", "team_id", "stint", "league_id"]
    )

    # Return batters data and batters key value pair to columns
    return data_filter.to_dict("records"), [{"name": x, "id": x} for x in data_filter]


# Call back to OBP Line Graph
@app.callback(
    Output("obp-line", "figure"),
    [
        Input("player-dropdown", "value"),
        Input("team-dropdown", "value"),
        Input("era-slider", "value"),
    ],
)
def update_figure6(player, selected_team, year_range):
    # Create filter dataframe of requested team data
    filter_team = batter_df[batter_df.team_id == selected_team]
    # Filter year range
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        year_range = [1903, 1919]
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    # Filter player
    filter_batter = filter_year[filter_year.player_id == player]

    # Make calculations for graph
    OBP = data.calculate_obp(filter_batter)
    # WOBA = data.calculate_woba(filter_batter)

    # Create Line char figure
    fig6 = go.Figure()
    # Add OBP line
    fig6.add_trace(
        go.Scatter(
            name="On-Base Percentage",
            x=filter_batter.year,
            y=OBP * 750,
            mode="lines+markers",
            marker_color="Orange",
            hovertemplate="OBP: %{text:.3f}<extra></extra><br>",
            text=["{}".format(i) for i in OBP],
        )
    )
    # WOBA line
    # fig6.add_trace(go.Scatter(name='Weighted On-Base Average', x=filter_batter.year, y=WOBA*750, mode='lines+markers', marker_color='orangered',
    #     hovertemplate = 'WOBA: %{text:.3f}<extra></extra><br>',
    #     text = ['{}'.format(i) for i in WOBA]))

    # add supporting bar charts, stats that are required for calculating OBP
    # exception of at-bats which would dwarf the other stats... not that hits dont :\
    fig6.add_trace(
        go.Bar(
            name="Hits",
            x=filter_batter.year,
            y=filter_batter.h,
            marker_color="saddlebrown",
            opacity=0.4,
        )
    )
    fig6.add_trace(
        go.Bar(
            name="Walks",
            x=filter_batter.year,
            y=filter_batter.bb,
            marker_color="tan",
            opacity=0.4,
        )
    )
    fig6.add_trace(
        go.Bar(
            name="Hit-By-Pitch",
            x=filter_batter.year,
            y=filter_batter.hbp,
            marker_color="Orange",
            opacity=0.4,
        )
    )
    fig6.add_trace(
        go.Bar(
            name="Sacrifice Flys",
            x=filter_batter.year,
            y=filter_batter.sf,
            marker_color="black",
            opacity=0.4,
        )
    )

    # Update figure
    fig6.update_xaxes(title="year", tickformat="d")
    fig6.update_yaxes(fixedrange=True)
    fig6.update_layout(
        barmode="stack",
        hovermode="x",
        title="On-Base Performance",
        font={"color": "darkslategray"},
        paper_bgcolor="white",
        plot_bgcolor="#f8f5f0",
        legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1.34),
    )

    return fig6


# Call back to SLG Line Graph
@app.callback(
    Output("slg-line", "figure"),
    [
        Input("player-dropdown", "value"),
        Input("team-dropdown", "value"),
        Input("era-slider", "value"),
    ],
)
def update_figure7(player, selected_team, year_range):
    # Create filter dataframe of requested team data
    filter_team = batter_df[batter_df.team_id == selected_team]
    # Set year range
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        year_range = [1903, 1919]
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    # Filter player
    filter_batter = filter_year[filter_year.player_id == player]

    # Set singles
    SNG = (
        filter_batter.h - filter_batter.double - filter_batter.triple - filter_batter.hr
    )

    # Calculate Slugging Average
    SLG = data.calculate_slg(filter_batter)

    # Create figure
    fig7 = go.Figure()
    # add SLG scatter plot
    fig7.add_trace(
        go.Scatter(
            name="Slugging Average",
            x=filter_batter.year,
            y=SLG * 325,
            mode="lines+markers",
            marker_color="dodgerblue",
            hovertemplate="SLG: %{text:.3f}<extra></extra><br>",
            text=["{}".format(i) for i in SLG],
        )
    )
    # add supporting bar charts, stats that are required for calculating SLG
    fig7.add_trace(
        go.Bar(
            name="Singles",
            x=filter_batter.year,
            y=SNG,
            marker_color="#005A9C",
            opacity=0.4,
        )
    )
    fig7.add_trace(
        go.Bar(
            name="Doubles",
            x=filter_batter.year,
            y=filter_batter.double,
            marker_color="#A5ACAF",
            opacity=0.4,
        )
    )
    fig7.add_trace(
        go.Bar(
            name="Triples",
            x=filter_batter.year,
            y=filter_batter.triple,
            marker_color="#EF3E42",
            opacity=0.4,
        )
    )
    fig7.add_trace(
        go.Bar(
            name="Home Runs",
            x=filter_batter.year,
            y=filter_batter.hr,
            marker_color="darkslategray",
            opacity=0.4,
        )
    )

    # Update figure
    fig7.update_xaxes(title="Year", tickformat="d")
    fig7.update_yaxes(fixedrange=True)
    fig7.update_layout(
        barmode="stack",
        hovermode="x",
        title="Slugging Performance",
        font={"color": "darkslategray"},
        paper_bgcolor="white",
        plot_bgcolor="#f8f5f0",
        legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1.15),
    )

    return fig7


# Call back to OPS Line Chart, Takes request data from dropdown menu
@app.callback(
    Output("ops-line", "figure"),
    [
        Input("player-dropdown", "value"),
        Input("team-dropdown", "value"),
        Input("era-slider", "value"),
    ],
)
def update_figure8(player, selected_team, year_range):
    # Create filter dataframe of requested team data
    filter_team = batter_df[batter_df.team_id == selected_team]
    # Set year range
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        year_range = [1903, 1919]
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    # Filter player
    filter_batter = filter_year[filter_year.player_id == player]

    # Calculate On-Base Percentage
    OBP = data.calculate_obp(filter_batter)
    # Calculate Slugging Average
    SLG = data.calculate_slg(filter_batter)
    # Calculate On-Base + Slugging
    OPS = OBP + SLG

    # Create line chart figure using calculations
    fig8 = go.Figure(
        data=[
            go.Scatter(
                name="OPS",
                x=filter_batter.year,
                y=OPS,
                mode="lines+markers",
                marker_color="green",
                hovertemplate="OPS: %{text:.3f}<extra></extra><br>",
                text=["{}".format(i) for i in OPS],
            ),
            go.Scatter(
                name="SLG",
                x=filter_batter.year,
                y=SLG,
                mode="lines+markers",
                marker_color="dodgerblue",
                hovertemplate="SLG: %{text:.3f}<extra></extra><br>",
                text=["{}".format(i) for i in SLG],
            ),
            go.Scatter(
                name="OBP",
                x=filter_batter.year,
                y=OBP,
                mode="lines+markers",
                marker_color="orange",
                hovertemplate="OBP: %{text:.3f}<extra></extra><br>",
                text=["{}".format(i) for i in OBP],
            ),
        ]
    )

    # Update figure
    fig8.update_xaxes(title="Year", tickformat="d")
    fig8.update_yaxes(fixedrange=True)
    fig8.update_layout(
        hovermode="x",
        title="On-Base + Slugging (OPS) Performance",
        font={"color": "darkslategray"},
        paper_bgcolor="white",
        plot_bgcolor="#f8f5f0",
        legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1.27),
    )

    return fig8


# fielding datatable
@app.callback(
    [Output("fieldTable", "data"), Output("fieldTable", "columns")],
    [
        Input("player-dropdown", "value"),
        Input("team-dropdown", "value"),
        Input("era-slider", "value"),
    ],
)
def update_field_table(player, selected_team, year_range):
    # Filter the selected team
    filter_team = field_df[field_df.team_id == selected_team]
    # Filter year range
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]

    # Apply filter player id
    filter_player = filter_year[filter_year.player_id == player]

    # drop unneccesary columns
    data_filter = filter_player.drop(
        columns=[
            "player_id",
            "team_id",
            "stint",
            "league_id",
            "inn_outs",
            "gs",
            "pb",
            "wp",
            "sb",
            "cs",
            "zr",
        ]
    )

    # Return batters dictionary to data and batters key value pair to columns
    return data_filter.to_dict("records"), [{"name": x, "id": x} for x in data_filter]


# pitchers datatable
@app.callback(
    [Output("pitch-data", "children")],
    [
        Input("player-dropdown", "value"),
        Input("team-dropdown", "value"),
        Input("era-slider", "value"),
    ],
)
def update_pitch_table(player, selected_team, year_range):
    # take in the selected team
    filter_team = pitching_df[pitching_df.team_id == selected_team]
    # Filter year range
    if year_range:
        filter_year = filter_team[
            (filter_team.year >= year_range[0]) & (filter_team.year <= year_range[1])
        ]
    else:
        filter_year = filter_team[
            (filter_team.year >= 1903) & (filter_team.year <= 1919)
        ]
    # Filter player id
    filter_batter = filter_year[filter_year.player_id == player]

    # drop unneccesary columns
    data_filter = filter_batter.drop(
        columns=["player_id", "team_id", "stint", "league_id", "ipouts", "baopp"]
    )

    # Set empty list
    data_note = []
    # if data filter is empty, append and return notice
    if data_filter.empty:
        data_note.append(
            html.Div(
                dbc.Alert(
                    "No pitching data is available for selected player.",
                    color="warning",
                )
            )
        )
        return data_note
    # else set and return datatable
    else:
        data_note.append(
            html.Div(
                dash_table.DataTable(
                    data=data_filter.to_dict("records"),
                    columns=[{"name": x, "id": x} for x in data_filter],
                    style_as_list_view=True,
                    editable=False,
                    style_table={
                        "overflowY": "scroll",
                        "width": "100%",
                        "minWidth": "100%",
                    },
                    style_header={"backgroundColor": "#f8f5f0", "fontWeight": "bold"},
                    style_cell={"textAlign": "center", "padding": "8px"},
                )
            )
        )
        return data_note


# Callback to position dropdown
@app.callback(
    [Output("pos-dropdown", "options"), Output("pos-dropdown", "value")],
    [
        Input("player-dropdown", "value"),
        Input("team-dropdown", "value"),
        Input("era-slider", "value"),
    ],
)
def update_pos_dropdown(player, selected_team, year_range):
    # Filter team
    team = field_df[field_df.team_id == selected_team]
    # Filter year range
    if year_range:
        filter_year = team[(team.year >= year_range[0]) & (team.year <= year_range[1])]
    else:
        filter_year = team[(team.year >= 1903) & (team.year <= 1919)]
    # filter player id
    filter_player = filter_year[filter_year.player_id == player]

    # Use position to set a key value pair for dropdown
    positions = [
        {"label": k, "value": v} for k, v in zip(filter_player.pos, filter_player.pos)
    ]
    # set empty list
    pos_list = []
    # append list if position is not in list
    [pos_list.append(x) for x in positions if x not in pos_list]

    # Return given name key value pair to options and value of dropdown
    return pos_list, pos_list[0]["value"]


# Call back to Fielding/Pitching Graphs
@app.callback(
    Output("pitch-graphs", "children"),
    [
        Input("player-dropdown", "value"),
        Input("pos-dropdown", "value"),
        Input("team-dropdown", "value"),
        Input("era-slider", "value"),
    ],
)
def update_figure9(player, position, selected_team, year_range):
    # Create filter dataframe of requested team data
    field_filter = field_df[field_df.team_id == selected_team]
    pitch_filter = pitching_df[pitching_df.team_id == selected_team]

    # Filter year range
    if year_range:
        filter_field_year = field_filter[
            (field_filter.year >= year_range[0]) & (field_filter.year <= year_range[1])
        ]
        filter_pitch_year = pitch_filter[
            (pitch_filter.year >= year_range[0]) & (pitch_filter.year <= year_range[1])
        ]
    else:
        filter_field_year = field_filter[
            (field_filter.year >= 1903) & (field_filter.year <= 1919)
        ]
        filter_pitch_year = pitch_filter[
            (pitch_filter.year >= 1903) & (pitch_filter.year <= 1919)
        ]

    # Set empty list
    data_note = []
    # if position is NOT pitcher create only a fielding percentage chart
    if position != "P":
        data_note.append(
            html.Div(
                dbc.Alert(
                    "No pitching graphs are available for selected player.",
                    color="warning",
                ),
                style={"padding-top": "1%"},
            )
        )
        return data_note
    # else if position is pitcher then create subplot with feilding and pitching performances
    elif position == "P":
        # Filter player
        filter_pitcher = filter_pitch_year[filter_pitch_year.player_id == player]
        filter_field_player = filter_field_year[filter_field_year.player_id == player]
        filter_pos = filter_field_player[filter_field_player.pos == position]

        # Create subplot, spacing and titles
        fig9 = make_subplots(
            rows=3,
            cols=1,
            subplot_titles=(
                "Walks + Hits over Innings Pitched",
                "Strikeout Rate with ERA Bubble",
                "Win-Loss Performance",
            ),
        )

        ###### WHIP #####
        # Set innings and calculate WHIP
        IP = filter_pitcher.ipouts / 3
        WHIP = (filter_pitcher.bb + filter_pitcher.h) / IP
        # set whip as line plot
        whip_fig = go.Scatter(
            name="WHIP",
            x=filter_pitcher.year,
            y=WHIP,
            legendgroup="group4",
            mode="lines+markers",
            marker_color="lightseagreen",
            hovertemplate="WHIP: %{text:.3f}<extra></extra><br>",
            text=["{}".format(i) for i in WHIP],
        )
        # apend figure to subplot 1,1
        fig9.append_trace(whip_fig, row=1, col=1)

        ##### K/BB ratio with ERA bubble #####
        # Create ERA variable and K/BB ratio calculation
        ERA = filter_pitcher.era
        RATIO = filter_pitcher.so / filter_pitcher.bb
        # Create line chart of K/BB ratio with ERA used for bubble size
        era_bubble_fig = go.Scatter(
            name="K/BB Ratio + ERA",
            x=filter_pitcher.year,
            y=RATIO,
            legendgroup="group3",
            mode="markers",
            marker=dict(
                symbol="circle-open-dot",
                size=ERA,
                sizemode="area",
                sizeref=2.0 * max(ERA) / (40.0 ** 2),
                sizemin=4,
                color="darkblue",
            ),
            hovertemplate="K/BB: %{y:.2f}<extra></extra><br>" + "%{text}",
            text=["ERA: {}".format(i) for i in ERA],
        )
        # append figure to subplot 2,1
        fig9.append_trace(era_bubble_fig, row=2, col=1)

        ##### Win Percentage #####
        # calculate winning percentage
        wl_pct = filter_pitcher.w / (filter_pitcher.w + filter_pitcher.l)
        # set Win Loss Performance
        pwin_bar = go.Bar(
            name="wPCT",
            x=filter_pitcher.year,
            y=wl_pct,
            legendgroup="group2",
            marker_color="#00A3E0",
            opacity=0.4,
            hovertemplate="wPCT: %{y:.3f}<extra></extra><br>",
        )
        # add figure to subplot 3,1
        fig9.append_trace(pwin_bar, row=3, col=1)

        # Update layout, set hover to x-axis, barmode, establish title and colors
        fig9.update_layout(
            height=1400,
            hovermode="x",
            barmode="stack",
            xaxis=dict(title="Year", tickformat="d"),
            yaxis=dict(title="WHIP"),
            xaxis2=dict(title="Year", tickformat="d"),
            yaxis2=dict(title="K/BB Ratio"),
            xaxis3=dict(title="Year", tickformat="d"),
            yaxis3=dict(title="W-L PCT"),
            font={"color": "darkslategray"},
            paper_bgcolor="white",
            plot_bgcolor="#f8f5f0",
            legend=dict(
                orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1
            ),
        )
        # return subplots
        graph = dcc.Graph(figure=fig9, config={"displayModeBar": False})
        return graph


# Call back to Fielding/Pitching Graphs
@app.callback(
    Output("field-graph", "figure"),
    [
        Input("player-dropdown", "value"),
        Input("pos-dropdown", "value"),
        Input("team-dropdown", "value"),
        Input("era-slider", "value"),
    ],
)
def update_figure9(player, position, selected_team, year_range):
    # Create filter dataframe of requested team data
    field_filter = field_df[field_df.team_id == selected_team]
    pitch_filter = pitching_df[pitching_df.team_id == selected_team]

    # Filter year range
    if year_range:
        filter_field_year = field_filter[
            (field_filter.year >= year_range[0]) & (field_filter.year <= year_range[1])
        ]
        filter_pitch_year = pitch_filter[
            (pitch_filter.year >= year_range[0]) & (pitch_filter.year <= year_range[1])
        ]
    else:
        filter_field_year = field_filter[
            (field_filter.year >= 1903) & (field_filter.year <= 1919)
        ]
        filter_pitch_year = pitch_filter[
            (pitch_filter.year >= 1903) & (pitch_filter.year <= 1919)
        ]

    # Filter player
    filter_field_player = filter_field_year[filter_field_year.player_id == player]
    # filter position
    filter_pos = filter_field_player[filter_field_player.pos == position]

    # Calculate feilding Average
    FLDP = (filter_pos.po + filter_pos.a) / (
        filter_pos.po + filter_pos.a + filter_pos.e
    )

    # Set figure
    fig10 = go.Figure()
    # add Fielding percentage line plot
    fig10.add_trace(
        go.Scatter(
            name="Fielding Percentage",
            x=filter_pos.year,
            y=FLDP * 1000,
            mode="lines+markers",
            marker_color="rebeccapurple",
            hovertemplate="FLD: %{text:.3f}<extra></extra><br>",
            text=["{}".format(i) for i in FLDP],
        )
    )
    # add supporting bar charts, stats that are required for calculating Fielding percentage
    fig10.add_trace(
        go.Bar(
            name="Put Outs",
            x=filter_pos.year,
            y=filter_pos.po,
            marker_color="#33006F",
            opacity=0.4,
        )
    )
    fig10.add_trace(
        go.Bar(
            name="Assists",
            x=filter_pos.year,
            y=filter_pos.a,
            marker_color="#C4CED4",
            opacity=0.4,
        )
    )
    fig10.add_trace(
        go.Bar(
            name="Errors",
            x=filter_pos.year,
            y=filter_pos.e,
            marker_color="black",
            opacity=0.4,
        )
    )
    # update figure
    fig10.update_xaxes(title="Year", tickformat="d")
    fig10.update_yaxes(fixedrange=True)
    fig10.update_layout(
        height=400,
        barmode="stack",
        hovermode="x",
        title="Fielding Performance",
        font={"color": "darkslategray"},
        paper_bgcolor="white",
        plot_bgcolor="#f8f5f0",
        legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1.27),
    )
    # Return figure
    return fig10
