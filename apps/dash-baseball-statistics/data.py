# Improtant note: This data file would ordinarily be used to connect with a proper database server
# more likely PostgreSQL, but thats me. I do plan on rewritting this in the future for such implementations.
# With that said, this file will be be very slow to run and only to demonstrate data processing using
# functions and pandas along with providing a central file for data references
#
# Import Pandas
import pandas as pd


# Import CSV data
# Import team historical statistics
# Some historical team names are correlated with their more modern counter part
# Custome CSV files where created from the original by combining data to allow
# for easier display of historical team data
teams = pd.read_csv("data/update_team.csv")
# Import Players batting data
batters = pd.read_csv("data/update_batting.csv")
# Import custom Fielding data
fielding = pd.read_csv("data/update_fielding.csv")
# Import custom pitching data
pitching = pd.read_csv("data/update_pitching.csv")
# Import Player profile data
players = pd.read_csv("data/update_player.csv")
# Import custom player and team id dataframe
team_players = pd.read_csv("data/player_team.csv")

# Hardcoded list of era names as key value pairs
era_list = [
    {"label": "Dead Ball ('03-'19)", "value": "Dead Ball"},
    {"label": "Live Ball ('20-'41)", "value": "Live Ball"},
    {"label": "Integration ('42-'60)", "value": "Integration"},
    {"label": "Expantion ('61-'76)", "value": "Expantion"},
    {"label": "Free Agency ('77-'93)", "value": "Free Agency"},
    {"label": "Steroid ('94-'05)", "value": "Steroid"},
    {"label": "Post-Steroid ('06-'15)", "value": "Post-Steroid"},
    {"label": "Statcast ('16-'20)", "value": "Statcast"},
]

# Era markers
era_marks = {
    1903: {"label": "1903"},
    1919: {"label": "1919"},
    1941: {"label": "1941"},
    1960: {"label": "1960"},
    1976: {"label": "1976"},
    1993: {"label": "1993"},
    2005: {"label": "2005"},
    2015: {"label": "2015"},
    2020: {"label": "2020"},
}


# Creates a dynamic list of team names based on era
def dynamicteams(x):
    # Hardcoded list of era time spans, wouldnt do it this way if the set where larger
    era_time = [
        (1903, 1919),
        (1920, 1941),
        (1942, 1960),
        (1961, 1976),
        (1977, 1993),
        (1994, 2005),
        (2006, 2015),
        (2016, 2020),
    ]
    # create a filter list of just years and team names
    filter_team_yr = teams[["year", "name", "team_id"]]
    # filter the above list by year span
    filter_year = filter_team_yr[
        (filter_team_yr.year >= era_time[x][0])
        & (filter_team_yr.year <= era_time[x][1])
    ]  # High Year
    # filter_year = filter_year[] # Low Year
    # Create a filter list of Team names based on years filtered
    filter_teams = filter_year["name"].unique()
    filter_team_ids = filter_year["team_id"].unique()
    # return unique list of team names as a list of key value pairs, rather than calling a function to create and return the list
    # list comp of key value pair
    # new is a list of names while x is the name in the list
    return [{"label": k, "value": v} for k, v in zip(filter_teams, filter_team_ids)]


def dynamicrange(x):
    # Hardcoded data is not typically what i do unless the set is small
    era_time = [
        (1903, 1919),
        (1920, 1941),
        (1942, 1960),
        (1961, 1976),
        (1977, 1993),
        (1994, 2005),
        (2006, 2015),
        (2016, 2020),
    ]
    return [era_time[x][0], era_time[x][1]]


# Calculate On-Base Percentage function
def calculate_obp(df):
    # Set lists of team data
    AB = df.ab
    Ht = df.h
    BB = df.bb
    HBP = df.hbp
    SF = df.sf
    # return On-Base Percentage
    return (Ht + BB + HBP) / (AB + BB + HBP + SF)


# Calculate Slugging Average
def calculate_slg(df):
    # Set lists of player data
    AB = df.ab
    Ht = df.h
    DBL = df.double
    TRP = df.triple
    HR = df.hr
    SNG = Ht - DBL - TRP - HR
    # return Slugging Average
    return (SNG + 2 * DBL + 3 * TRP + 4 * HR) / AB


# Calculate WOBA
def calculate_woba(df):
    # Selected players singles
    SNG = df.h - df.double - df.triple - df.hr
    # Weighted Plate Appearances
    WPA = df.ab + df.bb - df.ibb + df.sf + df.hbp
    # weighted on-base average, 2013 https://library.fangraphs.com/offense/woba/
    return (
        (0.690 * df.bb)
        + (0.722 * df.hbp)
        + (0.888 * SNG)
        + (1.271 * df.double)
        + (1.616 * df.triple)
        + (2.101 * df.hr)
    ) / WPA
    # weighted on-base average, 2019 https://en.wikipedia.org/wiki/WOBA#2019_Formula
    # return ((0.690 * df.bb) + (0.719 * df.hbp) + (0.87 * SNG) + (1.217 * df.double) + (1.529 * df.triple) + (1.94 * df.hr)) / WPA
