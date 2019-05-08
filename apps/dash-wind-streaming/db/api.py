import sqlite3
import pandas as pd


def get_wind_data(start, end):
    """
    Query wind data rows between two ranges

    :params start: start row id
    :params end: end row id
    :returns: pandas dataframe object 
    """

    con = sqlite3.connect("./db/wind-data.db")
    statement = f'SELECT Speed, SpeedError, Direction FROM Wind WHERE rowid > "{start}" AND rowid <= "{end}";'
    df = pd.read_sql_query(statement, con)
    return df


def get_wind_data_by_id(id):
    """
    Query a row from the Wind Table

    :params id: a row id
    :returns: pandas dataframe object 
    """

    con = sqlite3.connect("./db/wind-data.db")
    statement = f'SELECT * FROM Wind WHERE rowid = "{id}";'
    df = pd.read_sql_query(statement, con)
    return df
