import os

from sqlalchemy import create_engine
import pandas as pd

import utils

SYNAPSE_HOST = os.environ["SYNAPSE_HOST"]
SYNAPSE_USER = os.environ["SYNAPSE_USER"]
SYNAPSE_PW = os.environ["SYNAPSE_PW"]

driver = "ODBC+Driver+17+for+SQL+Server"
db_name = "loan"
port = "1433"

engine = create_engine(
    f"mssql+pyodbc://{SYNAPSE_USER}:{SYNAPSE_PW}@{SYNAPSE_HOST}:{port}/{db_name}?driver={driver}"
)


query = """
SELECT TOP 100 *
FROM loan.dbo.cleanLoan
"""
connection = engine.connect()
result = pd.read_sql(query, connection)
connection.close()


connection = engine.connect()

loan_min, loan_max = utils.get_range(connection, db_name, "cleanLoan", "loan_amnt")
inc_min, inc_max = utils.get_range(connection, db_name, "cleanLoan", "annual_inc")
app_types = utils.get_unique(connection, db_name, "cleanLoan", "application_type")
purposes = utils.get_unique(connection, db_name, "cleanLoan", "purpose")
ownerships = utils.get_unique(connection, db_name, "cleanLoan", "home_ownership")[1:]

connection.close()
