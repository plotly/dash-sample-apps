import os
from textwrap import dedent

from snowflake.sqlalchemy import URL
from sqlalchemy import create_engine
import pandas as pd


def create_or_replace_table(connection, df, name):
    # Load the actual csv file

    # Create SQL columns based on the columns of that dataframe
    types = (
        df.dtypes.copy()
        .replace("float64", "FLOAT")
        .replace("int64", "INT")
        .replace("object", "STRING")
    )

    sql_cols = ",\n".join(types.index + " " + types.values)

    cmd = f"""
CREATE OR REPLACE TABLE {name} (

{sql_cols}
    
)

STAGE_FILE_FORMAT = ( TYPE = 'csv' FIELD_DELIMITER= ',');
    """

    return connection.execute(cmd)


# Snowflake vars
FLAKE_ACCOUNT = os.getenv("FLAKE_ACCOUNT")
FLAKE_USER = os.getenv("FLAKE_USER")
FLAKE_PW = os.getenv("FLAKE_PW")

flake_warehouse = "snowflake_demos"
flake_db = "NEW_CLIENTS"
stage = "NEW_CLIENT_STAGE"

# New table vars
path = "owner.csv"
table_name = "OWNER"
abs_path = os.path.abspath(path)

# Create Engine and connect to DB
engine = create_engine(
    URL(
        account=FLAKE_ACCOUNT,
        user=FLAKE_USER,
        password=FLAKE_PW,
        database=flake_db,
        schema="public",
        warehouse=flake_warehouse,
        role="sysadmin",
    ),
    pool_size=5,
    pool_recycle=1800,
    pool_pre_ping=True,
)

connection = engine.connect()


# Create Table
df = pd.read_csv(abs_path, nrows=5)
create_or_replace_table(connection, df, name=table_name)

# First remove all the content of the stage
connection.execute(f"REMOVE @{stage}")

# Upload the CSV file
connection.execute(f"PUT file://{abs_path} @{stage}")

# Copy the CSV file into the table
connection.execute(
    f"COPY INTO {table_name} FROM @{stage} FILE_FORMAT=(TYPE=CSV SKIP_HEADER=1)"
)
