import pandas as pd

import utils

df = pd.read_csv("clean_loan.csv", nrows=5)
print(utils.get_column_strings(df))
