import pandas as pd
import numpy as np

# Load data
df = pd.read_csv("data/EKY_shale/EKY_shale_production/Well_ID-Table_1.csv")
df = df.applymap(lambda x: x.strip() if type(x) == str else x)

# Drop non-gas well
df = df[df["result"] == "GAS"]

# Drop non-identified well
ind = df[df["tdfm"] == "000"].index
df = df.drop(ind)

form_id = df["tdfm"].unique()

# Replace form id to form name
df_fm = pd.read_csv("data/EKY_shale/EKY_shale_production/FmCodes-Table_1.csv")
df_fm = df_fm.applymap(lambda x: x.strip() if type(x) == str else x)

df["fm_name"] = df["tdfm"]

for formation in form_id:
    formation_name = df_fm.loc[df_fm["fm_code"] == formation, "fm_name"].item().strip()
    df["fm_name"] = df["fm_name"].apply(
        lambda x: formation_name if x == formation else x
    )

# Generate random mineral composition for ternary plot
ternary_ax_title = ["Quartz", "Carbonate", "Clay"]
rand_composition = np.rint(np.random.dirichlet((10, 5, 3), len(df)) * 100).transpose()

df["Quartz"] = rand_composition[0]
df["Carbonate"] = rand_composition[1]
df["Clay"] = rand_composition[2]

df.to_csv("test_composition.csv")
