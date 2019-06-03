# In[]:
# Import required libraries
import pandas as pd
import pickle


# In[]:
# Load required dataframes
dataset = dict()

columns1 = ['API Well Number', 'Gas Produced, MCF',
            'Water Produced, bbl', 'Reporting Year']
columns2 = ['API Well Number', 'Gas Produced, MCF',
            'Water Produced, bbl', 'Oil Produced, bbl', 'Reporting Year']
df1 = pd.read_csv('data/Oil_and_Gas_Annual_Production__1985_-_2000.csv',
                  low_memory=False)[columns1]
df2 = pd.read_csv('data/Oil_and_Gas_Annual_Production__Beginning_2001.csv',
                  low_memory=False)[columns2]


# In[]:
# Concatenate dataframes
df = pd.concat([df1, df2])
df.fillna(0, inplace=True)


# In[]:
# Append production data to each API
columns = ['Gas Produced, MCF', 'Water Produced, bbl',
           'Oil Produced, bbl', 'Reporting Year']
for api, df_well in df.groupby('API Well Number'):
    # print(api)
    # print(df_well)
    df_well = df_well[columns]
    df_well.index = df_well['Reporting Year']
    df_well = df_well.to_dict(orient='index')
    dataset[api] = df_well

# print(dataset)

# In[]:
# Serialise to Pickle
with open('data/points.pkl', 'wb') as f:
    pickle.dump(dataset, f)
