from aix360.algorithms.rbm import LogisticRuleRegression, FeatureBinarizer
import pandas as pd
from sklearn.model_selection import train_test_split

col_map = {
    "trestbps": "Resting blood pressure (trestbps)",
    "chol": "Cerum cholestoral (chol)",
    "fbs": "Fasting blood sugar (fbs)",
    "restecg": "Resting electrocardiographic results (restecg)",
    "thalach": "Maximum heart rate achieved (thalach)",
    "exang": "Exercise induced angina (exang)",
    "oldpeak": "S-T depression induced by exercise relative to rest (oldpeak)",
    "age": "Age",
    "sex": "Sex",
    "cp": "Chest pain type (cp)",
    "slope": "Slope of peak exercise S-T segment (slope)",
    "ca": "Number of major vessels (ca)",
    "thal": "Defect type (thal)",
}

num2desc = {
    "sex": {0: "Female", 1: "Male"},
    "cp": {
        0: "typical angina",
        1: "atypical angina",
        2: "non-aginal pain",
        3: "asymptomatic",
    },
    "fbs": {0: "False", 1: "True"},
    "restecg": {
        0: "normal",
        1: "ST-T wave abnormality",
        2: "left ventricular hypertrophy",
    },
}

# Load and preprocess dataset
df = pd.read_csv("heart.csv")
for k, v in num2desc.items():
    df[k] = df[k].replace(v)

y = df.pop("target")
dfTrain, dfTest, yTrain, yTest = train_test_split(df, y, random_state=0, stratify=y)

fb = FeatureBinarizer(negations=True, returnOrd=True)
dfTrain, dfTrainStd = fb.fit_transform(dfTrain)
dfTest, dfTestStd = fb.transform(dfTest)

# Train model
lrr = LogisticRuleRegression(lambda0=0.005, lambda1=0.001, useOrd=True)
lrr.fit(dfTrain, yTrain, dfTrainStd)
