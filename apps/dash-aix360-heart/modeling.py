from aix360.algorithms.rbm import LogisticRuleRegression, FeatureBinarizer
import pandas as pd
from sklearn.model_selection import train_test_split

col_map = {
    "trestbps": "Resting blood pressure",
    "chol": "Cerum cholestoral",
    "fbs": "asting blood sugar",
    "restecg": "resting electrocardiographic results",
    "thalach": "maximum heart rate achieved",
    "exang": "exercise induced angina",
    "oldpeak": "ST depression due to exercise vs. rest",
    "age": "age",
    "sex": "sex",
    "cp": "chest pain type",
    "slope": "slope of peak exercise ST segment",
    "ca": "num. major vessels",
    "thal": "defect type",
}

# Load and preprocess dataset
df = pd.read_csv("heart.csv")
y = df.pop("target")
dfTrain, dfTest, yTrain, yTest = train_test_split(df, y, random_state=0, stratify=y)

fb = FeatureBinarizer(negations=True, returnOrd=True)
dfTrain, dfTrainStd = fb.fit_transform(dfTrain)
dfTest, dfTestStd = fb.transform(dfTest)

# Train model
lrr = LogisticRuleRegression(lambda0=0.005, lambda1=0.001, useOrd=True)
lrr.fit(dfTrain, yTrain, dfTrainStd)
