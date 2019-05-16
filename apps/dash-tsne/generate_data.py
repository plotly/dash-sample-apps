# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from keras.datasets import *
import os
import sys

# create the data folder if it doesn't exist yet
if not os.path.exists("data"):
    os.makedirs("data")

# Try to parse the dataset argument
dataset_name = sys.argv[1]

# Try to parse the sample size argument
try:
    sample_size = int(sys.argv[2])
except ValueError as e:
    print("The input value for sample size is invalid:", e)
    sys.exit(1)

# Dataset indices
mnist_idx = [
    "Digit 0",
    "Digit 1",
    "Digit 2",
    "Digit 3",
    "Digit 4",
    "Digit 5",
    "Digit 6",
    "Digit 7",
    "Digit 8",
    "Digit 9",
]

cifar_idx = [
    "airplane",
    "automobile",
    "bird",
    "cat",
    "deer",
    "dog",
    "frog",
    "horse",
    "ship",
    "truck",
]

fashion_idx = [
    "T-Shirt",
    "Trouser",
    "Pullover",
    "Dress",
    "Coat",
    "Sandal",
    "Shirt",
    "Sneaker",
    "Bag",
    "Ankle boot",
]

# Load Dataset
if dataset_name.lower() in ["mnist", "mnist_digits", "mnistdigits"]:
    X, y = mnist.load_data()[0]
    selected_idx = mnist_idx

elif dataset_name.lower() in ["cifar", "cifar10"]:
    X, y = cifar10.load_data()[0]
    selected_idx = cifar_idx

elif dataset_name.lower() in ["cifar_gray", "cifar10_gray", "cifargray"]:
    from skimage.color import rgb2gray

    X, y = cifar10.load_data()[0]
    selected_idx = cifar_idx
    X = np.moveaxis(X, 1, 3)
    X = rgb2gray(X) * 255
    print("Array Shape:", X.shape)

elif dataset_name.lower() in ["fashion", "fashion_mnist", "fashionmnist"]:
    X, y = fashion_mnist.load_data()[0]
    selected_idx = fashion_idx

else:
    print("Dataset not found.")
    sys.exit(1)

y = np.array([selected_idx[int(val)] for val in y])

print("Dataset loaded.")

# Flatten the array, and normalize it
X = X.reshape(X.shape[0], -1) / 255.0

# We will select the integer values to be the index
df = pd.DataFrame(X, index=y)

if sample_size > df.shape[0]:
    print("Sample size is too great.")
    sys.exit(1)

samples = df.sample(n=sample_size, random_state=1234)

samples.to_csv(f"data/{dataset_name}_{sample_size}_input.csv", index=False)
pd.DataFrame(samples.index).to_csv(
    f"data/{dataset_name}_{sample_size}_labels.csv", index=False
)

print("CSV files created.")
