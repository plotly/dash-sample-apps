import pickle

import numpy as np
from cuml.manifold import TSNE
from sklearn.linear_model import LinearRegression

from helpers import load_mnist


train_images, train_labels = load_mnist("./fashion", subset="train")
test_images, test_labels = load_mnist("./fashion", subset="test")

all_images = np.concatenate((train_images, test_images))

tsne = TSNE(
    n_components=2,
    method="barnes_hut",
    random_state=23,
    learning_rate=200,
    perplexity=50,
    n_iter=3000,
)

train_X_hat = tsne.fit_transform(train_images)
test_X_hat = tsne.fit_transform(test_images)

all_X_hat = tsne.fit_transform(all_images)

np.save("trained_data/train_tsne", train_X_hat)
np.save("trained_data/test_tsne", test_X_hat)
np.save("trained_data/all_images_tsne", all_X_hat)

#################################
# infer approximate embeddings for new images:

all_X_hat = np.load("trained_data/all_images_tsne.npy")

reg = LinearRegression()
reg.fit(all_images, all_X_hat)

filename = "trained_data/linear_model_embeddings.sav"
pickle.dump(reg, open(filename, "wb"))
