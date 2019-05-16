import os
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import numpy as np
import pandas as pd

datasets = ["wikipedia_3000", "twitter_3000", "crawler_3000"]
iterations_ls = [250, 500, 750, 1000]
perplexity_ls = [3, 10, 30, 50, 100]
pca_dim_ls = [25, 50, 100]
learning_rate_ls = [10, 50, 100, 200]


def generate_embedding(
    dataset, iterations, perplexity, pca_dim, learning_rate, verbose=1, mode="two_files"
):
    path = f"demo_embeddings/{dataset}/iterations_{iterations}/perplexity_{perplexity}/pca_{pca_dim}/learning_rate_{learning_rate}"

    def display(string):
        if verbose:
            print(string)

    if os.path.exists(path):
        if os.path.exists(path + f"/data.csv"):
            display(path + " already exists.")
            return
    else:
        os.makedirs(path)

    if mode == "two_files":
        data = pd.read_csv(f"data/{dataset}_input.csv")
        labels = pd.read_csv(f"data/{dataset}_labels.csv")
    elif mode == "one_file":
        data = pd.read_csv(f"data/{dataset}.csv", index_col=0, encoding="ISO-8859-1")
        labels = data.index

    nb_col = data.shape[1]

    pca = PCA(n_components=min(nb_col, pca_dim))
    data_pca = pca.fit_transform(data.values)

    tsne = TSNE(
        n_components=3,
        n_iter=iterations,
        learning_rate=learning_rate,
        perplexity=perplexity,
        random_state=1131,
    )

    embedding = tsne.fit_transform(data_pca)

    embedding_df = pd.DataFrame(embedding, columns=["x", "y", "z"])

    embedding_df.index = np.squeeze(labels.values)

    embedding_df.to_csv(path + f"/data.csv")

    display(f"{path} has been generated.")


for dataset in datasets:
    for iterations in iterations_ls:
        for perplexity in perplexity_ls:
            for pca_dim in pca_dim_ls:
                for learning_rate in learning_rate_ls:
                    generate_embedding(
                        dataset,
                        iterations,
                        perplexity,
                        pca_dim,
                        learning_rate,
                        mode="one_file",
                    )
