import numpy as np
import pandas as pd
from sklearn import datasets
from sklearn.model_selection import train_test_split

test_size = 0.25


def sampling(**kwargs):
    if kwargs['dataset'] == 'moons':
        X, y = datasets.make_moons(n_samples=kwargs['sample_size'],
                                   noise=kwargs['noise'],
                                   random_state=5)

        return train_test_split(X,
                                y.astype(str),
                                test_size=kwargs['test_size'],
                                random_state=5), X, y

    elif kwargs['dataset'] == 'circles':
        X, y = datasets.make_circles(n_samples=kwargs['sample_size'],
                                     noise=kwargs['noise'],
                                     factor=0.5,
                                     random_state=1)
        return train_test_split(X,
                                y.astype(str),
                                test_size=kwargs['test_size'],
                                random_state=5), X, y

    elif kwargs['dataset'] == 'LS':
        X, y = datasets.make_classification(n_samples=kwargs['sample_size'],
                                            n_features=2,
                                            n_redundant=0,
                                            n_informative=2,
                                            random_state=2,
                                            n_clusters_per_class=1)

        rng = np.random.RandomState(2)
        X += kwargs['noise'] * rng.uniform(size=X.shape)

        return train_test_split(X,
                                y.astype(str),
                                test_size=kwargs['test_size'],
                                random_state=5), X, y

    else:
        return ValueError('error!')


def df_split(**kwargs):
    _df = kwargs['df']

    return train_test_split(
        _df[['x', 'y']].to_numpy(),
        _df['c'].to_numpy().astype(str),
        test_size=kwargs['test_size'],
        random_state=5), _df[['x', 'y']].to_numpy(), _df['c'].to_numpy()


def data_split(**kwargs):

    return train_test_split(kwargs['X'],
                            kwargs['y'].astype(str),
                            test_size=kwargs['test_size'],
                            random_state=5), kwargs['X'], kwargs['y']
