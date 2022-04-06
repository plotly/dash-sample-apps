import numpy as np
from sklearn.svm import SVC

margin = 0.25
mesh_size = 0.02


def modeling(**kwargs):

    _data = kwargs['data']
    split_data = _data[0]

    

    x_min, x_max = _data[1][:, 0].min() - margin, _data[1][:, 0].max() + margin
    y_min, y_max = _data[1][:, 1].min() - margin, _data[1][:, 1].max() + margin

    xrange = np.arange(x_min, x_max, mesh_size)
    yrange = np.arange(y_min, y_max, mesh_size)

    xx, yy = np.meshgrid(xrange, yrange)

    clf = SVC(C=kwargs['cost'],
              kernel=kwargs['kernel'],
              degree=kwargs['degree'],
              gamma=kwargs['gamma'],
              shrinking=kwargs['shrinking'])

    clf.fit(split_data[0], split_data[2])

    if hasattr(clf, "decision_function"):
        Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
    else:
        Z = clf.predict_proba(np.c_[xx.ravel(), yy.ravel()])[:, 1]

    return clf, Z, xx, yy, xrange, yrange
