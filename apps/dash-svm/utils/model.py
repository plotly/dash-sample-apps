from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
import numpy as np

import utils.figures as figs
from utils.helper_functions import generate_data

def svm_prediction(kernel, degree, C_coef, C_power, gamma_coef, gamma_power, dataset, noise, shrinking, threshold, sample_size):
    h = 0.3  # step size in the mesh

    # Data Pre-processing
    X, y = generate_data(n_samples=sample_size, dataset=dataset, noise=noise)
    X = StandardScaler().fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.4, random_state=42
    )

    x_min = X[:, 0].min() - 0.5
    x_max = X[:, 0].max() + 0.5
    y_min = X[:, 1].min() - 0.5
    y_max = X[:, 1].max() + 0.5
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

    C = C_coef * 10 ** C_power
    gamma = gamma_coef * 10 ** gamma_power

    if shrinking == "True":
        flag = True
    else:
        flag = False

    # Train SVM
    clf = SVC(C=C, kernel=kernel, degree=degree, gamma=gamma, shrinking=flag)
    clf.fit(X_train, y_train)

    # Plot the decision boundary. For that, we will assign a color to each
    # point in the mesh [x_min, x_max]x[y_min, y_max].
    if hasattr(clf, "decision_function"):
        Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
    else:
        Z = clf.predict_proba(np.c_[xx.ravel(), yy.ravel()])[:, 1]

    prediction_figure = figs.serve_prediction_plot(
        model=clf,
        X_train=X_train,
        X_test=X_test,
        y_train=y_train,
        y_test=y_test,
        Z=Z,
        xx=xx,
        yy=yy,
        mesh_step=h,
        threshold=threshold,
    )

    roc_figure = figs.serve_roc_curve(model=clf, X_test=X_test, y_test=y_test)

    confusion_figure = figs.serve_pie_confusion_matrix(
        model=clf, X_test=X_test, y_test=y_test, Z=Z, threshold=threshold
    )
    return prediction_figure, roc_figure, confusion_figure
