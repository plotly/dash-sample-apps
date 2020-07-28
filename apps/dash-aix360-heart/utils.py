import pandas as pd
import numpy as np


def compute_plot_gam(model, Xorig, fb, features=None):
    """Plot generalized additive model component, which includes first-degree rules
    and linear functions of unbinarized ordinal features but excludes higher-degree rules.

    Args:
        Xorig (DataFrame): Original unbinarized features
        fb: FeatureBinarizer object used to binarize features
        features (list, optional): Subset of features to be plotted
    """
    # Number of ordinal features used
    if model.useOrd:
        nnzOrd = len(model.idxNonzeroOrd)
    else:
        nnzOrd = 0

    # Initialize terms Series and x values for plots
    terms = pd.Series(
        index=pd.MultiIndex.from_arrays([[], [], []], names=model.z.index.names)
    )
    xPlot = {}

    # Iterate over ordinal features
    for i in range(nnzOrd):
        # Restrict to specified features
        f = model.namesOrd[model.idxNonzeroOrd[i]]
        if (features is not None) and (f not in features):
            continue
        # Append term
        terms = terms.append(pd.Series(model.lr.coef_[0, i], index=[(f, "", "")]))
        # Initialize x values with min and max
        xPlot[f] = [Xorig[f].min(), Xorig[f].max()]

    # Iterate over first-degree rules
    for i in range(model.z.shape[1]):
        if model.z.iloc[:, i].sum() > 1:
            continue
        # MultiIndex of rule
        idxTerm = model.z.index[model.z.iloc[:, i] > 0]
        (f, o, v) = idxTerm[0]
        # Restrict to specified features
        if (features is not None) and (f not in features):
            continue
        # Append new term
        terms = terms.append(pd.Series(model.lr.coef_[0, i + nnzOrd], index=idxTerm))
        # Update x values
        if f not in xPlot:
            if o in ["<=", ">"]:
                # Ordinal feature, initialize with min and max
                xPlot[f] = [Xorig[f].min(), Xorig[f].max()]
            else:
                # Categorical feature, use all values
                xPlot[f] = np.sort(Xorig[f].unique())
        if o in ["<=", ">"]:
            # Append values around threshold
            xPlot[f].extend([v - model.eps, v + model.eps])

    # Initialize y values for plots and variance calculation
    yPlot = {}
    yVar = pd.DataFrame(0.0, index=Xorig.index, columns=xPlot.keys())
    plotLine = pd.Series(False, index=xPlot.keys())
    # Iterate over GAM features
    for f in xPlot.keys():
        # Sort x values
        xPlot[f] = np.sort(np.array(xPlot[f]))
        yPlot[f] = np.zeros_like(xPlot[f], dtype=float)
        # Iterate over terms involving feature
        for ((o, v), c) in terms[f].iteritems():
            if o == "":
                if model.useOrd and (f in fb.ordinal):
                    # Add linear function of standardized feature with same factor of 0.4
                    idxf = fb.ordinal.index(f)
                    yPlot[f] += (
                        0.4
                        * c
                        * (xPlot[f] - fb.scaler.mean_[idxf])
                        / fb.scaler.scale_[idxf]
                    )
                    yVar[f] += (
                        0.4
                        * c
                        * (Xorig[f] - fb.scaler.mean_[idxf])
                        / fb.scaler.scale_[idxf]
                    )
                    plotLine[f] = True
                else:
                    # Binary feature, add indicator function
                    yPlot[f] += c * (xPlot[f] == fb.maps[f].index[1])
                    yVar[f] += c * (Xorig[f] == fb.maps[f].index[1])
            elif o == "<=":
                # Add step function
                yPlot[f] += c * (xPlot[f] <= v)
                yVar[f] += c * (Xorig[f] <= v)
                plotLine[f] = True
            elif o == ">":
                # Add step function
                yPlot[f] += c * (xPlot[f] > v)
                yVar[f] += c * (Xorig[f] > v)
                plotLine[f] = True
            elif o == "==":
                # Add indicator function
                yPlot[f] += c * (xPlot[f].astype(str) == v)
                yVar[f] += c * (Xorig[f].astype(str) == v)
            elif o == "!=":
                # Add indicator function
                yPlot[f] += c * (xPlot[f].astype(str) != v)
                yVar[f] += c * (Xorig[f].astype(str) != v)
            elif o == "not":
                # Binary feature, add indicator function
                yPlot[f] += c * (xPlot[f] == fb.maps[f].index[0])
                yVar[f] += c * (Xorig[f] == fb.maps[f].index[0])

    return xPlot, yPlot, plotLine
