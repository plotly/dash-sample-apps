import plotly.graph_objs as go
import plotly.express as px
from sklearn.metrics import roc_curve, confusion_matrix, roc_auc_score, accuracy_score
import numpy as np
import pandas as pd

mesh_size = 0.02


def prediction_plot(**kwargs):

    _data = kwargs['data']
    split_data = _data[0]
    model = kwargs['model']

    y_pred_train = (model[0].decision_function(split_data[0]) >
                    kwargs['threshold']).astype(int)
    y_pred_test = (model[0].decision_function(split_data[1]) >
                   kwargs['threshold']).astype(int)
    train_score = accuracy_score(y_true=split_data[2].astype(int),
                                 y_pred=y_pred_train)
    test_score = accuracy_score(y_true=split_data[3].astype(int),
                                y_pred=y_pred_test)

    scaled_threshold = kwargs['threshold'] * (model[1].max() -
                                              model[1].min()) + model[1].min()

    range = max(abs(scaled_threshold - model[1].min()),
                abs(scaled_threshold - model[1].max()))

    trace0 = go.Contour(x=np.arange(model[2].min(), model[2].max(), mesh_size),
                        y=np.arange(model[3].min(), model[3].max(), mesh_size),
                        z=model[1].reshape(model[2].shape),
                        zmin=scaled_threshold - range,
                        zmax=scaled_threshold + range,
                        hoverinfo='none',
                        showscale=False,
                        contours=dict(showlines=False),
                        colorscale='rdgy',
                        opacity=0.6)

    trace1 = go.Contour(x=np.arange(model[2].min(), model[2].max(), mesh_size),
                        y=np.arange(model[3].min(), model[3].max(), mesh_size),
                        z=model[1].reshape(model[2].shape),
                        showscale=False,
                        hoverinfo='none',
                        contours=dict(
                            showlines=False,
                            type='constraint',
                            operation='=',
                            value=scaled_threshold,
                        ),
                        name=f'Threshold ({scaled_threshold:.3f})',
                        line=dict(color='#454545'))

    trace2 = go.Scatter(x=split_data[0][:, 0],
                        y=split_data[0][:, 1],
                        mode='markers',
                        name=f'Training Data (accuracy={train_score:.3f})',
                        marker=dict(size=10,
                                    color=split_data[2].astype(int),
                                    colorscale='tealrose',
                                    line=dict(width=1)))

    trace3 = go.Scatter(x=split_data[1][:, 0],
                        y=split_data[1][:, 1],
                        mode='markers',
                        name=f'Test Data (accuracy={test_score:.3f})',
                        marker=dict(
                            size=10,
                            symbol='triangle-up',
                            color=split_data[3].astype(int),
                            colorscale='tealrose',
                            line=dict(width=1),
                        ))

    layout = go.Layout(xaxis=dict(
        ticks='',
        showticklabels=False,
        showgrid=False,
        zeroline=False,
    ),
                       transition=dict(easing='exp-in-out',
                                       ordering="traces first",
                                       duration=500),
                       yaxis=dict(
                           ticks='',
                           showticklabels=False,
                           showgrid=False,
                           zeroline=False,
                       ),
                       plot_bgcolor="#fff",
                       paper_bgcolor="#fff",
                       hovermode='closest',
                       legend=dict(x=0, y=-0.01, orientation="h"),
                       margin=dict(l=0, r=0, t=0, b=0))

    fig = go.Figure(data=[trace0, trace1, trace2, trace3], layout=layout)

    return fig


def roc_curve_plot(**kwargs):

    _data = kwargs['data']
    split_data = _data[0]
    model = kwargs['model']

    y_score = model[0].decision_function(_data[1])
    fpr, tpr, thresholds = roc_curve(_data[2], y_score)

    auc_score = roc_auc_score(y_true=_data[2], y_score=y_score)

    fig = px.line(x=fpr, y=tpr)
    fig.update_traces(hovertemplate=None, line_color='rgb(49,130,189)')
    fig.update_layout(
        title={
            'text': f'ROC Curve (AUC = {auc_score:.3f})',
            'y': 0.5,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'bottom'
        },
        transition=dict(easing='cubic-in-out', duration=500),
        yaxis=dict(  #range=[0, 1],
            title='True Positive Rate',
            scaleanchor="x",
            scaleratio=1),
        xaxis=dict(  #range=[0, 1],
            title='False Positive Rate', constrain='domain'),
        hovermode='closest',
        height=400,
        showlegend=False,
        plot_bgcolor='#FAF9DE',
        margin=dict(l=10, r=0, t=50, b=20))

    return fig


def confusion_matrix_plot(**kwargs):

    _data = kwargs['data']
    split_data = _data[0]
    model = kwargs['model']

    scaled_threshold = kwargs['threshold'] * (model[1].max() -
                                              model[1].min()) + model[1].min()
    y_pred_test = (model[0].decision_function(split_data[1]) >
                   scaled_threshold).astype(int).astype(str)

    matrix = confusion_matrix(y_true=split_data[3], y_pred=y_pred_test)
    mtx = matrix / matrix.sum()

    label_text = [["True Negative", "False Positive"],
                  ["False Negative", "True Positive"]]

    fig = px.imshow(mtx,
                    x=['X', 'y'],
                    y=['X', 'y'],
                    color_continuous_scale='sunsetdark',
                    zmin=0,
                    zmax=1,
                    aspect="auto")

    fig.update_traces(text=label_text,
                      texttemplate="%{text}",
                      name='',
                      customdata=matrix,
                      hovertemplate='%{customdata:,}')

    fig.update_layout(xaxis_title="TRAIN",
                      yaxis_title="TEST",
                      transition=dict(easing='sin-in-out', duration=500),
                      height=400,
                      margin=dict(l=10, r=20, t=50, b=20))

    return fig
