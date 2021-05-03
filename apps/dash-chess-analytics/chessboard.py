import plotly.graph_objects as go
import pandas as pd
import numpy as np
import plotly.express as px

# Define function to output an 8*8 dataframe based on a vector of coordinates.
def board_output(vector):
    brd = np.zeros((8, 8))
    for tup in vector:
        brd[tup] += 1

    return pd.DataFrame(brd)


def getStackedBar(dictionary):
    print(dictionary)
    fig = px.bar(
        pd.DataFrame({"Games": dictionary}).T * 100,
        height=50,
        orientation="h",
        barmode="stack",
        color_discrete_map={"BLACK": "black", "WHITE": "white", "DRAW": "gray"},
    )

    margin = 0
    fig.update_layout(
        xaxis_title="",
        yaxis_title="",
        xaxis_visible=False,
        yaxis_visible=False,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        legend_orientation="h",
        legend_itemwidth=50,
        legend_title="",
        legend_xanchor="auto",
        legend_x=0.5,
        margin=dict(l=margin, r=margin, t=margin, b=margin, pad=0),
        legend_font=dict(family="Arial", size=12, color="black"),
    )

    fig.update_traces(marker_line_width=0, hovertemplate="%{x}%")
    return fig


def getChessboard(dimensions: int = 600, margin: int = 50):
    row = [0, 1] * 4
    chessboard = go.Figure(
        layout=dict(
            margin=dict(l=margin, r=margin, t=margin, b=margin, pad=10),
            width=dimensions,
            height=dimensions,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            font_color="#303030",
            coloraxis_showscale=False,
            showlegend=False,
            yaxis=dict(
                range=[-0.5, 7.5],
                tickfont_size=12,
                fixedrange=True,
                tickmode="array",
                tickvals=list(range(0, 8)),
                ticktext=list(range(1, 9)),
            ),
            xaxis=dict(
                range=[-0.5, 7.5],
                tickfont_size=12,
                fixedrange=True,
                tickmode="array",
                tickvals=list(range(0, 8)),
                ticktext=["A", "B", "C", "D", "E", "F", "G", "H"],
            ),
        )
    )
    getBoard(chessboard)
    return chessboard


def getBoard(fig):
    size = 87.5
    x_black = []
    for i in range(8):
        if (i % 2) == 0:
            x_black += list(range(0, 8, 2))
        else:
            x_black += list(range(1, 8, 2))

    x_white = []
    for i in range(8):
        if (i % 2) == 0:
            x_white += list(range(1, 8, 2))
        else:
            x_white += list(range(0, 8, 2))
    fig.add_trace(
        go.Scatter(
            x0=0,
            y0=0,
            dx=0,
            x=x_black,
            y=sorted([i for i in range(0, 8, 1)] * 4),
            name="black_squares",
            mode="markers",
            opacity=1,
            marker_symbol="square",
            marker_line_color="black",
            marker_size=size,
            marker_sizemode="diameter",
            marker_opacity=1,
            marker_color="#303030",
            hoverinfo="none",
        )
    )

    fig.add_trace(
        go.Scatter(
            x0=0,
            y0=0,
            dx=0,
            x=x_white,
            y=sorted([i for i in range(0, 8, 1)] * 4),
            name="white_squares",
            mode="markers",
            opacity=1,
            marker_symbol="square",
            marker_size=size,
            marker_sizemode="diameter",
            marker_opacity=1,
            marker_color="white",
            hoverinfo="none",
        )
    )


def getHeatmap(dataframe: pd.DataFrame):
    """DataFrame must have columns named:
    rows => 1 to 8
    letters => A to H
    freq => frequency"""
    if dataframe["freq"].sum() == 0:
        freq = dataframe["freq"]
    else:
        freq = np.round(dataframe["freq"] / dataframe["freq"].sum() * 100, 2)
    heatmap = go.Scatter(
        x0=0,
        y0=0,
        dx=0,
        x=dataframe["cols"],
        y=dataframe["rows"],
        name="",
        mode="markers",
        opacity=1,
        marker_symbol="square",
        marker_line_color="#c12917",
        marker_size=freq,
        marker_sizeref=freq.max() / 60,
        marker_sizemin=0,
        marker_sizemode="diameter",
        marker_opacity=1,
        marker_color="#c12917",
        hovertext=dataframe["freq"],
        hoverinfo="all",
        # TODO
        hovertemplate="<b># Games:</b> %{hovertext}<extra></extra>",
    )
    return heatmap
