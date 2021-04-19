import plotly.graph_objects as go

COLOR_WHITE = "rgba(0,0,0,0)"  # plot background
COLOR_GRID = "#35363E"
COLOR_APPSILON_1 = "#E87272"  # unspecyfied / compare left
COLOR_APPSILON_2 = "#7a0091"  # navigation
COLOR_APPSILON_3 = "#11498A"  # fishing
COLOR_APPSILON_4 = "#1A6D9B"  # tug
COLOR_APPSILON_5 = "#12A5B0"  # passenger
COLOR_APPSILON_6 = "#3A9971"  # Cargo
COLOR_APPSILON_7 = "#79BD00"  # tanker
COLOR_APPSILON_8 = "#DBB657"  # pleasure / compare right


def generate_plot_layout(x_title: str, y_title: str, bar_mode: str) -> go.Layout:
    """
    Creates a layout code for the charts

    :param x_title: str, X axis title
    :param y_title: str, Y axis title
    :param bar_mode: str, stacked or grouped bar chart
    :return: Plotly layout object
    """
    if bar_mode:
        return go.Layout(
            barmode=bar_mode,
            xaxis_title=x_title,
            yaxis_title=y_title,
            xaxis={"showgrid": False},
            yaxis={"gridcolor": COLOR_GRID},
            paper_bgcolor=COLOR_WHITE,
            plot_bgcolor=COLOR_WHITE,
            margin={"l": 0, "r": 0, "t": 0, "b": 0},
        )
    return go.Layout(
        xaxis_title=x_title,
        yaxis_title=y_title,
        xaxis={"showgrid": False},
        yaxis={"gridcolor": COLOR_GRID},
        paper_bgcolor=COLOR_WHITE,
        plot_bgcolor=COLOR_WHITE,
        margin={"l": 0, "r": 0, "t": 0, "b": 0},
    )
