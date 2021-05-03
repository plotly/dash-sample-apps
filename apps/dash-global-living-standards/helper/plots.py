import plotly.graph_objects as go
from plotly.subplots import make_subplots

import re


# from app import city_info_num, city_info_num_agg

# city_info_num = city_info_num
# city_info_num_agg = city_info_num_agg


def custom_dims_plot_deprecated(
    location, dims_selected, city_info_num, city_info_num_agg
):
    dims_selected = [re.sub("high |low ", "", dim) for dim in dims_selected]

    vals_city = city_info_num.loc[
        city_info_num["City"] == location, dims_selected
    ].values.tolist()[0]
    vals_agg = city_info_num_agg[dims_selected].values.tolist()

    data = [
        go.Bar(
            y=dims_selected,
            x=vals_agg,
            name="Median",
            orientation="h",
            marker=dict(color="rgba(58, 71, 80, 0.6)"),
        ),
        go.Bar(
            y=dims_selected,
            x=vals_city,
            name="Selected",
            orientation="h",
            marker=dict(color="rgba(246, 78, 139, 0.6)"),
        ),
    ]

    fig = go.Figure(data)

    fig.update_layout(barmode="group")
    return fig


def custom_dims_plot(
    location,
    dims_selected,
    city_info_num,
    city_info_num_agg,
    dimension_mapper,
    width,
    height,
):
    dims_selected = [re.sub("high |low ", "", dim) for dim in dims_selected]

    dimnames = [dimension_mapper[dim] for dim in dims_selected]

    vals_city = city_info_num.loc[
        city_info_num["City"] == location, dims_selected
    ].values.tolist()[0]
    vals_agg = city_info_num_agg[dims_selected].values.tolist()

    fig = make_subplots(rows=len(dims_selected), cols=1, subplot_titles=dimnames)
    legend = [False for _ in range(len(dims_selected))]
    legend[0] = True

    median_name = "Average"
    # selected_name = 'Selected city'

    for idx, dim in enumerate(dimnames):
        # crate traces
        trace1 = go.Bar(
            y=[dim],
            x=[vals_agg[idx]],
            name=median_name,
            legendgroup=median_name,
            orientation="h",
            marker=dict(color="#986EA8"),  # dark grey #333333
            showlegend=legend[idx],
            opacity=0.7,
        )
        trace2 = go.Bar(
            y=[dim],
            x=[vals_city[idx]],
            name=location,
            legendgroup=location,
            orientation="h",
            marker=dict(color="#F3D576"),  # purple: #986EA8
            showlegend=legend[idx],
            opacity=0.7,
        )

        fig.add_trace(trace2, row=idx + 1, col=1)
        fig.add_trace(trace1, row=idx + 1, col=1)

    fig.update_layout(
        barmode="overlay",  # overlay, group
        legend=dict(
            orientation="h"
            # ,
            # yanchor="bottom",
            # y=1.02,
            # xanchor="right",
            # x=1
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="Open sans", size=10, color="White"),
        height=int(height * 0.35),
        width=int(width * 0.20),
        margin=dict(
            l=25,  # left margin
            r=25,  # right margin
            b=0,  # bottom margin
            t=25,  # top margin
        ),
    )
    fig.update_yaxes(visible=False, gridcolor="rgba(0,0,0,0)")
    fig.update_xaxes(gridcolor="rgba(0,0,0,0)")

    return fig
