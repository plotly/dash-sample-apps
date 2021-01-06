import numpy as np
import pandas as pd
from nilearn import connectome, plotting
from plotly import graph_objects as go


def compute_connectivity_matrix(masker, fc_data):
    time_series = []
    connectome_measure = connectome.ConnectivityMeasure(kind="correlation")

    for func, confounds in zip(fc_data.func, fc_data.confounds):
        time_series.append(masker.fit_transform(func, confounds=confounds))
    correlation_matrices = connectome_measure.fit_transform(time_series)
    mean_correlation_matrix = connectome_measure.mean_

    return mean_correlation_matrix


def make_circos(df, gap=0.1, focus_region=None):
    conn_mat = ~df.isnull().values
    node_names = df.columns
    n_nodes = len(node_names)
    if focus_region is not None:
        focus_idx = list(node_names).index(focus_region)
        conn_mask = np.ones(conn_mat.shape, dtype=bool)
        conn_mask[focus_idx, :] = False
        conn_mask[:, focus_idx] = False
        conn_mat[conn_mask] = False

    np.fill_diagonal(conn_mat, False)
    n_conn_nodes = np.sum(conn_mat, 0)
    node_layout = [
        {
            "id": region_id,
            "label": str(region_name),
            "color": "grey",
            "len": n_conn_nodes[region_id] + (n_conn_nodes[region_id] - 1) * gap
            if not n_conn_nodes[region_id] == 0
            else 1,
        }
        for region_id, region_name in enumerate(node_names)
    ]
    highlight_data = [
        {
            "name": str(region_id),
            "block_id": str(region_id),
            "start": 0,
            "end": n_conn_nodes[region_id] + (n_conn_nodes[region_id] - 1) * gap
            if not n_conn_nodes[region_id] == 0
            else 1,
            "color": "blue",
        }
        for region_id in range(n_nodes)
    ]

    # Compute the connections
    # TODO: explain why we go through all this trouble to get the connections to point to the right place
    conn_thr = np.argwhere(conn_mat)
    conn_order = np.zeros(conn_mat.shape)
    for row_idx in range(n_nodes):
        # Get the indices of above threshold connections for the current row
        conn_idx = conn_thr[conn_thr[:, 0] == row_idx, :][:, 1]
        # Replace these values with their order in the row
        conn_order[row_idx, conn_idx] = np.arange(len(conn_idx))

    # Make the connection data
    connection_data = []
    for (row_id, column_id) in zip(*np.tril_indices_from(df.values, -1)):
        conn_value = df.iloc[row_id, column_id]
        # Compute order of source and target
        source_pos = conn_order[row_id, column_id] + conn_order[row_id, column_id] * gap
        target_pos = (
            conn_order.T[row_id, column_id] + conn_order.T[row_id, column_id] * gap
        )
        if conn_mat[row_id, column_id]:
            conn_dict = {
                "color": "red",
                "source": {
                    "id": str(row_id),
                    "start": source_pos,
                    "end": source_pos + 1,
                },
                "target": {
                    "id": str(column_id),
                    "start": target_pos,
                    "end": target_pos + 1,
                },
                "value": df.iloc[row_id, column_id],
                "name": f"r = {conn_value:.3f} between {node_names[row_id]} and {node_names[column_id]}",
            }
            connection_data.append(conn_dict)
    return node_layout, highlight_data, connection_data


def make_heatmap(df, colormap=None):
    axis_ticks = df.columns
    n_elements = len(axis_ticks)
    gap_width = 20 / n_elements
    if n_elements > 100:
        gap_width = 0
    layout = go.Layout(
        yaxis=dict(
            autorange="reversed",
            scaleanchor="x",
            constrain="domain",
            scaleratio=1,
            showgrid=True,
            showticklabels=False,
            tickvals=np.arange(-0.5, df.shape[0], 1),
        ),
        xaxis=dict(
            constrain="domain",
            tickvals=np.arange(-0.5, df.shape[0], 1),
            showgrid=True,
            showticklabels=False,
        ),
        margin=dict(l=0, r=0, b=0, t=0, pad=0),
    )
    heatmap = go.Heatmap(
        z=df,
        zmin=0,
        zmax=1,
        x=axis_ticks,
        y=axis_ticks,
        xgap=gap_width,
        ygap=gap_width,
        hoverongaps=False,
        autocolorscale=True if colormap is None else False,
        colorscale=colormap,
    )
    fig = go.Figure(heatmap, layout)
    fig.update_traces(showscale=True)

    return fig


def add_region_shape(fig, region_id, line=None, orientation="h"):
    # We just use the position along the y axis here but because the matrix is symmetrical, there is no difference
    draw_pos = np.where(np.array(fig.data[0]["y"]) == region_id)[0][0]
    fig.add_shape(
        type="rect",
        xref="x",
        yref="y",
        x0=-0.5 if orientation == "h" else draw_pos - 0.5,
        y0=draw_pos - 0.5 if orientation == "h" else -0.5,
        x1=len(fig.data[0]["x"]) - 0.5 if orientation == "h" else draw_pos + 0.5,
        y1=draw_pos + 0.5 if orientation == "h" else len(fig.data[0]["y"]) - 0.5,
        line=line if line is not None else dict(color="LightSeaGreen", width=3),
    )
    return fig


def compute_region_table(parcellation_img, parcellation_labels):
    # This step alone takes about 1 second. Might need to move outside the function in the future.
    coordinates, region_ids = plotting.find_parcellation_cut_coords(
        parcellation_img, return_label_names=True
    )
    if len(parcellation_labels) > len(region_ids):
        parcellation_labels = parcellation_labels[1:]
    table = pd.DataFrame(
        {
            "Region ID": region_ids,
            "Region Name": parcellation_labels,
            "X": np.round(coordinates[:, 0], 2),
            "Y": np.round(coordinates[:, 1], 2),
            "Z": np.round(coordinates[:, 2], 2),
        }
    )
    return table
