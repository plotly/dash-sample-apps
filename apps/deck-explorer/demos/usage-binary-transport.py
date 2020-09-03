"""
Adapted from: https://pydeck.gl/gallery/binary_transport.html

Example of binary transport in pydeck. This notebook renders 10k points
via the web sockets within a Jupyter notebook if you run with
``generate_vis(notebook_display=True)``

Since binary transfer relies on Jupyter's kernel communication,
note that the .html in the pydeck documentation does not use binary transfer
and is just for illustration.
"""
import os

import dash
import dash_deck
import dash_html_components as html
import pydeck
import pandas as pd

mapbox_api_token = os.getenv("MAPBOX_ACCESS_TOKEN")


NODES_URL = "https://raw.githubusercontent.com/ajduberstein/geo_datasets/master/social_nodes.csv"


def generate_graph_data(num_nodes, random_seed):
    """Generates a graph of 10k nodes with a 3D force layout

    This function is unused but serves as an example of how the data in
    this visualization was generated
    """
    import networkx as nx  # noqa

    g = nx.random_internet_as_graph(num_nodes, random_seed)
    node_positions = nx.fruchterman_reingold_layout(g, dim=3)

    force_layout_df = pd.DataFrame.from_records(node_positions).transpose()
    force_layout_df["group"] = [d[1]["type"] for d in g.nodes.data()]
    force_layout_df.columns = ["x", "y", "z", "group"]
    return force_layout_df


def make_renderer(nodes, use_binary_transport=False):
    """Creates the pydeck visualization for rendering"""
    view_state = pydeck.ViewState(
        offset=[0, 0], latitude=None, longitude=None, bearing=None, pitch=None, zoom=10,
    )

    views = [pydeck.View(type="OrbitView", controller=True)]

    nodes_layer = pydeck.Layer(
        "PointCloudLayer",
        nodes,
        get_position="position",
        get_normal=[10, 100, 10],
        get_color="color",
        pickable=True,
        # Set use_binary_transport to `True`
        use_binary_transport=use_binary_transport,
        auto_highlight=True,
        highlight_color=[255, 255, 0],
        radius=50,
    )

    return pydeck.Deck(
        layers=[nodes_layer],
        initial_view_state=view_state,
        views=views,
        map_provider=None,
    )


nodes = pd.read_csv(NODES_URL)

colors = pydeck.data_utils.assign_random_colors(nodes["group"])
# Divide by 255 to normalize the colors
# Specify positions and colors as columns of lists
nodes["color"] = nodes.apply(
    lambda row: [c / 255 if False else c for c in colors.get(row["group"])], axis=1
)
nodes["position"] = nodes.apply(lambda row: [row["x"], row["y"], row["z"]], axis=1)

# Remove all unused columns
del nodes["x"]
del nodes["y"]
del nodes["z"]
del nodes["group"]

r = make_renderer(nodes, use_binary_transport=False)

app = dash.Dash(__name__)

app.layout = html.Div(
    dash_deck.DeckGL(r.to_json(), id="deck-gl", style={"background-color": "charcoal"})
)


if __name__ == "__main__":
    app.run_server(debug=True)
