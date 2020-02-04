import math
import dash
from dash.dependencies import Input, Output
import dash_cytoscape as cyto
import dash_html_components as html
import pathlib

from Bio import Phylo

app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)
server = app.server


def generate_elements(tree, xlen=30, ylen=30, grabbable=False):
    def get_col_positions(tree, column_width=80):
        """Create a mapping of each clade to its column position."""
        taxa = tree.get_terminals()

        # Some constants for the drawing calculations
        max_label_width = max(len(str(taxon)) for taxon in taxa)
        drawing_width = column_width - max_label_width - 1

        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
            # Potential drawing overflow due to rounding -- 1 char per tree layer
        fudge_margin = int(math.ceil(math.log(len(taxa), 2)))
        cols_per_branch_unit = (drawing_width - fudge_margin) / float(
            max(depths.values())
        )
        return dict(
            (clade, int(blen * cols_per_branch_unit + 1.0))
            for clade, blen in depths.items()
        )

    def get_row_positions(tree):
        taxa = tree.get_terminals()
        positions = dict((taxon, 2 * idx) for idx, taxon in enumerate(taxa))

        def calc_row(clade):
            for subclade in clade:
                if subclade not in positions:
                    calc_row(subclade)
            positions[clade] = (
                positions[clade.clades[0]] + positions[clade.clades[-1]]
            ) // 2

        calc_row(tree.root)
        return positions

    def add_to_elements(clade, clade_id):
        children = clade.clades

        pos_x = col_positions[clade] * xlen
        pos_y = row_positions[clade] * ylen

        cy_source = {
            "data": {"id": clade_id},
            "position": {"x": pos_x, "y": pos_y},
            "classes": "nonterminal",
            "grabbable": grabbable,
        }
        nodes.append(cy_source)

        if clade.is_terminal():
            cy_source["data"]["name"] = clade.name
            cy_source["classes"] = "terminal"

        for n, child in enumerate(children):
            # The "support" node is on the same column as the parent clade,
            # and on the same row as the child clade. It is used to create the
            # 90 degree angle between the parent and the children.
            # Edge config: parent -> support -> child

            support_id = clade_id + "s" + str(n)
            child_id = clade_id + "c" + str(n)
            pos_y_child = row_positions[child] * ylen

            cy_support_node = {
                "data": {"id": support_id},
                "position": {"x": pos_x, "y": pos_y_child},
                "grabbable": grabbable,
                "classes": "support",
            }

            cy_support_edge = {
                "data": {
                    "source": clade_id,
                    "target": support_id,
                    "sourceCladeId": clade_id,
                }
            }

            cy_edge = {
                "data": {
                    "source": support_id,
                    "target": child_id,
                    "length": clade.branch_length,
                    "sourceCladeId": clade_id,
                }
            }

            if clade.confidence and clade.confidence.value:
                cy_source["data"]["confidence"] = clade.confidence.value

            nodes.append(cy_support_node)
            edges.extend([cy_support_edge, cy_edge])

            add_to_elements(child, child_id)

    col_positions = get_col_positions(tree)
    row_positions = get_row_positions(tree)

    nodes = []
    edges = []

    add_to_elements(tree.clade, "r")

    return nodes, edges


PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()

# Define elements, stylesheet and layout
tree = Phylo.read(DATA_PATH.joinpath("apaf.xml"), "phyloxml")
nodes, edges = generate_elements(tree)
elements = nodes + edges


stylesheet = [
    {
        "selector": ".nonterminal",
        "style": {
            "label": "data(confidence)",
            "background-opacity": 0,
            "text-halign": "left",
            "text-valign": "top",
        },
    },
    {"selector": ".support", "style": {"background-opacity": 0}},
    {
        "selector": "edge",
        "style": {
            "source-endpoint": "inside-to-node",
            "target-endpoint": "inside-to-node",
        },
    },
    {
        "selector": ".terminal",
        "style": {
            "label": "data(name)",
            "width": 10,
            "height": 10,
            "text-valign": "center",
            "text-halign": "right",
            "background-color": "#222222",
        },
    },
]

app.layout = html.Div(
    [
        html.Img(className="logo", src=app.get_asset_url("dash-logo.png")),
        html.Div(
            className="header",
            children=[
                html.Div(
                    className="div-info",
                    children=[
                        html.H2(className="title", children="Cytoscape Phylogeny"),
                        html.P(
                            """
                            Dash Cytoscape is a graph visualization component for creating easily customizable,
                            high-performance interactive, and web-based networks.
                            """
                        ),
                        html.A(
                            children=html.Button("Learn More", className="button"),
                            href="https://dash.plot.ly/cytoscape",
                            target="_blank",
                        ),
                    ],
                ),
                html.H4("Phylogeny"),
                cyto.Cytoscape(
                    id="cytoscape",
                    className="cytoscape",
                    elements=elements,
                    stylesheet=stylesheet,
                    layout={"name": "preset", "fit": True, "animate": True},
                    minZoom=0.25,
                ),
            ],
        ),
    ]
)


@app.callback(
    Output("cytoscape", "stylesheet"), [Input("cytoscape", "mouseoverEdgeData")]
)
def color_children(edgeData):
    if edgeData is None:
        return stylesheet

    if "s" in edgeData["source"]:
        val = edgeData["source"].split("s")[0]
    else:
        val = edgeData["source"]

    children_style = [
        {"selector": f'edge[source *= "{val}"]', "style": {"line-color": "#3ed6d2"}}
    ]

    return stylesheet + children_style


if __name__ == "__main__":
    app.run_server(debug=True)
