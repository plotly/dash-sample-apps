import os
import json


def load_elements(path_name):
    dirname = os.path.dirname(__file__)
    path = os.path.join(dirname, path_name)

    with open(path, "r") as f:
        elements = json.loads(f.read())

    return elements


gene_elements = load_elements("../data/gene.json")
compound_elements = load_elements("../data/compound.json")
wine_and_cheese_elements = load_elements("../data/wine-and-cheese.json")
basic_elements = [
    {"data": {"id": "one", "label": "Node 1"}, "position": {"x": 50, "y": 50}},
    {"data": {"id": "two", "label": "Node 2"}, "position": {"x": 200, "y": 200}},
    {"data": {"id": "three", "label": "Node 3"}, "position": {"x": 100, "y": 150}},
    {"data": {"id": "four", "label": "Node 4"}, "position": {"x": 400, "y": 50}},
    {"data": {"id": "five", "label": "Node 5"}, "position": {"x": 250, "y": 100}},
    {
        "data": {"id": "six", "label": "Node 6", "parent": "three"},
        "position": {"x": 150, "y": 150},
    },
    {"data": {"source": "one", "target": "two", "label": "Edge from Node1 to Node2"}},
    {
        "data": {
            "source": "one",
            "target": "five",
            "label": "Edge from Node 1 to Node 5",
        }
    },
    {
        "data": {
            "source": "two",
            "target": "four",
            "label": "Edge from Node 2 to Node 4",
        }
    },
    {
        "data": {
            "source": "three",
            "target": "five",
            "label": "Edge from Node 3 to Node 5",
        }
    },
    {
        "data": {
            "source": "three",
            "target": "two",
            "label": "Edge from Node 3 to Node 2",
        }
    },
    {
        "data": {
            "source": "four",
            "target": "four",
            "label": "Edge from Node 4 to Node 4",
        }
    },
    {
        "data": {
            "source": "four",
            "target": "six",
            "label": "Edge from Node 4 to Node 6",
        }
    },
    {
        "data": {
            "source": "five",
            "target": "one",
            "label": "Edge from Node 5 to Node 1",
        }
    },
]

ELEMENTS = {
    "Basic": basic_elements,
    "Compound": compound_elements,
    "Gene": gene_elements,
    "WineandCheese": wine_and_cheese_elements,
}
ARROW_POSITIONS = ("source", "mid-source", "target", "mid-target")
LABEL_ELEMENT_TYPES = ("node", "edge")
LABEL_ELEMENT_TYPES_ALL = ("node", "edge", "source", "target")
