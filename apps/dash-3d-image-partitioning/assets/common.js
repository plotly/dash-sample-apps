function json_copy(obj) {
    return JSON.parse(JSON.stringify(obj));
}

function graph_relayout_get_shapes(graph_relayout_data) {
    if (graph_relayout_data && ("shapes" in graph_relayout_data)) {
        return graph_relayout_data.shapes;
    }
}

function zip(a,b) {
    return a.map((x,i) => [x,b[i]]);
}
