// Track shapes in many "slices" using undo

// This takes the input from some graph.relayoutData events or undo/redo button
// press events and returns an array of arrays each containing an array of
// shapes. The intention is to be compatible with a system where there are
// multiple figures each showing a view of a slice of a 3D tensor. For example,
// if we have 2 figures and the first has 3 slices and the second 4, then the
// empty datastructure will look like [[[],[],[]],[[],[],[],[]]] where each
// entry is a variable-length array depending on the number of shapes in that
// slice+figure combination.
function undo_track_slice_figure_shapes (
// an array of relayoutData passed from graphs
graph_relayoutData,
// their corresponding DOM IDs
graph_relayoutData_ids,
// number of clicks of an undo button
undo_n_clicks,
// number of clicks of a redo button
redo_n_clicks,
// The undo data (ususally from a State variable wrapping a store)
undo_data,
// a store wrapped in a State variable storing the shapes for each figure and each slice
shapes,
// the active slice index for each figure
slice_indices,
// a function that takes a list of shapes and returns those that we want to
// track (for example if some shapes are to show some attribute but should not
// be tracked by undo/redo)
shape_filter,
) {
    let triggered = dash_clientside.callback_context.triggered.map(t => t['prop_id']);
    console.log("triggered"); console.log(triggered);
    console.log("graph_relayoutData"); console.log(graph_relayoutData);
    console.log("undo_data"); console.log(undo_data);
    console.log("shapes"); console.log(shapes);
    // Things that could happen:
    // Shape drawn / deleted in one of the graphs
    // undo or redo pressed
    // TODO we don't need the new figures, because we don't yet return the whole
    // figures, we just use the shapes already in the figure if we didn't get
    // any new figures passed in through the inputs to the callback
    // NOTE: the shapes in the figures (new or old) only correspond to a
    // single slice, but when undo or redo is pressed, we want to update the
    // shapes for all the slices.
    if (graph_relayoutData.concat([undo_n_clicks,redo_n_clicks]).map(
            x=>x!=undefined).reduce((a,v)=>a&&v,true)) {
        if (UndoState_undo_clicked(undo_data,undo_n_clicks)) {
            console.log("undo clicked");
            shapes = UndoState_apply_undo(undo_data);
        } else if (UndoState_redo_clicked(undo_data,redo_n_clicks)) {
            console.log("redo clicked");
            shapes = UndoState_apply_redo(undo_data);
        } else {
            console.log("shape drawn, no?");
            // When returning the new graph figure with the shapes, this doesn't
            // seem to change the contents of the relayoutData of the graph that
            // was not drawn on.
            // Also if the graph is zoomed, there won't be any shapes in graphA
            // or graphB, but we want the shapes to persist. So what we do is:
            // shapesA and shapesB are set to the current shapes in graphA_fig
            // and graphB_fig respectively. If graphA or graphB triggered the
            // callback and this data contains shapes, then we populate with the
            // new shapes. Otherwise the shapes from before are used.
            zip(graph_relayoutData_ids,graph_relayoutData).forEach(
                    function(t,i){
                    if ((triggered == t[0]) && ("shapes" in t[1])) {
                    shapes[i][slice_indices[i]] = shape_filter(t[1].shapes);
                    }
                    }
                    );
            UndoState_track_changes(undo_data,shapes);
        }
    }
    console.log("shapes",shapes);
    return [undo_data,shapes];
}
