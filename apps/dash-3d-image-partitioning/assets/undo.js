function* _range_gen(N) {
    let n=N; while(n>0){yield --n;}
}

// We cannot use this syntax because we want UndoState to be JSON serializable.
// All the methods below will work on objects with the following fields.
//function UndoState(n_shape_lists) {
//    self.undo_n_clicks=0;
//    self.redo_n_clicks=0;
//    self.undo_shapes=[];
//    self.redo_shapes=[];
//    // This has to be set depending on the shape. Here the example is for two
//    // figures, the shapes of the left figure in the first item of the list, the
//    // shapes of the right in the second item.
//    self.empty_shapes=[[],[]];
//}

// Check if clicked and update the number of clicks
function UndoState_undo_clicked(self, undo_n_clicks) {
    ret = (undo_n_clicks > self.undo_n_clicks);
    self.undo_n_clicks = undo_n_clicks;
    return ret;
}

// Check if clicked and update the number of clicks
function UndoState_redo_clicked(self, redo_n_clicks) {
    ret = (redo_n_clicks > self.redo_n_clicks);
    self.redo_n_clicks = redo_n_clicks;
    return ret;
}

// Return true if shapes are different from the last shapes stored in
// undo_shapes
function UndoState_shapes_changed (self,shapes) {
    return !(JSON.stringify(shapes)
             == JSON.stringify(self.undo_shapes[self.undo_shapes.length-1]));
}

// Add new shapes and empty the redo list (because if redos exist,
// they are now irrelevant after adding a new shape)
function UndoState_add_new_shapes (self,shapes) {
    self.redo_shapes = [];
    self.undo_shapes.push(shapes);
}

// Add new shapes to the undo list if shapes are new
function UndoState_track_changes (self,shapes) {
    if (UndoState_shapes_changed(self,shapes)) {
        UndoState_add_new_shapes(self,shapes);
    }
}

function UndoState__return_last_undo (self) {
    if (!self.undo_shapes.length) {
        return json_copy(self.empty_shapes);
    }
    let ret = self.undo_shapes[self.undo_shapes.length-1];
    return ret;
}    

// Undo returns previous shapes and stores the current shapes in the redo list
function UndoState_apply_undo (self) {
    if (self.undo_shapes.length) {
        self.redo_shapes.push(self.undo_shapes.pop());
    }
    return UndoState__return_last_undo(self);
}

// Redo puts last redo on top of the undo list and returns the set of shapes at
// the end of the undo list
function UndoState_apply_redo (self) {
    if (self.redo_shapes.length) {
        self.undo_shapes.push(self.redo_shapes.pop());
    }
    return UndoState__return_last_undo(self);
}
