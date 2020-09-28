
// function that updates the figures as the slice idices, show checkbox or found
// segmentations changes
// returns a list of figures
function figure_display_update (
    // an array of integers indicating the slice index for each figure
    slice_indices,
    // a string that will show the shapes and predicted segmentation when set to
    // "show", and no shapes or segmentation otherwise
    show_seg_check,
    // the found (ones the user has selected) segments to show, same size as
    // image_slices_data
    found_segs_data,
    // an array of length equal to the number of figures, each containing an
    // array equal to the number of slices for that view. Each item of this
    // array is a uri encoded image
    image_slices_data,
    // an array containing the figures, in the same order as slice_indices
    image_display_figures,
    // the segments to show (as determined by a superpixel algorithm), same size
    // as image_slices_data
    seg_slices_data,
    // an array the same size as image_slices_data but each element contains an
    // array of shapes to draw on the figure (if show_seg_check=='show')
    drawn_shapes_data) {
    console.log("image_display_figures");
    console.log(image_display_figures);
    let image_display_figures_ = image_display_figures.map(function(f,i) {
        var f_ = json_copy(f),
            image_select_value = ((slice_indices[i]==undefined)?0:slice_indices[i]);
        f_.layout.images[0].source=image_slices_data[i][image_select_value];
        if (show_seg_check == 'show') {
            f_.layout.images[1].source=seg_slices_data[i][image_select_value];
            f_.layout.images[2].source=found_segs_data[i][image_select_value];
        } else {
            f_.layout.images[1].source="";
            f_.layout.images[2].source="";
        }
        f_.layout.shapes=drawn_shapes_data[i][image_select_value];
        return f_;
    });
    console.log("image_display_figures_");
    console.log(image_display_figures_);
    return image_display_figures_;
}

// Return a triangle shape
// x,y are coordinates of center of triangle
// w is width, h height,
// o is orientation, which is the direction the point is pointing
// (perpendicular to width dimension), possible values are 'up', 'down',
// 'left', 'right'
function tri_shape (x,y,w=2.5,h=2.5,o='down',path_args={
       'line':{'width':2, 'color':"DarkOrange"},
       'fillcolor':"DarkOrange"}) {
    // TODO: I believe the positive direction for the y-coordinates is downward
    // on the screen because we're working with images (the y-axis is reversed
    // otherwise the image would be upside-down)
    let m0,c1,c2;
    if (o=='up') {
        m0=[x-w*0.5,y+h*0.5];
        c1=[x,y-0.5*h];
        c2=[x+w*0.5,y+0.5*h];
    } else if (o=='down') {
        m0=[x-w*0.5,y-h*0.5];
        c1=[x,y+0.5*h];
        c2=[x+w*0.5,y-0.5*h];
    } else if (o=='right') {
        m0=[x-h*0.5,y-w*0.5];
        c1=[x+h*0.5,y];
        c2=[x-h*0.5,y+0.5*w];
    } else if (o=='left') {
        m0=[x+h*0.5,y-w*0.5];
        c1=[x-h*0.5,y];
        c2=[x+h*0.5,y+0.5*w];
    } else {
        throw 'bad orientation: ' + o;
    }
    path=`M ${m0[0]} ${m0[1]} `;
    [c1,c2].forEach(c=> path+=`L ${c[0]} ${c[1]} `);
    path+="Z";
    return {
       'type':"path",
       'path':path,
       ...path_args
    };
}
