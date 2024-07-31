import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import pydeck as pdk
from PIL import Image

from constants import NAME2COLOR

def build_deck(mode, pc_df, polygon_data):
    if mode == "first_person":
        view = pdk.View(type="FirstPersonView", controller=True)
        view_state = pdk.ViewState(latitude=0, longitude=0, bearing=-90, pitch=15)
        point_size = 10
    elif mode == "orbit":
        view = pdk.View(type="OrbitView", controller=True)
        view_state = pdk.ViewState(
            target=[0, 0, 1e-5],
            controller=True,
            zoom=23,
            rotation_orbit=-90,
            rotation_x=15,
        )
        point_size = 3

    else:
        view_state = pdk.ViewState(
            latitude=0,
            longitude=0,
            bearing=45,
            pitch=50,
            zoom=20,
            max_zoom=30,
            position=[0, 0, 1e-5],
        )
        view = pdk.View(type="MapView", controller=True)
        point_size = 1

    pc_layer = pdk.Layer(
        "PointCloudLayer",
        data=pc_df,
        get_position=["x", "y", "z"],
        get_color=[255, 255, 255],
        auto_highlight=True,
        pickable=False,
        point_size=point_size,
        coordinate_system=2,
        coordinate_origin=[0, 0],
    )

    box_layer = pdk.Layer(
        "PolygonLayer",
        data=polygon_data,
        stroked=True,
        pickable=True,
        filled=True,
        extruded=True,
        opacity=0.2,
        wireframe=True,
        line_width_min_pixels=1,
        get_polygon="polygon",
        get_fill_color="color",
        get_line_color=[255, 255, 255],
        get_line_width=0,
        coordinate_system=2,
        get_elevation="elevation",
    )

    tooltip = {"html": "<b>Label:</b> {name}"}

    r = pdk.Deck(
        [pc_layer, box_layer],
        initial_view_state=view_state,
        views=[view],
        tooltip=tooltip,
        map_provider=None,
    )

    return r


def compute_pointcloud_for_image(
    lv5,
    sample_token: str,
    dot_size: int = 2,
    pointsensor_channel: str = "LIDAR_TOP",
    camera_channel: str = "CAM_FRONT",
    out_path: str = None,
):
    """Scatter-plots a point-cloud on top of image.
    Args:
        sample_token: Sample token.
        dot_size: Scatter plot dot size.
        pointsensor_channel: RADAR or LIDAR channel name, e.g. 'LIDAR_TOP'.
        camera_channel: Camera channel name, e.g. 'CAM_FRONT'.
        out_path: Optional path to save the rendered figure to disk.
    Returns:
        tuple containing the points, array of colors and a pillow image
    """
    sample_record = lv5.get("sample", sample_token)

    # Here we just grab the front camera and the point sensor.
    pointsensor_token = sample_record["data"][pointsensor_channel]
    camera_token = sample_record["data"][camera_channel]

    points, coloring, im = lv5.explorer.map_pointcloud_to_image(
        pointsensor_token, camera_token
    )

    return points, coloring, im


def render_box_in_image(lv5, im, sample: str, camera_channel: str):
    camera_token = sample["data"][camera_channel]
    data_path, boxes, camera_intrinsic = lv5.get_sample_data(
        camera_token, flat_vehicle_coordinates=False
    )

    arr = np.array(im)

    for box in boxes:
        c = NAME2COLOR[box.name]
        box.render_cv2(arr, normalize=True, view=camera_intrinsic, colors=(c, c, c))

    new = Image.fromarray(arr)
    return new


def build_figure(lv5, sample, lidar, camera, overlay):
    points, coloring, im = compute_pointcloud_for_image(
        lv5, sample["token"], pointsensor_channel=lidar, camera_channel=camera
    )

    if "boxes" in overlay:
        im = render_box_in_image(lv5, im, sample, camera_channel=camera)

    fig = px.imshow(im, binary_format="jpeg", binary_compression_level=2)

    if "pointcloud" in overlay:
        fig.add_trace(
            go.Scattergl(
                x=points[0,],
                y=points[1,],
                mode="markers",
                opacity=0.4,
                marker_color=coloring,
                marker_size=3,
            )
        )

    fig.update_layout(
        margin=dict(l=10, r=10, t=0, b=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        hovermode=False,
    )
    fig.update_xaxes(showticklabels=False, showgrid=False, range=(0, im.size[0]))
    fig.update_yaxes(showticklabels=False, showgrid=False, range=(im.size[1], 0))

    return fig
