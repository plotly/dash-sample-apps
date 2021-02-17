import os

import dash
import dash_html_components as html
import dash_bootstrap_components as dbc
import itk

import dash_vtk
from dash_vtk.utils import to_volume_state


def dcm_to_volume(dir_path):
    itk_image = itk.imread(dir_path)
    vtk_image = itk.vtk_image_from_image(itk_image)
    volume_state = to_volume_state(vtk_image)

    return volume_state


# Place a DICOM series (a set of per-file slices) in a directory. ITK sorts, sets spatial metadata, etc.
demo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
volume_state = dcm_to_volume(os.path.join(demo_dir, "data", "mri_pancreas"))

app = dash.Dash(__name__)
server = app.server

vtk_view = dash_vtk.View(
    dash_vtk.VolumeRepresentation(
        [dash_vtk.VolumeController(), dash_vtk.Volume(state=volume_state)]
    )
)

app.layout = html.Div(
    style={"height": "calc(100vh - 50px)", "width": "100%"},
    children=[html.Div(style={"height": "100%", "width": "100%"}, children=vtk_view)],
)


if __name__ == "__main__":
    app.run_server(debug=True)
