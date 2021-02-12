import os

import dash
import dash_html_components as html
import itk

import dash_vtk
from dash_vtk.utils import to_volume_state

# Place a DICOM series (a set of per-file slices) in a directory. ITK sorts, sets spatial metadata, etc.
demo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
itk_image = itk.imread(os.path.join(demo_dir, "data", "ct_lung"))

# Convert itk.Image to vtkImageData
vtk_image = itk.vtk_image_from_image(itk_image)
volume_state = to_volume_state(vtk_image)


vtk_view = dash_vtk.View(
    dash_vtk.VolumeRepresentation(
        children=[dash_vtk.VolumeController(), dash_vtk.Volume(state=volume_state),]
    )
)

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(
    style={"height": "calc(100vh - 50px)", "width": "100%"},
    children=[html.Div(vtk_view, style={"height": "100%", "width": "100%"})],
)


if __name__ == "__main__":
    app.run_server(debug=True)
