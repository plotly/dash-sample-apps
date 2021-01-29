import os
import dash
import dash_html_components as html

import dash_vtk
from dash_vtk.utils import to_volume_state

import vtk

# Data file path
demo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
head_vti = os.path.join(demo_dir, "data", "head.vti")

# Load dataset from dist
reader = vtk.vtkXMLImageDataReader()
reader.SetFileName(head_vti)
reader.Update()

volume_state = to_volume_state(reader.GetOutput())

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(
    style={"height": "calc(100vh - 30px)", "width": "100%",},
    children=[
        dash_vtk.View(
            dash_vtk.VolumeRepresentation(
                children=[
                    dash_vtk.VolumeController(),
                    dash_vtk.Volume(state=volume_state),
                ]
            )
        )
    ],
)

if __name__ == "__main__":
    app.run_server(debug=True)
