import vtk
import os
import dash
import dash_vtk
from dash_vtk.utils import to_mesh_state
import dash_html_components as html

head_vti = os.path.join("demos", "data", "head.vti")

# Load dataset from dist
reader = vtk.vtkXMLImageDataReader()
reader.SetFileName(head_vti)
reader.Update()

# Extract iso-mesh from image
contour = vtk.vtkContourFilter()
contour.SetInputConnection(reader.GetOutputPort())
contour.SetNumberOfContours(1)
contour.SetValue(0, 1500)
contour.Update()

# Get mesh to dash_vtk
mesh_state = to_mesh_state(contour.GetOutput())

content = dash_vtk.View(
    [dash_vtk.GeometryRepresentation([dash_vtk.Mesh(state=mesh_state)]),]
)

# Dash setup
app = dash.Dash(__name__)
server = app.server
app.layout = html.Div(
    style={"width": "100%", "height": "calc(100vh - 15px)"}, children=[content],
)

if __name__ == "__main__":
    app.run_server(debug=True)
