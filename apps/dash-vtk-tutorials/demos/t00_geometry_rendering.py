import dash
import dash_html_components as html

import dash_vtk
from dash_vtk.utils import to_mesh_state

from vtk.vtkImagingCore import vtkRTAnalyticSource

# Use VTK to get some data
data_source = vtkRTAnalyticSource()
data_source.Update()  # <= Execute source to produce an output
dataset = data_source.GetOutput()

# Use helper to get a mesh structure that can be passed as-is to a Mesh
# RTData is the name of the field
mesh_state = to_mesh_state(dataset)

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
