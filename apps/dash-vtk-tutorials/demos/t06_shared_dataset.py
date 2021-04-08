import dash
import dash_html_components as html

import dash_vtk
from dash_vtk.utils import to_volume_state

from vtk.vtkImagingCore import vtkRTAnalyticSource

# Use VTK to get some data
data_source = vtkRTAnalyticSource()
data_source.Update()  # <= Execute source to produce an output
dataset = data_source.GetOutput()

# Use helper to get a volume structure that can be passed as-is to a Volume
volume_state = to_volume_state(dataset)  # No need to select field

content = dash_vtk.View(
    [
        dash_vtk.VolumeRepresentation(
            [
                # GUI to control Volume Rendering
                # + Setup good default at startup
                dash_vtk.VolumeController(),
                # Actual volume
                dash_vtk.ShareDataSet([dash_vtk.Volume(state=volume_state),]),
            ]
        ),
        dash_vtk.SliceRepresentation(iSlice=10, children=[dash_vtk.ShareDataSet(),]),
        dash_vtk.SliceRepresentation(jSlice=10, children=[dash_vtk.ShareDataSet(),]),
        dash_vtk.SliceRepresentation(kSlice=10, children=[dash_vtk.ShareDataSet(),]),
    ]
)

# Dash setup
app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(
    style={"width": "100%", "height": "calc(100vh - 15px)"}, children=[content],
)

if __name__ == "__main__":
    app.run_server(debug=True)
