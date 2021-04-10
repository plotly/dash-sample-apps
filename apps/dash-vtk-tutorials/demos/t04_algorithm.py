import os
import dash
import dash_html_components as html

import dash_vtk

content = dash_vtk.View(
    [
        dash_vtk.GeometryRepresentation(
            mapper={
                "colorByArrayName": "layer",
                "scalarMode": 4,
                "interpolateScalarsBeforeMapping": False,
            },
            colorMapPreset="jet",
            colorDataRange=[0.2, 0.9],
            children=[
                dash_vtk.Algorithm(
                    vtkClass="vtkConcentricCylinderSource",
                    state={
                        "height": 0.25,
                        "radius": [0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1],
                        "cellFields": [0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1],
                        "mask": [1, 0, 1, 0, 1, 0, 1, 1],
                        "resolution": 60,
                        "skipInnerFaces": True,
                        "startTheta": 45,
                        "endTheta": 315,
                        "center": [0, 0, 0.5],
                    },
                ),
            ],
        ),
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
