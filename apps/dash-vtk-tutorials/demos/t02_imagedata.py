import dash
import dash_html_components as html

import dash_vtk

content = dash_vtk.View(
    [
        dash_vtk.VolumeRepresentation(
            [
                # GUI to control Volume Rendering
                # + Setup good default at startup
                dash_vtk.VolumeController(),
                # Actual Imagedata
                dash_vtk.ImageData(
                    dimensions=[5, 5, 5],
                    origin=[-2, -2, -2],
                    spacing=[1, 1, 1],
                    children=[
                        dash_vtk.PointData(
                            [
                                dash_vtk.DataArray(
                                    registration="setScalars",
                                    values=list(range(5 * 5 * 5)),
                                )
                            ]
                        )
                    ],
                ),
            ]
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
