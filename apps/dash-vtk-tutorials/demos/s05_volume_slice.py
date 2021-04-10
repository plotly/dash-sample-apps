import dash
import random
import dash_html_components as html
import dash_vtk

field = [random.random() * i for i in range(10 * 10 * 10)]

content = dash_vtk.View(
    [
        dash_vtk.VolumeRepresentation(
            [
                dash_vtk.VolumeController(),
                dash_vtk.ShareDataSet(
                    [
                        dash_vtk.ImageData(
                            dimensions=[10, 10, 10],
                            spacing=[1, 1, 1],
                            origin=[-4, -4, -4],
                            children=[
                                dash_vtk.PointData(
                                    [
                                        dash_vtk.DataArray(
                                            registration="setScalars", values=field,
                                        )
                                    ]
                                )
                            ],
                        ),
                    ]
                ),
            ]
        ),
        dash_vtk.SliceRepresentation(
            property={"colorWindow": 500, "colorLevel": 200},
            iSlice=5,
            children=[dash_vtk.ShareDataSet()],
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
