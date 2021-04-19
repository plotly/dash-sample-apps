import dash
import dash_html_components as html

import dash_vtk

content = dash_vtk.View(
    [
        dash_vtk.GeometryRepresentation(
            children=[
                dash_vtk.PolyData(
                    points=[0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0,],
                    lines=[3, 1, 3, 2],
                    polys=[3, 0, 1, 2],
                    children=[
                        dash_vtk.PointData(
                            [
                                dash_vtk.DataArray(
                                    # registration='setScalars', # To activate field
                                    name="onPoints",
                                    values=[0, 0.33, 0.66, 1],
                                )
                            ]
                        ),
                        dash_vtk.CellData(
                            [
                                dash_vtk.DataArray(
                                    # registration='setScalars', # To activate field
                                    name="onCells",
                                    values=[0, 1],
                                )
                            ]
                        ),
                    ],
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
