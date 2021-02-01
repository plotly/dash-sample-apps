import random

import dash
import dash_vtk
import dash_html_components as html

random.seed(42)

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(
    style={"height": "calc(100vh - 30px)", "width": "100%",},
    children=[
        dash_vtk.View(
            children=dash_vtk.VolumeDataRepresentation(
                spacing=[1, 1, 1],
                dimensions=[10, 10, 10],
                origin=[0, 0, 0],
                scalars=[random.random() for _ in range(1000)],
                rescaleColorMap=False,
            )
        )
    ],
)

if __name__ == "__main__":
    app.run_server(debug=True)
