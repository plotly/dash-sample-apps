import dash
import dash_html_components as html
import dash_vtk
import random

xyz = []
for i in range(10000):
    xyz.append(random.random())  # x
    xyz.append(random.random())  # y
    xyz.append(random.random() * 0.01)  # z

content = dash_vtk.View(
    children=[
        dash_vtk.GeometryRepresentation(
            property={"pointSize": 3},
            children=[dash_vtk.PolyData(points=xyz, connectivity="points")],
        ),
    ],
)

# Dash setup
app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(
    style={"width": "100%", "height": "calc(100vh - 15px)"}, children=[content],
)

if __name__ == "__main__":
    app.run_server(debug=True)
