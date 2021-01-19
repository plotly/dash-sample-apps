import dash

from editor.callbacks import assign_callbacks
from editor.layout import layout as cytoscape_layout

app = dash.Dash(__name__)
server = app.server

app.layout = cytoscape_layout
assign_callbacks(app)


if __name__ == "__main__":
    app.run_server(debug=True)
