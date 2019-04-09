import dash
import dash_html_components as html
from shared_code import hello

hello()

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div([
    html.H1('Simple')
])
