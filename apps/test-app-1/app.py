import dash
import dash_html_components as html

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div([
    html.H1('test. this is a test.')
])
