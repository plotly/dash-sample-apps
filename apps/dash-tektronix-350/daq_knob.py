import dash
import dash_daq as daq
import dash_html_components as html

app = dash.Dash(__name__)
app.config['suppress_callback_exceptions'] = True
server = app.server

app.layout = html.Div(
    children=[
        daq.Knob(
            value=1E5,
            id="frequency-input",
            label="Frequency (Hz)",
            labelPosition="bottom",
            scale={'interval': 1000000},
            max=2500000,
            min=1E5,
        )
    ]
)

if __name__ == '__main__':
    app.run_server(debug=True)
