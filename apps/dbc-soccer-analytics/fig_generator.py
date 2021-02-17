import plotly.io as pio

# This script takes the motion_graph output json file and converts it back into a figure which can be
# displayed as an animated graph using the main app.py Dash app

def fig_from_json(filename):
    with open(filename, 'r') as f:
        fig = pio.from_json(f.read())
    return fig
