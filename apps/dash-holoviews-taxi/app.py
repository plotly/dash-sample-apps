# Dash import
import dash
import dash_html_components as html
import dash_bootstrap_components as dbc

# HoloViews imports
import holoviews as hv
from holoviews.operation import histogram
from holoviews.operation.datashader import datashade
from holoviews.plotting.plotly.dash import to_dash
from holoviews.selection import link_selections

import pandas as pd
import plotly.io as pio
from plotly import colors

pio.templates.default = "plotly_white"

df = pd.read_parquet(
    "https://github.com/plotly/dash-holoviews-taxi/releases/download/v0.0.1a1/nyc_taxi_small.parq"
)
ds = hv.Dataset(df)

# Uncomment for CUDF support
# import cudf
# ds = hv.Dataset(cudf.from_pandas(df))

# Add more descriptive labels
ds = ds.redim.label(fare_amount="Fare Amount")

points = hv.Points(ds, ["dropoff_x", "dropoff_y"])
shaded = datashade(points, cmap=colors.sequential.Plasma)
tiles = hv.Tiles().opts(
    mapboxstyle="light",
    accesstoken="pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNrOWJqb2F4djBnMjEzbG50amg0dnJieG4ifQ.Zme1-Uzoi75IaFbieBDl3A",
    height=500,
    width=500,
    padding=0,
)

hist = histogram(
    ds, dimension="fare_amount", normed=False, num_bins=20, bin_range=(0, 30.0)
).opts(color=colors.qualitative.Plotly[0], height=500)

lnk_sel = link_selections.instance()
linked_map = lnk_sel(tiles * shaded)

linked_hist = lnk_sel(hist)

# Use plot hook to set the default drag mode to box selection
def set_dragmode(plot, element):
    fig = plot.state
    fig["layout"]["dragmode"] = "select"
    fig["layout"]["selectdirection"] = "h"


linked_hist.opts(hv.opts.Histogram(hooks=[set_dragmode]))

# Set plot margins, ordered (left, bottom, right, top)
linked_hist.opts(margins=(60, 40, 30, 30))
linked_map.opts(margins=(30, 30, 30, 30))

# Create Dash App
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

components = to_dash(
    app, [linked_map, linked_hist], reset_button=True, button_class=dbc.Button,
)

app.layout = dbc.Container(
    [
        html.H1("NYC Taxi Demo", style={"padding-top": 40}),
        html.H3("Crossfiltering 10 million trips with Dash, Datashader, and HoloViews"),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    children=[
                        dbc.Card(
                            [
                                dbc.CardHeader("Drop off locations"),
                                dbc.CardBody(children=[components.graphs[0],]),
                            ]
                        )
                    ]
                ),
                dbc.Col(
                    children=[
                        dbc.Card(
                            [
                                dbc.CardHeader("Fair Amount"),
                                dbc.CardBody(children=[components.graphs[1]]),
                            ]
                        )
                    ]
                ),
            ]
        ),
        html.Div(style={"margin-top": 10}, children=components.resets[0]),
        components.store,
    ]
)

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    app.run_server(debug=True, port=8068)
