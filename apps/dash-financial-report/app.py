# coding: utf-8

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go

from pages import (
    overview,
    pricePerformance,
    portfolioManagement,
    feesMins,
    distributions,
    newsReviews,
)
from utils import Header, make_dash_table

import pandas as pd

# import pathlib

app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)
server = app.server

# # get relative data folder
# PATH = pathlib.Path(__file__).parent
# DATA_PATH = PATH.joinpath("data").resolve()


# # read data for tables (one df per table)
# df_fund_facts = pd.read_csv(DATA_PATH.joinpath("df_fund_facts.csv"))
# df_price_perf = pd.read_csv(DATA_PATH.joinpath("df_price_perf.csv"))
# df_current_prices = pd.read_csv(DATA_PATH.joinpath("df_current_prices.csv"))
# df_hist_prices = pd.read_csv(DATA_PATH.joinpath("df_hist_prices.csv"))
# df_avg_returns = pd.read_csv(DATA_PATH.joinpath("df_avg_returns.csv"))
# df_after_tax = pd.read_csv(DATA_PATH.joinpath("df_after_tax.csv"))
# df_recent_returns = pd.read_csv(DATA_PATH.joinpath("df_recent_returns.csv"))
# df_equity_char = pd.read_csv(DATA_PATH.joinpath("df_equity_char.csv"))
# df_equity_diver = pd.read_csv(DATA_PATH.joinpath("df_equity_diver.csv"))
# df_expenses = pd.read_csv(DATA_PATH.joinpath("df_expenses.csv"))
# df_minimums = pd.read_csv(DATA_PATH.joinpath("df_minimums.csv"))
# df_dividend = pd.read_csv(DATA_PATH.joinpath("df_dividend.csv"))
# df_realized = pd.read_csv(DATA_PATH.joinpath("df_realized.csv"))
# df_unrealized = pd.read_csv(DATA_PATH.joinpath("df_unrealized.csv"))
# df_graph = pd.read_csv(DATA_PATH.joinpath("df_graph.csv"))

# Describe the layout/ UI, of the app
app.layout = html.Div(
    [dcc.Location(id="url", refresh=False), html.Div(id="page-content")]
)

# Update page
@app.callback(
    dash.dependencies.Output("page-content", "children"),
    [dash.dependencies.Input("url", "pathname")],
)
def display_page(pathname):
    if (
        pathname == "/dash-financial-report"
        or pathname == "/"
        or pathname == "/dash-financial-report/overview"
    ):
        return overview.create_layout(app)
    elif pathname == "/dash-financial-report/price-performance":
        return pricePerformance.create_layout(app)
    elif pathname == "/dash-financial-report/portfolio-management":
        return portfolioManagement.create_layout(app)
    elif pathname == "/dash-financial-report/fees":
        return feesMins.create_layout(app)
    elif pathname == "/dash-financial-report/distributions":
        return distributions.create_layout(app)
    elif pathname == "/dash-financial-report/news-and-reviews":
        return newsReviews.create_layout(app)
    elif pathname == "/dash-financial-report/full-view":
        return (
            overview.create_layout(app),
            pricePerformance.create_layout(app),
            portfolioManagement.create_layout(app),
            feesMins.create_layout(app),
            distributions.create_layout(app),
            newsReviews.create_layout(app),
        )
    else:
        return overview.create_layout(app)


# detail the way that external_css and external_js work and link to alternative method locally hosted
# external_css = ["https://codepen.io/bcd/pen/KQrXdb.css"]

external_css = [
    "https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
    "https://codepen.io/bcd/pen/KQrXdb.css",
    "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
]

for css in external_css:
    app.css.append_css({"external_url": css})

if __name__ == "__main__":
    app.run_server(debug=True)
