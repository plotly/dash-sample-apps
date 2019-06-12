import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go
# import pandas as pd

from utils import Header, make_dash_table

# df_fund_facts = pd.read_csv('../data/df_fund_facts.csv')
# df_price_perf = pd.read_csv('../data/df_price_perf.csv')
# df_current_prices = pd.read_csv('../data/df_current_prices.csv')
# df_hist_prices = pd.read_csv('../data/df_hist_prices.csv')
# df_avg_returns = pd.read_csv('../data/df_avg_returns.csv')
# df_after_tax = pd.read_csv('../data/df_after_tax.csv')
# df_recent_returns = pd.read_csv('../data/df_recent_returns.csv')
# df_equity_char = pd.read_csv('../data/df_equity_char.csv')
# df_equity_diver = pd.read_csv('../data/df_equity_diver.csv')
# df_expenses = pd.read_csv('../data/df_expenses.csv')
# df_minimums = pd.read_csv('../data/df_minimums.csv')
# df_dividend = pd.read_csv('../data/df_dividend.csv')
# df_realized = pd.read_csv('../data/df_realized.csv')
# df_unrealized = pd.read_csv('../data/df_unrealized.csv')
# df_graph = pd.read_csv("../data/df_graph.csv")

layout = html.Div(
    [  # page 6
        html.Div(
            [
                Header(),
                # Row 1
                html.Div(
                    [
                        html.Div(
                            [
                                html.H6(
                                    "News", className="gs-header gs-text-header padded"
                                ),
                                html.Br([]),
                                html.P(
                                    "10/25/16    The rise of indexing and the fall of costs"
                                ),
                                html.Br([]),
                                html.P(
                                    "08/31/16    It's the index mutual fund's 40th anniversary: Let the low-cost, passive party begin"
                                ),
                            ],
                            className="six columns",
                        ),
                        html.Div(
                            [
                                html.H6(
                                    "Reviews",
                                    className="gs-header gs-table-header padded",
                                ),
                                html.Br([]),
                                html.Li("Launched in 1976."),
                                html.Li(
                                    "On average, has historically produced returns that have far outpaced the rate of inflation.*"
                                ),
                                html.Li(
                                    "Quantitative Equity Group, the fund's advisor, is among the world's largest equity index managers."
                                ),
                                html.Br([]),
                                html.P(
                                    "Did you know? The fund launched in 1976 as First Index Investment Trust—the nation's first index fund available to individual investors."
                                ),
                                html.Br([]),
                                html.P(
                                    "* The performance of an index is not an exact representation of any particular investment, as you cannot invest directly in an index."
                                ),
                                html.Br([]),
                                html.P(
                                    "Past performance is no guarantee of future returns. See performance data current to the most recent month-end."
                                ),
                            ],
                            className="six columns",
                        ),
                    ],
                    className="row ",
                ),
            ],
            className="sub_page",
        )
    ],
    className="page",
)
