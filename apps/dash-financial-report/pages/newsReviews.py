import dash_html_components as html
from utils import Header


def create_layout(app):
    return html.Div(
        [
            Header(app),
            # page 6
            html.Div(
                [
                    # Row 1
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H6("News", className="subtitle padded"),
                                    html.Br([]),
                                    html.Div(
                                        [
                                            html.P(
                                                "10/25/16    The rise of indexing and the fall of costs"
                                            ),
                                            html.P(
                                                "08/31/16    It's the index mutual fund's 40th anniversary: Let the low-cost, passive party begin"
                                            ),
                                        ],
                                        style={"color": "#7a7a7a"},
                                    ),
                                ],
                                className="row",
                            ),
                            html.Div(
                                [
                                    html.H6("Reviews", className="subtitle padded"),
                                    html.Br([]),
                                    html.Div(
                                        [
                                            html.Li("Launched in 1976."),
                                            html.Li(
                                                "On average, has historically produced returns that have far outpaced the rate of inflation.*"
                                            ),
                                            html.Li(
                                                "Quantitative Equity Group, the fund's advisor, is among the world's largest equity index managers."
                                            ),
                                        ],
                                        id="reviews-bullet-pts",
                                    ),
                                    html.Div(
                                        [
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
                                        style={"color": "#7a7a7a"},
                                    ),
                                ],
                                className="row",
                            ),
                        ],
                        className="row ",
                    )
                ],
                className="sub_page",
            ),
        ],
        className="page",
    )
