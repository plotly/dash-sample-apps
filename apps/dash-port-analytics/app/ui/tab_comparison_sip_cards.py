import dash_html_components as html
from config import strings


def make_port_comparison_sip_cards(
    port1: str,
    port1_inc: int,
    port1_dec: int,
    port1_trend: str,
    port2: str,
    port2_inc: int,
    port2_dec: int,
    port2_trend: str,
):
    """
    Returns HTML div of cards found on the bottom of the Compare tab.

    :param port1: str, current first port
    :param port1_inc: int, number of increases in the first port
    :param port1_dec: int, number of decreases in the first port
    :param port1_trend: str, trend of the first port
    :param port2: str, current second port
    :param port2_inc: int, number of increases in the second port
    :param port2_dec: int, number of decreases in the second port
    :param port2_trend: str, trend of the second port
    :return: HTML div
    """

    return html.Div(
        className="port-comparison-sip-cards-container",
        children=[
            html.Div(
                className="sip-card-stats-single",
                children=[
                    html.Div(
                        className="sip-card-title",
                        children=[
                            html.Img(src="assets/img/insights.svg"),
                            html.H3(
                                children=[f"{port1} {strings.COMPARISON_TITLE_END}"]
                            ),
                        ],
                    ),
                    html.Div(
                        className="sip-card-body",
                        children=[
                            html.Ul(
                                children=[
                                    html.Li(
                                        children=[
                                            f"{strings.COMPARISON_NUM_INCREASE} {port1_inc}"
                                        ]
                                    ),
                                    html.Li(
                                        children=[
                                            f"{strings.COMPARISON_NUM_DECREASE} {port1_dec}"
                                        ]
                                    ),
                                    html.Li(
                                        children=[
                                            f"{strings.COMPARISON_TREND} {port1_trend}"
                                        ]
                                    ),
                                ]
                            )
                        ],
                    ),
                ],
            ),
            html.Div(
                className="sip-card-stats-single",
                children=[
                    html.Div(
                        className="sip-card-title",
                        children=[
                            html.Img(src="assets/img/insights.svg"),
                            html.H3(
                                children=[f"{port2} {strings.COMPARISON_TITLE_END}"]
                            ),
                        ],
                    ),
                    html.Div(
                        className="sip-card-body",
                        children=[
                            html.Ul(
                                children=[
                                    html.Li(
                                        children=[
                                            f"{strings.COMPARISON_NUM_INCREASE} {port2_inc}"
                                        ]
                                    ),
                                    html.Li(
                                        children=[
                                            f"{strings.COMPARISON_NUM_DECREASE} {port2_dec}"
                                        ]
                                    ),
                                    html.Li(
                                        children=[
                                            f"{strings.COMPARISON_TREND} {port2_trend}"
                                        ]
                                    ),
                                ]
                            )
                        ],
                    ),
                ],
            ),
        ],
    )
