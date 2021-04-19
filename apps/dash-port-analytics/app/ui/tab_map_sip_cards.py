import pandas as pd
import dash_html_components as html
from app import helpers
from config import strings


def make_tab_port_map_sip_card(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
) -> html.Div:
    """
    Returns a container filled with card for current filter condition. Shows the number of vessels in port,
    number of sailing, number moored, and when the data was last updated.

    :param df: pandas DataFrame, data source
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: HTML div, a container of 4 individual cards
    """

    def make_single_sip_card(icon_class: str, title: str, value: str) -> html.Div:
        """
        Makes a single card.

        :param icon_class: str, FontAwesome icon class
        :param title: str, name of the card
        :param value: str, value of the card
        :return: HTML div, a single card
        """
        return html.Div(
            className="sip-card-single",
            children=[
                html.Div(
                    className="sip-card-single-left",
                    children=[html.I(className=f"fas {icon_class}")],
                ),
                html.Div(
                    className="sip-card-single-right",
                    children=[
                        html.P(children=[f"{title}"]),
                        html.H3(children=[f"{value}"]),
                    ],
                ),
            ],
        )

    data = helpers.filter_by_port_vessel_and_time(
        df=df, port=port, vessel_type=vessel_type, year=year, month=month
    )
    n_ships = len(data)
    n_ships_stop = len(data[data["is_parked"] == 1])
    n_ships_sailing = len(data[data["is_parked"] != 1])
    date = data["DATETIME"].max()

    return html.Div(
        className="map-card-container",
        children=[
            html.Div(
                className="sip-card-stats-single",
                children=[
                    html.Div(
                        className="sip-card-title",
                        children=[
                            html.Img(src="assets/img/ship-2.svg"),
                            html.H3(
                                children=[
                                    f"{strings.SIP_CARD_MAP_IN_PORT}: {str(n_ships)}"
                                ]
                            ),
                        ],
                    ),
                    html.Div(
                        className="sip-card-body",
                        children=[
                            html.P(
                                children=[
                                    f"{strings.SIP_CARD_MAP_SAILING}: {str(n_ships_sailing)}, {strings.SIP_CARD_MAP_MOORED}: {str(n_ships_stop)}"
                                ]
                            )
                        ],
                    ),
                    html.Div(
                        className="sip-card-footer",
                        children=[
                            html.P(
                                children=[f"{strings.SIP_CARD_MAP_DATE}: {str(date)}"]
                            )
                        ],
                    ),
                ],
            )
        ],
    )
