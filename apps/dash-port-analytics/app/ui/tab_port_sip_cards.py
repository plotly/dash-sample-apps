import pandas as pd
import dash_html_components as html
from config import strings
from app import tab_stats


def make_tab_port_stats_sip_cards(
    df: pd.DataFrame,
    df_stop: pd.DataFrame,
    port: str,
    vessel_type: str,
    year: int,
    month: int,
) -> html.Div:
    """
    Returns HTML code for the cards in the Stats tab.

    :param df: Pandas DataFrame, input data
    :param df_stop: Pandas DataFrame, dataset with stop durations
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: HTML div
    """

    def make_single_sip_card(
        heading: str, pct: int, direction: str, port: str
    ) -> html.Div:
        """
        Returns HTML code for a single card in the Stats tab.

        :param heading: str, card heading text
        :param pct: int, percentage
        :param direction: str, lower or higher
        :param port: str, port of interest
        :return: HTML div
        """
        return html.Div(
            className="sip-card-stats-single",
            children=[
                html.Div(
                    className="sip-card-title",
                    children=[
                        html.Img(src="assets/img/marker.svg"),
                        html.H3(children=[f"{heading}: {pct}%"]),
                    ],
                ),
                html.Div(
                    className="sip-card-body",
                    children=[
                        html.P(children=[f"{direction} {strings.TAB1_CARD_BODY_TEXT}"])
                    ],
                ),
                html.Div(
                    className="sip-card-footer",
                    children=[
                        html.P(children=[f"{strings.TAB1_CARD_FOOTER_TEXT} {port}"])
                    ],
                ),
            ],
        )

    data_card_1 = tab_stats.get_stats_card1_data(
        df=df, port=port, vessel_type=vessel_type, year=year, month=month
    )
    data_card_2 = tab_stats.get_stats_card2_data(
        df=df_stop, port=port, vessel_type=vessel_type, year=year, month=month
    )
    data_card_3 = tab_stats.get_stats_card3_data(
        df=df, port=port, vessel_type=vessel_type, year=year, month=month
    )

    return html.Div(
        className="stats-card-container",
        children=[
            make_single_sip_card(
                heading=strings.TAB1_CARD1_HEADING,
                pct=data_card_1[0],
                direction=data_card_1[1],
                port=port,
            ),
            make_single_sip_card(
                heading=strings.TAB1_CARD2_HEADING,
                pct=data_card_2[0],
                direction=data_card_2[1],
                port=port,
            ),
            make_single_sip_card(
                heading=strings.TAB1_CARD3_HEADING,
                pct=data_card_3[0],
                direction=data_card_3[1],
                port=port,
            ),
        ],
    )
