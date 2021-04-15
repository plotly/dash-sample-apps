import dash_core_components as dcc
import dash_html_components as html
from config import strings


def make_tab_port_map_controls(
    port_arr: list,
    port_val: str,
    vessel_types_arr: list,
    vessel_type_val: str,
    year_arr: list,
    year_val: int,
    month_arr: list,
    month_val: int,
) -> html.Div:
    """
    Returns a HTML div of user controls found on top of the map tab.

    :param port_arr: list, all possible ports
    :param port_val: str, current port value
    :param vessel_types_arr: list, all possible vessel types
    :param vessel_type_val: str, current vessel type value
    :param year_arr: list, all possible years
    :param year_val: str, current year value
    :param month_arr: list, all possible months
    :param month_val: str, current month value
    :return: HTML div
    """
    return html.Div(
        className="tab-port-map-controls",
        children=[
            html.Div(
                className="tab-port-map-single-control-container area-a",
                children=[
                    html.Label(
                        className="control-label", children=[strings.LABEL_PORT]
                    ),
                    dcc.Dropdown(
                        id="port-map-dropdown-port",
                        clearable=False,
                        options=[{"label": port, "value": port} for port in port_arr],
                        value=port_val,
                    ),
                ],
            ),
            html.Div(className="tab-port-map-single-control-separator area-b"),
            html.Div(
                className="tab-port-map-single-control-container area-c",
                children=[
                    html.Label(
                        className="control-label", children=[strings.LABEL_VESSEL]
                    ),
                    dcc.Dropdown(
                        id="port-map-dropdown-vessel-type",
                        clearable=False,
                        options=[
                            {"label": vessel_type, "value": vessel_type}
                            for vessel_type in vessel_types_arr
                        ],
                        value=vessel_type_val,
                    ),
                ],
            ),
            html.Div(className="tab-port-map-single-control-separator area-d"),
            html.Div(
                className="tab-port-map-single-control-container date-grid area-e",
                children=[
                    html.Div(
                        className="tab-port-map-single-control-container-date",
                        children=[
                            html.Label(
                                className="control-label", children=[strings.LABEL_YEAR]
                            ),
                            dcc.Dropdown(
                                id="port-map-dropdown-year",
                                clearable=False,
                                options=[
                                    {"label": year, "value": year} for year in year_arr
                                ],
                                value=year_val,
                            ),
                        ],
                    ),
                    html.Div(
                        className="tab-port-map-single-control-separator smaller-line"
                    ),
                    html.Div(
                        className="tab-port-map-single-control-container-date",
                        children=[
                            html.Label(
                                className="control-label",
                                children=[strings.LABEL_MONTH],
                            ),
                            dcc.Dropdown(
                                id="port-map-dropdown-month",
                                clearable=False,
                                options=[
                                    {"label": month, "value": month}
                                    for month in month_arr
                                ],
                                value=month_val,
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )
