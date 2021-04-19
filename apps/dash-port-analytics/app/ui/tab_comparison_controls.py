import dash_core_components as dcc
import dash_html_components as html
from config import strings


def make_port_comparison_controls(
    port_1_arr: list,
    port_1_val: str,
    port_2_arr: list,
    port_2_val: str,
    vessel_types_arr: list,
    vessel_type_val: str,
) -> html.Div:
    """
    Returns dropdown controls options for the Compare tab.

    :param port_1_arr: list, possible values for the first port
    :param port_1_val: str, current value for the first port
    :param port_2_arr: list, possible values for the second port
    :param port_2_val: str, current value for the second port
    :param vessel_types_arr: list, possible values for the vessel types
    :param vessel_type_val: str, current value for the vessel type
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
                        id="port-compare-port-1-dpd",
                        clearable=False,
                        options=[{"label": port, "value": port} for port in port_1_arr],
                        value=port_1_val,
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
                        id="port-compare-vessel-type-dpd",
                        clearable=False,
                        options=[
                            {"label": vessel, "value": vessel}
                            for vessel in vessel_types_arr
                        ],
                        value=vessel_type_val,
                    ),
                ],
            ),
            html.Div(className="tab-port-map-single-control-separator area-d"),
            html.Div(
                className="tab-port-map-single-control-container area-e",
                children=[
                    html.Label(
                        className="control-label", children=[strings.LABEL_PORT_COMPARE]
                    ),
                    dcc.Dropdown(
                        id="port-compare-port-2-dpd",
                        clearable=False,
                        options=[{"label": port, "value": port} for port in port_2_arr],
                        value=port_2_val,
                    ),
                ],
            ),
        ],
    )
