import dash_core_components as dcc
import dash_html_components as html
from config import strings


def make_header() -> html.Header:
    """
    Returns a HTML Header element for the application Header.

    :return: HTML Header
    """
    return html.Header(
        children=[
            # Icon and title container
            html.Div(
                className="dash-title-container",
                children=[
                    html.Img(className="dash-icon", src="assets/img/ship-1.svg"),
                    html.H1(className="dash-title", children=["Dash ports analytics"]),
                ],
            ),
            # create navigator with buttons
            html.Nav(
                children=[
                    dcc.Tabs(
                        id="navigation-tabs",
                        value="tab-port-map",
                        children=[
                            dcc.Tab(
                                label=strings.TAB1_NAME,
                                value="tab-port-map",
                                className="dash-tab",
                                selected_className="dash-tab-selected",
                            ),
                            dcc.Tab(
                                label=strings.TAB2_NAME,
                                value="tab-port-stats",
                                className="dash-tab",
                                selected_className="dash-tab-selected",
                            ),
                            dcc.Tab(
                                label=strings.TAB3_NAME,
                                value="tab-port-compare",
                                className="dash-tab",
                                selected_className="dash-tab-selected",
                            ),
                        ],
                    ),
                    # TODO Dario - remove below button
                    # html.Button(id='btn-sidebar-request-port', className='btn-sidebar-request-port', children=[strings.BTN_REQUEST_PORT])
                ]
            ),
        ]
    )
