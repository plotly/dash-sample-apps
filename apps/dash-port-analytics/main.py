import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
from app import helpers
from app.ui import (
    header,
    contact_modal,
    tab_comparison_controls,
    tab_comparison_sip_cards,
    tab_port_sip_cards,
    tab_map_controls,
    tab_map_sip_cards,
)
from app import tab_map, tab_stats, tab_compare
from config import strings, constants
from dash.dependencies import Input, Output


# DATASET LOADING AND TRANSFORMATION
df_tab_1 = pd.read_csv("data/first_tab_dataset.csv")
df_ships_duration = pd.read_csv("data/ships_stop_duration.csv")
df_ships_duration["min_time"] = pd.to_datetime(df_ships_duration["min_time"])
df_ships_duration["year"] = df_ships_duration["min_time"].dt.year
df_ships_duration["month"] = df_ships_duration["min_time"].dt.month
df_tab_1 = pd.merge(left=df_tab_1, right=df_ships_duration, on="SHIP_ID")
df_tab_1 = df_tab_1.drop(["port_y", "ship_type_y"], axis=1)
df_tab_1 = df_tab_1.rename({"port_x": "port", "ship_type_x": "ship_type"}, axis=1)
df_tab_1["year"] = df_tab_1["date"].astype("datetime64").dt.year
df_tab_1["month"] = df_tab_1["date"].astype("datetime64").dt.month
df_tab_2 = pd.read_csv("data/second_tab_dataset.csv")
df_tab_2["date"] = pd.to_datetime(df_tab_2["date"])
df_tab_2["year"] = df_tab_2["date"].dt.year
df_tab_2["month"] = df_tab_2["date"].dt.month
df_tab_3 = pd.read_csv("data/third_tab_dataset.csv", index_col=[0])

# EXTERNAL SCRIPTS AND STYLES
external_scripts = ["https://kit.fontawesome.com/0bb0d79500.js"]
external_stylesheets = [
    "https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css"
]

app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    external_scripts=external_scripts,
    external_stylesheets=external_stylesheets,
)
app.title = strings.APP_NAME
app.config["suppress_callback_exceptions"] = True
server = app.server

# MAP TAB OPTIONS
dpd_options_port = helpers.get_dropdown_items(df=df_tab_1, attribute="port")
dpd_options_vessel = helpers.get_dropdown_items(df=df_tab_1, attribute="vessel_type")
dpd_options_year = helpers.get_dropdown_items(df=df_tab_1, attribute="year")
dpd_options_month = helpers.get_dropdown_items(df=df_tab_1, attribute="month")

# TAB COMPARE OPTIONS
dpd_options_port1, dpd_options_port2 = helpers.get_port_dropdown_values(
    curr_port_1=strings.CITY_GDANSK, curr_port_2=strings.CITY_GDYNIA
)

curr_port = strings.CITY_GDANSK
curr_vessel = strings.CITY_ALL
curr_year = 2016
curr_month = 12


# GENERAL LAYOUT
app.layout = html.Div(
    [
        header.make_header(),
        html.Div(
            className="wrapper",
            children=[
                html.Div(id="main-area", className="main-area"),
                contact_modal.make_contact_modal(),
            ],
        ),
        html.Footer(
            className="logo-footer",
            children=[
                html.Footer(
                    className="logo-footer-centering",
                    children=[
                        html.A(
                            href="https://appsilon.com/",
                            target="_blank",
                            className="logo-footer-container",
                            children=[
                                html.H4(
                                    children=["Powered by"], className="footer-element"
                                ),
                                html.Img(
                                    src="assets/img/appsilon-icon.svg",
                                    id="appsilon-icon",
                                    className="footer-element",
                                ),
                                html.Img(
                                    src="assets/img/appsilon-text.svg",
                                    id="appsilon-text",
                                    className="footer-element",
                                ),
                            ],
                        ),
                    ],
                ),
            ],
        ),
    ]
)


# TAB RENDERER
@app.callback(Output("main-area", "children"), [Input("navigation-tabs", "value")])
def render_tab(tab):
    """
    Renders content depending on the tab selected.

    :param tab: tab option selected
    :return: HTML div
    """
    if tab == "tab-port-map":
        return [
            html.Div(
                id="tab-port-map-container",
                className="tab-port-map-container",
                children=[
                    tab_map_controls.make_tab_port_map_controls(
                        port_arr=dpd_options_port,
                        port_val=curr_port,
                        vessel_types_arr=dpd_options_vessel,
                        vessel_type_val=curr_vessel,
                        year_arr=dpd_options_year,
                        year_val=curr_year,
                        month_arr=dpd_options_month,
                        month_val=curr_month,
                    )
                ],
            )
        ]
    elif tab == "tab-port-stats":
        return [
            html.Div(
                id="tab-port-stats-container",
                className="tab-port-stats-container",
                children=[
                    tab_map_controls.make_tab_port_map_controls(
                        port_arr=dpd_options_port,
                        port_val=curr_port,
                        vessel_types_arr=dpd_options_vessel,
                        vessel_type_val=curr_vessel,
                        year_arr=dpd_options_year,
                        year_val=curr_year,
                        month_arr=dpd_options_month,
                        month_val=curr_month,
                    )
                ],
            )
        ]
    elif tab == "tab-port-compare":
        return [
            html.Div(
                id="tab-port-compare-container",
                className="tab-port-compare-container",
                children=[
                    tab_comparison_controls.make_port_comparison_controls(
                        port_1_arr=dpd_options_port1,
                        port_1_val=dpd_options_port1[0],
                        port_2_arr=dpd_options_port2,
                        port_2_val=dpd_options_port2[0],
                        vessel_types_arr=dpd_options_vessel,
                        vessel_type_val=dpd_options_vessel[0],
                    )
                ],
            )
        ]


# MAP RENDERER (TAB 1)
@app.callback(
    Output("tab-port-map-container", "children"),
    [
        Input("port-map-dropdown-port", "value"),
        Input("port-map-dropdown-vessel-type", "value"),
        Input("port-map-dropdown-year", "value"),
        Input("port-map-dropdown-month", "value"),
    ],
)
def update_port_map_tab(port, vessel_type, year, month):
    """
    Renders content for the Map tab.

    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: HTML div
    """
    global curr_port
    global curr_vessel
    global curr_year
    global curr_month
    curr_port = port
    curr_vessel = vessel_type
    curr_year = year
    curr_month = month

    return html.Div(
        children=[
            tab_map_controls.make_tab_port_map_controls(
                port_arr=dpd_options_port,
                port_val=port,
                vessel_types_arr=dpd_options_vessel,
                vessel_type_val=vessel_type,
                year_arr=dpd_options_year,
                year_val=year,
                month_arr=dpd_options_month,
                month_val=month,
            ),
            tab_map_sip_cards.make_tab_port_map_sip_card(
                df=df_tab_1, port=port, vessel_type=vessel_type, year=year, month=month
            ),
            tab_map.make_tab_port_map_map(
                df=df_tab_1, port=port, vessel_type=vessel_type, year=year, month=month
            ),
            tab_map.make_tab_port_map_table(
                df=df_tab_1, port=port, vessel_type=vessel_type, year=year, month=month
            ),
        ]
    )


# STATS RENDERER (TAB 2)
@app.callback(
    Output("tab-port-stats-container", "children"),
    [
        Input("port-map-dropdown-port", "value"),
        Input("port-map-dropdown-vessel-type", "value"),
        Input("port-map-dropdown-year", "value"),
        Input("port-map-dropdown-month", "value"),
    ],
)
def update_port_stats_tab(port, vessel_type, year, month):
    """
    Renders content for the Stats tab.

    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: HTML div
    """
    global curr_port
    global curr_vessel
    global curr_year
    global curr_month
    curr_port = port
    curr_vessel = vessel_type
    curr_year = year
    curr_month = month

    return html.Div(
        children=[
            tab_map_controls.make_tab_port_map_controls(
                port_arr=dpd_options_port,
                port_val=port,
                vessel_types_arr=dpd_options_vessel,
                vessel_type_val=vessel_type,
                year_arr=dpd_options_year,
                year_val=year,
                month_arr=dpd_options_month,
                month_val=month,
            ),
            tab_port_sip_cards.make_tab_port_stats_sip_cards(
                df=df_tab_2,
                df_stop=df_ships_duration,
                port=port,
                vessel_type=vessel_type,
                year=year,
                month=month,
            ),
            html.Div(
                className="graphGrid",
                children=[
                    html.Div(
                        className="ui-card",
                        children=[
                            html.H3(strings.CHART_STATS_TOTAL_VESSELS_TITLE),
                            dcc.Graph(
                                config={"displayModeBar": False},
                                figure=tab_stats.plot_stats_total_num_vessels(
                                    df=df_tab_2,
                                    port=port,
                                    vessel_type=vessel_type,
                                    year=year,
                                    month=month,
                                ),
                                style={"height": "370px"},
                            ),
                        ],
                    ),
                    html.Div(
                        className="ui-card",
                        children=[
                            html.H3(strings.CHART_STATS_STOP_DUR_TITLE),
                            dcc.Graph(
                                config={"displayModeBar": False},
                                figure=tab_stats.plot_avg_vessel_stop_duration(
                                    df=df_ships_duration,
                                    port=port,
                                    vessel_type=vessel_type,
                                    year=year,
                                    month=month,
                                ),
                                style={"height": "370px"},
                            ),
                        ],
                    ),
                    html.Div(
                        className="ui-card",
                        children=[
                            html.H3(strings.CHART_STATS_TOTAL_CAP_VESSELS_TITLE),
                            dcc.Graph(
                                config={"displayModeBar": False},
                                figure=tab_stats.plot_total_capacity_of_vessels(
                                    df=df_tab_2,
                                    port=port,
                                    vessel_type=vessel_type,
                                    year=year,
                                    month=month,
                                ),
                                style={"height": "370px"},
                            ),
                        ],
                    ),
                ],
            ),
        ]
    )


# COMPARE RENDERER (TAB 3)
@app.callback(
    Output("tab-port-compare-container", "children"),
    [
        Input("port-compare-port-1-dpd", "value"),
        Input("port-compare-port-2-dpd", "value"),
        Input("port-compare-vessel-type-dpd", "value"),
    ],
)
def update_compare_tab(port1, port2, vessel_type):
    """
    Renders content for the Stats tab.

    :param port1: str, a port to compare
    :param port2:  str, a port to compare
    :param vessel_type: str, vessel type of interest
    :return: HTML div
    """
    port_comparison = tab_compare.get_compare_tab_insights(
        df=df_tab_3, port1=port1, port2=port2, vessel_type=vessel_type
    )
    return html.Div(
        children=[
            tab_comparison_controls.make_port_comparison_controls(
                port_1_arr=dpd_options_port1,
                port_1_val=port1,
                port_2_arr=dpd_options_port2,
                port_2_val=port2,
                vessel_types_arr=dpd_options_vessel,
                vessel_type_val=vessel_type,
            ),
            html.Div(
                className="graphGrid",
                children=[
                    html.Div(
                        className="ui-card",
                        children=[
                            html.H3(strings.CHART_COMPARE_NUM_VESSELS_TITLE),
                            dcc.Graph(
                                config={"displayModeBar": False},
                                figure=tab_compare.plot_num_vessels_comparison(
                                    df=df_tab_3,
                                    port1=port1,
                                    port2=port2,
                                    vessel_type=vessel_type,
                                ),
                                style={"height": "370px"},
                            ),
                        ],
                    ),
                    html.Div(
                        className="ui-card",
                        children=[
                            html.H3(strings.CHART_COMPARE_AVG_DURATION_TITLE),
                            dcc.Graph(
                                config={"displayModeBar": False},
                                figure=tab_compare.plot_avg_stop_duration_comparison(
                                    df=df_ships_duration,
                                    port1=port1,
                                    port2=port2,
                                    vessel_type=vessel_type,
                                ),
                                style={"height": "370px"},
                            ),
                        ],
                    ),
                    html.Div(
                        className="ui-card",
                        children=[
                            html.H3(strings.CHART_COMPARE_CAPACITY_TITLE),
                            dcc.Graph(
                                config={"displayModeBar": False},
                                figure=tab_compare.plot_avg_sum_capacity_comparison(
                                    df=df_tab_3,
                                    port1=port1,
                                    port2=port2,
                                    vessel_type=vessel_type,
                                ),
                                style={"height": "370px"},
                            ),
                        ],
                    ),
                ],
            ),
            tab_comparison_sip_cards.make_port_comparison_sip_cards(
                port1=port1,
                port1_inc=port_comparison["Increase_P1"],
                port1_dec=port_comparison["Decrease_P1"],
                port1_trend=port_comparison["Trend_P1"],
                port2=port2,
                port2_inc=port_comparison["Increase_P2"],
                port2_dec=port_comparison["Decrease_P2"],
                port2_trend=port_comparison["Trend_P2"],
            ),
        ]
    )


# COMPARE TAB DROPDOWNS
@app.callback(
    Output("port-compare-port-1-dpd", "options"),
    [Input("port-compare-port-2-dpd", "value")],
)
def update_options_dpd1(dpd2_val) -> list:
    """
    Updates the contents of the first dropdown menu based of the value of the second dropdown.

    :param dpd2_val: str, second dropdown value
    :return: list of dictionaries, labels and values
    """
    all_options = [
        strings.CITY_GDANSK,
        strings.CITY_GDYNIA,
        strings.CITY_KALINGRAD,
        strings.CITY_KLAIPEDA,
        strings.CITY_STPETERBURG,
    ]
    all_options.remove(dpd2_val)
    options = [{"label": opt, "value": opt} for opt in all_options]
    return options


@app.callback(
    Output("port-compare-port-2-dpd", "options"),
    [Input("port-compare-port-1-dpd", "value")],
)
def update_options_dpd2(dpd1_val):
    """
    Updates the contents of the second dropdown menu based of the value of the first dropdown.

    :param dpd1_val: str, first dropdown value
    :return: list of dictionaries, labels and values
    """
    all_options = [
        strings.CITY_GDANSK,
        strings.CITY_GDYNIA,
        strings.CITY_KALINGRAD,
        strings.CITY_KLAIPEDA,
        strings.CITY_STPETERBURG,
    ]
    all_options.remove(dpd1_val)
    options = [{"label": opt, "value": opt} for opt in all_options]
    return options


if __name__ == "__main__":
    # app.run_server(host='0.0.0.0', port=9000) # production
    app.run_server(host="0.0.0.0", port=9000, debug=True)  # development
    # app.run_server()
