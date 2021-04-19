import folium
import dash_table
import numpy as np
import pandas as pd
import dash_html_components as html
from app import helpers
from config import strings


def make_tab_port_map_map(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
) -> html.Div:
    """
    Makes the interactive map for the Map tab.

    :param df: Pandas DataFrame, input data
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: HTML div with the embedded map
    """

    def generate_popup(row: pd.Series) -> folium.Popup:
        """
        Makes a popup for a single map marker.

        :param row: Pandas Series, a single row of the dataset
        :return: Folium Popup object
        """
        html = f"""<div class="map-popup">
                        <b>{strings.MAP_POPUP_VESSEL}</b> {row.SHIPNAME}<br>
                        <b>{strings.MAP_POPUP_LATITUDE}</b>{np.round(row.LAT, 5)}<br>
                        <b>{strings.MAP_POPUP_LONGITUDE}</b>{np.round(row.LON, 5)}<br>
                        <b>{strings.MAP_POPUP_STOP_DUR}</b>{row.len_stop}
                    </div>"""
        iframe = folium.IFrame(html, width=300, height=90)
        return folium.Popup(iframe)

    def generate_center_coordinates(df: pd.DataFrame) -> list:
        """
        Returns a list of center coordinates for the filtered dataset

        :param df: Pandas DataFrame, input data
        :return: list - [latitude, longitude]
        """
        if len(df) > 0:
            return [df["center_lat"].median(), df["center_lon"].median()]
        return [-1, -1]

    data = helpers.filter_by_port_vessel_and_time(
        df=df, port=port, vessel_type=vessel_type, year=year, month=month
    )
    geolocation = generate_center_coordinates(df=data)
    if geolocation == [-1, -1]:
        geolocation = helpers.get_lat_long_for_port(port=port)
    port_map = folium.Map(
        location=geolocation, zoom_start=11.5, tiles="CartoDB positron"
    )
    for row in data.itertuples(index=False):
        folium.Circle(
            location=[row.LAT, row.LON],
            popup=generate_popup(row),
            color=helpers.generate_color(row.ship_type),
            weight=20,
            opacity=0.85,
        ).add_to(port_map)
        map_legend = helpers.generate_map_legend()
        port_map.get_root().add_child(map_legend)

    port_map.save("data/index.html")
    return html.Div(
        className="map-container",
        children=[html.Iframe(srcDoc=open("data/index.html", "r").read())],
    )


def make_tab_port_map_table(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
) -> html.Div:
    """
    Makes a table shown below the map on the Map tab.

    :param df: Pandas DataFrame, input data
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: HTML div with embedded table
    """
    data = helpers.filter_by_port_vessel_and_time(
        df=df, port=port, vessel_type=vessel_type, year=year, month=month
    )
    columns = [
        "SHIPNAME",
        "ship_type",
        "DATETIME",
        "LAT",
        "LON",
        "DWT",
        "WIDTH",
        "LENGTH",
        "SPEED",
    ]
    data = data[columns]
    data.columns = [
        "Name",
        "Type",
        "Date",
        "Latitude",
        "Longitude",
        "Deadweight",
        "Width",
        "Length",
        "Speed",
    ]
    data["Latitude"] = data["Latitude"].apply(lambda x: np.round(x, 5))
    data["Longitude"] = data["Longitude"].apply(lambda x: np.round(x, 5))

    return html.Div(
        className="map-table-container",
        children=[
            dash_table.DataTable(
                id="table",
                columns=[{"name": i, "id": i} for i in data.columns],
                data=data.to_dict("records"),
                page_size=10,
                sort_action="native",
                filter_action="native",
                style_cell={"padding": "15px 5px", "boxShadow": "0 0",},
                style_data={"border": "0px", "textAlign": "center"},
                style_header={
                    "padding": "2px 5px",
                    "fontWeight": "bold",
                    "textAlign": "center",
                    "border": "none",
                    "backgroundColor": "transparent",
                },
                style_table={"overflowX": "auto", "width": "calc(100% - 26px)",},
                style_data_conditional=[
                    {
                        "if": {"state": "selected"},
                        "backgroundColor": "transparent",
                        "border": "0px solid transparent",
                    }
                ],
            )
        ],
    )
