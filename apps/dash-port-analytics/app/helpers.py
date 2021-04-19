import pandas as pd
from branca.element import Template, MacroElement
from config import constants, strings, styles


def filter_by_port_vessel_and_time(
    df: pd.DataFrame, port: str, vessel_type: str, year: int, month: int
) -> pd.DataFrame:
    """
    Repetitive function, used to perform initial filtering for data displayed on the dashboard.

    :param df:  Pandas DataFrame, input data
    :param port: str, port of interest
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: Pandas DataFrame, filtered
    """
    if vessel_type == strings.CITY_ALL:
        return df[
            (df["port"] == port) & (df["year"] == year) & (df["month"] == month)
        ].copy()
    return df[
        (df["port"] == port)
        & (df["ship_type"] == vessel_type)
        & (df["year"] == year)
        & (df["month"] == month)
    ].copy()


def filter_by_vessel_and_time(
    df: pd.DataFrame, vessel_type: str, year: int, month: int
) -> pd.DataFrame:
    """
    Repetitive function, used to perform initial filtering for data displayed on the dashboard.

    :param df: Pandas DataFrame, input data
    :param vessel_type: str, vessel type of interest
    :param year: int, year of interest
    :param month: int, month of interest
    :return: Pandas DataFrame, filtered
    """
    if vessel_type == strings.CITY_ALL:
        return df[(df["year"] == year) & (df["month"] == month)].copy()
    return df[
        (df["ship_type"] == vessel_type) & (df["year"] == year) & (df["month"] == month)
    ].copy()


def filter_by_vessel_and_port(
    df: pd.DataFrame, port1: str, port2: str, vessel_type: str
) -> pd.DataFrame:
    """
    Repetitive function, used to perform initial data filtering for comparison tab.

    :param df: Pandas DataFrame, input data
    :param port1: str, first port to compare
    :param port2: str, second port to compare
    :param vessel_type: str, vessel type of interest
    :return: Pandas DataFrame, filtered
    """
    if vessel_type == strings.CITY_ALL:
        return df[df["port"].isin([port1, port2])].copy()
    return df[
        (df["ship_type"] == vessel_type) & (df["port"].isin([port1, port2]))
    ].copy()


def get_dropdown_items(df: pd.DataFrame, attribute: str) -> list:
    """
    Returns a list of dropdown elements for a given attribute name.

    :param df: Pandas DataFrame object which contains the attribute
    :param attribute: str, can be either port, vessel_type, year, or month
    :return: list of unique attribute values
    """
    if attribute == "port":
        return df["port"].unique().tolist()
    elif attribute == "vessel_type":
        return ["All", *sorted(df["ship_type"].unique().tolist())]
    elif attribute == "year":
        return df["year"].unique().tolist()
    elif attribute == "month":
        return df["month"].unique().tolist()
    else:
        raise KeyError("Invalid value for `argument`")


def get_port_dropdown_values(curr_port_1: str, curr_port_2: str) -> list:
    """
    Returns list of possible port options for both dropdowns.

    :param curr_port_1: str, current value for the first port dropdown
    :param curr_port_2: str, current value for the second port dropdown
    :return: list, possible port values for both dropdowns
    """
    ports_dpd_1 = [
        strings.CITY_GDANSK,
        strings.CITY_GDYNIA,
        strings.CITY_KALINGRAD,
        strings.CITY_KLAIPEDA,
        strings.CITY_STPETERBURG,
    ]
    ports_dpd_2 = [
        strings.CITY_GDANSK,
        strings.CITY_GDYNIA,
        strings.CITY_KALINGRAD,
        strings.CITY_KLAIPEDA,
        strings.CITY_STPETERBURG,
    ]
    ports_dpd_1.remove(curr_port_2)
    ports_dpd_2.remove(curr_port_1)
    return [ports_dpd_1, ports_dpd_2]


def get_lat_long_for_port(port: str) -> list:
    """
    Returns a list of lat lon parameters for the map on the first tab.

    :param port: str, a location you want coordinates for
    :return: list - [latitude, longitude]
    """
    mappings = {
        strings.CITY_GDANSK: constants.LOC_GDANSK,
        strings.CITY_GDYNIA: constants.LOC_GDYNIA,
        strings.CITY_KALINGRAD: constants.LOC_KALINGRAD,
        strings.CITY_KLAIPEDA: constants.LOC_KLAIPEDA,
        strings.CITY_STPETERBURG: constants.LOC_STPETRSBURG,
    }
    return mappings[port]


def generate_map_legend() -> MacroElement:
    """
    Generates a legend for the map.

    :return: MacroElement, html added to the map
    """
    template = """
    {% macro html(this, kwargs) %}

    <!doctype html>
    <html lang="en">
    <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <title>Dashboard</title>
      <link rel="stylesheet" href="//code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">
      
      <link rel="preconnect" href="https://fonts.gstatic.com">
      <link href="https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap" rel="stylesheet">

      <script src="https://code.jquery.com/jquery-1.12.4.js"></script>
      <script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>


    </head>
    <body>


    <div id='maplegend' class='maplegend' 
      style='
        position: absolute; 
        z-index:9999; 
        background-color:rgba(255, 255, 255, 0.8);
        border: 1px solid #D4D4D4;
        border-radius:6px; 
        padding: 10px; 
        font-size:14px; 
        right: 10px; 
        bottom: 23px;
      '>

    <div class='legend-scale'>
      <ul class='legend-labels'>
        <li><span style='background:#E87272;'></span>Unspecified</li>
        <li><span style='background:#7a0091;'></span>Navigation</li>
        <li><span style='background:#11498A;'></span>Fishing</li>
        <li><span style='background:#1A6D9B;'></span>Tug</li>
        <li><span style='background:#12A5B0;'></span>Passenger</li>
        <li><span style='background:#3A9971;'></span>Cargo</li>
        <li><span style='background:#79BD00;'></span>Tanker</li>
        <li><span style='background:#DBB657;'></span>Pleasure</li>
      </ul>
    </div>
    </div>

    </body>
    </html>

    <style type='text/css'>
      * {
        font-family: "Roboto", sans-serif;
      }

      .maplegend .legend-title {
        text-align: left;
        margin-bottom: 5px;
        font-weight: bold;
        font-size: 90%;
        }
      .maplegend .legend-scale ul {
        margin: 0;
        margin-bottom: 5px;
        padding: 0;
        float: left;
        list-style: none;
        }
      .maplegend .legend-scale ul li {
        font-size: 80%;
        list-style: none;
        margin-left: 0;
        line-height: 18px;
        margin-bottom: 2px;
        }
      .maplegend ul.legend-labels li span {
        display: block;
        float: left;
        height: 16px;
        width: 30px;
        margin-right: 5px;
        margin-left: 0;
        border: 1px solid #999;
        }
      .maplegend .legend-source {
        font-size: 80%;
        color: #777;
        clear: both;
        }
      .maplegend a {
        color: #777;
        }

      .maplegend .legend-scale ul:last-child { 
        margin-bottom: 0px;
      }
      .maplegend .legend-scale ul li:last-child { 
        margin-bottom: 0px;
      }

    </style>
    {% endmacro %}"""

    macro = MacroElement()
    macro._template = Template(template)
    return macro


def generate_color(category: str) -> str:
    """
    Returns a color hex code for a given ship category.

    :param category: str, ship category (vessel type)
    :return: str, hex code for the color
    """
    mappings = {
        strings.STYPE_UNSPECIFIED: styles.COLOR_APPSILON_1,
        strings.STYPE_NAVIGATION: styles.COLOR_APPSILON_2,
        strings.STYPE_FISHING: styles.COLOR_APPSILON_3,
        strings.STYPE_TUG: styles.COLOR_APPSILON_4,
        strings.STYPE_PASSENGER: styles.COLOR_APPSILON_5,
        strings.STYPE_CARGO: styles.COLOR_APPSILON_6,
        strings.STYPE_TANKER: styles.COLOR_APPSILON_7,
        strings.STYPE_PLEASURE: styles.COLOR_APPSILON_8,
    }
    return mappings[category]
