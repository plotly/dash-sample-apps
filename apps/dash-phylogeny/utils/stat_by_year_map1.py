import pandas as pd
from geopy.geocoders import Nominatim
import os
import re


def read_metadata(filename):
    df = pd.read_csv(filename)
    return df


def get_lon_lat(city):
    '''
    Example:
    location = geolocator.geocode("Chicago Illinois")
    return:
    Chicago, Cook County, Illinois, United States of America
    location.address    location.altitude   location.latitude   location.longitude  location.point      location.raw
    '''
    geolocator = Nominatim()
    location = geolocator.geocode(city, timeout=None)
    return location.longitude, location.latitude


def get_lon(city):
    '''
    Example:
    location = geolocator.geocode("Chicago Illinois")
    return:
    Chicago, Cook County, Illinois, United States of America
    location.address    location.altitude   location.latitude   location.longitude  location.point      location.raw
    '''
    geolocator = Nominatim()
    location = geolocator.geocode(city, timeout=1000)
    return location.longitude


def get_lat(city):
    '''
    Example:
    location = geolocator.geocode("Chicago Illinois")
    return:
    Chicago, Cook County, Illinois, United States of America
    location.address    location.altitude   location.latitude   location.longitude  location.point      location.raw
    '''
    geolocator = Nominatim()
    location = geolocator.geocode(city, timeout=1000)
    return location.latitude


def prepare_stat(path, filename, df, coordonnee_dic):
    lat = 0
    lon = 0
    ISO2 = 0
    ISO3 = 0
    data_metadata_stat_csv = (df.groupby(['Country', 'Year', 'Month'])).count()
    stat_file = open(path + "/stat_year_" + filename, "w")

    line = "Country,Month,Year,Lat,Lon,Value,ISO2,ISO3\n"
    stat_file.write(line)

    for index_val, series_val in data_metadata_stat_csv.iterrows():
        Country = index_val[0]
        Year = index_val[1]
        Month = index_val[2]
        Value = int(round(series_val[0], 0))

        if Country in coordonnee_dic:
            if len(coordonnee_dic[Country]) == 4:
                lat, lon, ISO2, ISO3 = coordonnee_dic[Country]
        else:
            print(Country + " not match.")
            #lat, lon = get_lon_lat(Country)
            #coordonnee_dic[Country] = (lat, lon)
        
        lat = round(lat, 2)
        lon = round(lon, 2)
        line = str(Country) + "," + str(int(Month)) + "," + str(Year) + "," + str(lat) + "," + str(lon) + "," + str(Value) + "," + str(ISO2) + "," + str(ISO3) + "\n"
        stat_file.write(line)
    stat_file.close()


def country_lon_lat_ISO(coordonnee_dic):
    """
    Convert daframe to dictionary
    :param path:
    :param filename of metada of each virus
    :return: dictionary
    """
    country_iso_list_long_lat2 = "data/country_iso_list_long_lat.csv"
    df_iso = read_metadata(country_iso_list_long_lat2)
    for index, row in df_iso.iterrows():
        coordonnee_dic[row[3]] = (row[1], row[2], row[0])

    country_iso_list_long_lat3 = "data/country_iso_list_long_lat_iso3letters.csv"
    df_iso = read_metadata(country_iso_list_long_lat3)
    for index, row in df_iso.iterrows():
        if row[1] in coordonnee_dic:
            coordonnee_dic[row[1]] = coordonnee_dic[row[1]] + (row[0],)
        else:
            print(row[1] + " did not find.")


def data_plus(path, filename):
    """
    Function to prepare correctly the date
    :param path:
    :param filename:
    :return: dataframe
    """
    metadata_file = path + "/" + filename
    df = read_metadata(metadata_file)
    df['Year'], df['Month'], df['Day'] = zip(*df['Date'].map(lambda x: x.split('-')))
    return df



folder_path = ".."
os.chdir(folder_path)
coordonnee_dic = {}

for path, dirs, files in os.walk(folder_path):
    for filename in files:
        p = re.search(r'[\w_]*metadata.csv', filename)
        if p:
            p = re.search(r'stat[\w_]*.csv', filename)
            if not p:
                print(path)
                print(filename)
                country_lon_lat_ISO(coordonnee_dic)
                df = data_plus(path, filename)
                prepare_stat(path, filename, df, coordonnee_dic)
