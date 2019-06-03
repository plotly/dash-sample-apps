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
    location = geolocator.geocode(city, timeout=3)
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


def create_fig(path, filename):
    metadata_file = path + "/" + filename
    df = read_metadata(metadata_file)
    data_metadata_stat_csv = df.groupby('Country')['Strain'].count()

    statFile = open(path + "/stat_" + filename, "w")

    line = "name,pop,lat,lon\n"
    # print(line)
    statFile.write(line)

    for index_val, series_val in data_metadata_stat_csv.iteritems():
        lon, lat = get_lon_lat(index_val)
        # line = str(index_val) + "," + str(series_val) + "," + str(get_lat(index_val)) + "," + str(get_lon(index_val)) + "\n"
        line = str(index_val) + "," + str(series_val) + "," + str(lat) + "," + str(lon) + "\n"
        # print(line)
        statFile.write(line)

    statFile.close()


folder_path = "."

for path, dirs, files in os.walk(folder_path):
    for filename in files:
        p = re.search(r'[\w_]*metadata.csv', filename)
        if p:
            print(path)
            print(filename)
            # print(p)
            # os.chdir(path)
            create_fig(path, filename)
