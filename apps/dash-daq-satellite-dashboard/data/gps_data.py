from random import randint
import math

#######################################################################################################################
# Setup
#######################################################################################################################
# Default data
file_m = open('gps_data_m.csv', 'w')
file_h = open('gps_data_h.csv', 'w')

# Satellite H45-k1 data
file_m_0 = open('gps_data_m_0.csv', 'w')
file_h_0 = open('gps_data_h_0.csv', 'w')

# Satellite L12-5 data
file_m_1 = open('gps_data_m_1.csv', 'w')
file_h_1 = open('gps_data_h_1.csv', 'w')

#######################################################################################################################
# Data
#######################################################################################################################
file_m.write('lat,lon\n')
file_h.write('lat,lon\n')

file_m_0.write('lat,lon\n')
file_h_0.write('lat,lon\n')

file_m_1.write('lat,lon\n')
file_h_1.write('lat,lon\n')

for i in range(3600):
    file_m.write('0,0\n')

for i in range(60):
    file_h.write('0,0\n')


for i in range(3600):
    lat = 60 * math.cos(i * 2 * math.pi / 3600)
    lon = 360/3600 * i
    file_m_0.write('%f,%f\n' % (lat, lon))

for i in range(60):
    lat = 60 * math.cos(i * 2 * math.pi / 60)
    lon = 360/60 * i
    file_h_0.write('%f,%f\n' % (lat, lon))

for i in range(3600):
    lat = 20 + 50 * math.cos(i * 2 * math.pi / 3600 - math.pi/3)
    lon = 360/3600 * i
    file_m_1.write('%f,%f\n' % (lat, lon))

for i in range(60):
    lat = 20 + 50 * math.cos(i * 2 * math.pi / 60 - math.pi/3)
    lon = 360/60 * i
    file_h_1.write('%f,%f\n' % (lat, lon))
