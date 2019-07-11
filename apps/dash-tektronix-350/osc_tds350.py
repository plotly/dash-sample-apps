import visa
import numpy as np

oscilloscope = None

# Adapted from code seen here:
# https://github.com/baroobob/TektronixTDS2024B/blob/master/TektronixTDS2024B.py

def get_data():
    global oscilloscope
    rm = visa.ResourceManager()
    oscilloscope = rm.open_resource("GPIB0::1::INSTR")

    write("DATA:SOURCE CH1")
    write("DATA:WIDTH 2")
    write("DATa:ENCdg SRIbinary")

    ymult = float(query("WFMPRE:CH1:YMULT?"))
    yzero = float(query("WFMPRE:CH1:YZERO?"))
    yoff = float(query('WFMPRE:CH1:YOFF?'))
    xincr = float(query('WFMPRE:CH1:XINCR?'))

    write('AUTOSET EXECUTE')

    write("CURVE?")
    data = oscilloscope.read_raw()
    headerlen = 2 + int(data[1])
    header = data[:headerlen]
    ADC_wave = data[headerlen:-1]
    ADC_wave = np.fromstring(ADC_wave, dtype = np.int16)

    y = (ADC_wave - yoff) * ymult  + yzero
    x = np.arange(0, xincr * len(y), xincr)

    oscilloscope.close()

    return [{'x': x,
             'y': y,
             'type': 'line',
             'showscale': False,
             'colorscale': [[0, 'rgba(255, 255, 255,0)'], [1, 'rgba(0,0,255,1)']]}]

def get_data_tuple():
    global oscilloscope
    rm = visa.ResourceManager()
    oscilloscope = rm.open_resource("GPIB0::1::INSTR")

    write("DATA:SOURCE CH1")
    write("DATA:WIDTH 2")
    write("DATa:ENCdg SRIbinary")

    ymult = float(query("WFMPRE:CH1:YMULT?"))
    yzero = float(query("WFMPRE:CH1:YZERO?"))
    yoff = float(query('WFMPRE:CH1:YOFF?'))
    xincr = float(query('WFMPRE:CH1:XINCR?'))

    write('AUTOSET EXECUTE')

    write("CURVE?")
    data = oscilloscope.read_raw()
    headerlen = 2 + int(data[1])
    header = data[:headerlen]
    ADC_wave = data[headerlen:-1]
    ADC_wave = np.fromstring(ADC_wave, dtype = np.int16)

    y = (ADC_wave - yoff) * ymult  + yzero
    x = np.arange(0, xincr * len(y), xincr)
    return (x, y)

def query(command):
    return oscilloscope.query(command)

def write(command):
    oscilloscope.write(command)
#page 204
