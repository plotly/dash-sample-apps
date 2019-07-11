import visa

# Adapted from code seen here:
# https://github.com/baroobob/TektronixAFG3021B/blob/master/TektronixAFG3021B.py

rm = visa.ResourceManager()
fgenerator = None

def open_port(port = 'USB::0x0699::0x0340::C012268::INSTR'):
    global fgenerator
    try:
        fgenerator = rm.open_resource(port)

        # write("++mode 1")    # put Prologix in controller mode
        # write("++auto 0")    # turn off Prologix Read-After-Write mode
        # write("++addr 11")  # set GPIB address to the AFG3021B
        # write("*RST")        # Reset instrument
        device = fgenerator.query("*IDN?")      # ask instrument to identify itself
        write("++read 10")

        if not "TEKTRONIX,AFG3021" in device:
            print("Incompatible device")

    except visa.VisaIOError:
        print("ERROR: Unable to connect to AFG3021 function generator.")

def set_amplitude(amplitude):
    amplitude = isnumber(amplitude)
    if amplitude:
        offset = float(fgenerator.query("VOLTAGE:OFFSET?"))    # read present offset voltage

        if (amplitude < 10e-3):
            print('Warning: The minimum peak to peak amplitude for the AFG3021B '\
            'is 10 mV.')
        # # amplitude = 10e-3
        if (abs(amplitude/2) + abs(offset) > 5):
            print('Warning: The offset plus peak amplitude for the AFG3021B '\
            'cannot exceed +/-5 V.')
            amplitude = 2*(5 - abs(offset))

        write("VOLTAGE:AMPLITUDE " + str(amplitude))

def set_offset(offset):
    offset = isnumber(offset)
    if offset:
        write("VOLTAGE:AMPLITUDE?")      # read present amplitude

        if (abs(amplitude/2) + abs(offset) > 5):
            print('Warning: The offset plus peak amplitude for the AFG3021 '\
            'cannot exceed +/-5 V.')
        if (offset > 0):
            offset = 5 - amplitude/2
        else:
            offset = amplitude/2 - 5
        write("VOLTAGE:OFFSET " + str(offset))

def get_offset():
    return fgenerator.query("VOLTAGE:OFFSET?")

def get_frequency():
    # read present offset voltage
    return fgenerator.query("FREQUENCY?")

def get_amplitude():
    return fgenerator.query("VOLTAGE:AMPLITUDE?")

def set_frequency(frequency):
    write("FREQUENCY " + str(frequency))

def set_wave(wave):
    if wave in ['SIN', 'SQUARE', 'RAMP', 'PULSE']:
        write("FUNC " + wave)

# CHECK THIS - not tested
def get_wave(wave):
    return fgenerator.query("FUNC?")

def write(command):
    fgenerator.write(command)

def enable_output():
  write("OUTP ON")

def disable_output():
  write("OUTP OFF")

def toggle():
    on = int(get_output())
    if on:
        disable_output()
    else:
        enable_output()

def get_output():
  return fgenerator.query("OUTP?")

def isnumber(str):
    try:
        return float(str)
    except ValueError:
        return False
