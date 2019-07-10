# -*- coding: utf-8 -*-
"""
Created on Tue Jan 05 11:53:50 2018
Useful fonctions and classes to connect to instruments in general using GPIB,
COM port or Prologix controller

"""
import logging

import glob
import sys

try:
    import serial
    serial_available = True
except ImportError:
    serial_available = False
    print("pyserial package not installed, run 'pip install pyserial'")

try:
    import visa
    visa_available = True
except ImportError:
    visa_available = False
    print("pyvisa package not installed, run 'pip install pyvisa'")

INTF_VISA = 'pyvisa'
INTF_PROLOGIX = 'prologix'
INTF_GPIB = INTF_PROLOGIX
INTF_SERIAL = 'serial'

INTF_NONE = 'None'

PROLOGIX_COM_PORT = "COM3"


def list_gpib_ports():
    """ use pyvisa to list the GPIB ports """

    rm = visa.ResourceManager()
    available_ports = rm.list_resources()
    temp_ports = []
    for port in available_ports:
        if "GPIB" in port:
            temp_ports.append(str(port))
    available_ports = temp_ports

    return available_ports


def list_serial_ports(max_port_num=20):
    """ Lists serial port names from COM1 to COM20 in windows platform and
    lists all serial ports on linux platform

        :raises EnvironmentError:
            On unsupported or unknown platforms
        :returns:
            A list of the serial ports available on the system
    """
    if sys.platform.startswith('win'):
        ports = ['COM%s' % (i + 1) for i in range(max_port_num)]
    elif sys.platform.startswith('linux') or sys.platform.startswith('cygwin'):
        # this excludes your current terminal "/dev/tty"
        ports = glob.glob('/dev/tty[A-Za-z]*')
    elif sys.platform.startswith('darwin'):
        ports = glob.glob('/dev/tty.*')
    else:
        raise EnvironmentError('Unsupported platform')

    result = []
    for port in ports:
        try:
            s = serial.Serial(port, 9600, timeout=0.1)
            result.append(str(port))
            s.close()
        except (OSError):
            pass
    return result


def refresh_device_port_list(debug=False):
    """ Load VISA resource list for use in combo boxes """

    if debug:

        return ["GPIB0::%i" % i for i in range(30)]

    else:

        return list_gpib_ports() + list_serial_ports()


def find_prologix_ports():
    """
        go through the serial ports and wee which one returns the prologix
        version command
    """
    serial_ports = list_serial_ports()
    result = []
    for port in serial_ports:
        try:
            s = serial.Serial(port, 9600, timeout=0.2)
            s.write("++mode 1\n++auto 0\n++ver\n".encode())
            answer = s.readline()
            prologix_controller = "Prologix GPIB-USB Controller"
            prologix_controller = prologix_controller.encode()
            if prologix_controller in answer:
                result.append(port)
        except(OSError, serial.SerialException):
            pass

    return result


def test_prologix_controller_creation_with_com(com_port=None):
    if com_port is None:
        com_port = "COM3"

    pc = PrologixController(com_port)
    print(pc)
    print(pc.get_open_gpib_ports())


def test_prologix_controller_creation_with_wrong_com():
    pc = PrologixController("COM1")
    print(pc)
    print(pc.get_open_gpib_ports())


def test_prologix_controller_creation_with_no_arg_conflict():
    pc = PrologixController()
    print(pc)
    pc2 = PrologixController()
    print(pc2)


class PrologixController(object):
    connection = None

    def __init__(
            self,
            com_port=None,
            mock=False,
            auto=1,
            baud_rate=9600,
            timeout=5,
            **kwargs
    ):

        self.mock = mock

        if not self.mock:
            if com_port is None:
                # the user didn't provide a COM port, so we look for one
                com_port = find_prologix_ports()

                if com_port != []:
                    if len(com_port) > 1:
                        logging.warning("There is more than one Prologix \
                         controller, we are connecting to %s" % (com_port[0]))

                    com_port = com_port[0]
                    print(
                        "... found a Prologix controller on the port '%s'" %
                        com_port
                    )

                    self.connection = serial.Serial(
                        com_port,
                        baud_rate,
                        xonxoff=True,
                        stopbits=serial.STOPBITS_TWO,
                        timeout=timeout
                    )
                else:
                    self.connection = None
                    print("There is no Prologix controller to connect to")
            else:

                try:
                    self.connection = serial.Serial(
                        com_port,
                        baud_rate,
                        xonxoff=True,
                        stopbits=serial.STOPBITS_TWO,
                        timeout=timeout
                    )
                except serial.serialutil.SerialException:
                    self.connection = None
                    print(
                        "The port %s is not attributed to any device"
                        % com_port
                    )

            if self.connection is not None:
                # set the connector in controller mode and let the user
                self.write("++mode 1")
                # auto == 1 : ask for read without sending another command.
                # auto == 0 : simply send the command.
                self.auto = auto
                self.write("++auto %i" % self.auto)

                # check the version
                self.write("++ver")
                # important not to use the self.readline() method
                # at this stage if auto == 0, otherwise the connector is going
                # to prompt the instrument for a reading an generate an error
                version_number = (self.connection.readline()).decode()

                if "Prologix GPIB-USB Controller" not in version_number:
                    self.connection = None
                    print(
                        "The port %s isn't related to a Prologix controller "
                        "(try to plug and unplug the cable if it is there "
                        "nevertheless)" % (com_port))
                    print(version_number)

                print(
                    "%s is connected on the port '%s'"
                    % (version_number[:-2], com_port)
                )
            else:
                print("The connection to the Prologix connector failed")

    def __str__(self):
        if self.connection is not None:
            self.write("++ver")
            return (self.connection.readline()).decode()
        else:
            return ""

    def controller_id(self):
        return self.__str__()

    def write(self, cmd):
        """use serial.write"""
        # add a new line if the command didn't have one already
        if not cmd.endswith('\n'):
            cmd += '\n'
        if self.connection is not None:
            #  print("Prologix in : ", cmd)
            self.connection.write(cmd.encode())

    def read(self, num_bit):
        """use serial.read"""
        if self.connection is not None:
            if not self.auto:
                self.write('++read eoi')
            answer = self.connection.read(num_bit)
            # print("Prologix out (read) : ", answer)
            return (answer).decode()
        else:
            return ""

    def readline(self):
        """use serial.readline"""
        if self.connection is not None:
            if not self.auto:
                self.write('++read eoi')
            answer = self.connection.readline()
            # print("Prologix out (readline): ", answer)
            return answer.decode()
        else:
            return ""

    def timeout(self, new_timeout=None):
        """
        query the timeout setting of the serial port if no argument provided
        change the
        """
        if self.connection is not None:
            if new_timeout is None:
                return self.connection.timeout
            else:
                old_timeout = self.connection.timeout
                self.connection.timeout = new_timeout
                return old_timeout

    def get_open_gpib_ports(self, num_ports=30):
        """Finds out which GPIB ports are available for prologix controller"""
        open_ports = []

        if not self.mock:

            # sets the timeout to quite fast
            old_timeout = self.timeout(0.1)
            # iterate the ports number
            for i in range(num_ports + 1):

                # change the GPIB address on the prologix controller
                # prove if an instrument is connected to the port
                self.write('++addr %i\n*IDN?\n' % i)

                # probe the answer
                s = self.readline()

                # if it is longer than zero it is an instrument
                # we store the GPIB address
                if len(s) > 0:
                    open_ports.append("GPIB0::%s" % (i))

            # resets the timeout to its original value
            self.timeout(old_timeout)
        #        print "Time out is", self.timeout()

        return open_ports
