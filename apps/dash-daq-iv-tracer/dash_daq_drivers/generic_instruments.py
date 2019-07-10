# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 18:12:03 2018

@author: Pierre-Francois Duc
"""

import serial
import visa

from .communication_utils import PrologixController

# names to manage the different interfaces used to connect to an instrument
INTF_VISA = 'pyvisa'
INTF_PROLOGIX = 'prologix'
INTF_SERIAL = 'serial'
INTF_INTERNAL = 'internal'


class Instrument(object):
    """generic instrument class"""

    def __init__(
        self,
        instr_port_name='',
        instr_id_name='default',
        instr_user_name='default',
        mock_mode=False,
        instr_intf=None,
        instr_mesurands=None,
        **kwargs
    ):

        # should identify the instrument class in a unique way
        self.instr_id_name = instr_id_name

        # provided by the user to have a more meaningful name
        self.instr_user_name = instr_user_name

        # debug mode trigger
        self.mock_mode = mock_mode

        # Instrument measurement attributes

        # contains the names of measurement channels of the instrument
        self.measure_params = []
        # contains the display name of the channels, by default same as above
        self.params_names = {}
        # contains the units of the channels
        self.params_units = instr_mesurands
        # value of the last measure indexed per measurement channel
        self.last_measure = {}
        # array of all measures indexed per measurement channel
        self.measured_data = {}

        # Instrument connexion attributes

        # the name of the port to connect to the instrument
        self.instr_port_name = instr_port_name
        # the name of the interface used to connect to the instrument
        self.instr_intf = instr_intf
        # connexion handle
        self.instr_connexion = None
        # terminaison characters used to communicate with the instrument
        self.term_chars = ""

        for param in instr_mesurands:
            # initializes the first measured value to 0 and the channels'
            # names
            self.measure_params.append(param)
            self.measured_data[param] = []
            self.last_measure[param] = 0
            self.params_names[param] = param

        if self.instr_intf == INTF_VISA:
            # pyvisa version > 1.6
            self.rm = visa.ResourceManager()

        if self.instr_intf == INTF_PROLOGIX and not self.mock_mode:
            # there is only one COM port that the prologix has, then we go
            # through that for all GPIB communications

            if INTF_PROLOGIX in kwargs:

                # the connection is passed as an argument
                if isinstance(kwargs[INTF_PROLOGIX], str):
                    # if it was the COM PORT number we initiate an instance
                    # of prologix controller
                    if "COM" in kwargs[INTF_PROLOGIX]:
                        self.instr_connexion = PrologixController(
                            com_port=kwargs[INTF_PROLOGIX],
                            **kwargs
                        )
                else:
                    # it was the PrologixController instance
                    self.instr_connexion = kwargs[INTF_PROLOGIX]

                    if "Prologix GPIB-USB Controller" \
                            in self.instr_connexion.controller_id():
                        pass
                    else:
                        print(
                            "The controller passed as an argument is not "
                            "the good one"
                        )

            else:
                # the connection doesn't exist so we create it
                print('Searching for Prologix Controller...')
                self.instr_connexion = PrologixController(**kwargs)

        if not self.mock_mode and instr_port_name is not '':
            self.connect(instr_port_name, **kwargs)

    def __str__(self):
        """returns display name for instrument"""

        return "%s" % self.instr_user_name

    def unique_id(self):
        """returns a unique identifier for instruments
            The instrument id combined with the physical port used to
            connect to the instrument should provide a reliable unique
            identifier for instrument connected to the same computer
        """
        return "%s(%s)" % (self.instr_id_name, self.instr_port_name)

    def measure(self, instr_param='', **kwargs):
        """initiate a measure by the instrument
            Should be redefined in children classes
        """
        pass

    def read(self, num_bytes=None):
        """reads data available on the port"""

        if not self.mock_mode:
            if self.instr_intf == INTF_VISA:
                answer = self.instr_connexion.read()
            elif self.instr_intf in (INTF_SERIAL, INTF_PROLOGIX):
                if num_bytes is not None:
                    answer = self.instr_connexion.read(num_bytes)
                else:
                    answer = self.instr_connexion.readline()
            # the provided instrument interface is unknown
            else:
                answer = None
        # in mock mode
        else:
            answer = 'mock_mode_read'

        return answer

    def write(self, msg):
        """writes command to the instrument but does not require a response"""

        if not self.mock_mode:
            if self.instr_intf == INTF_PROLOGIX:
                # make sure the address is the right one
                self.instr_connexion.write(
                    "++addr %s" % self.instr_port_name)
            if self.instr_connexion is not None:
                answer = self.instr_connexion.write(msg + self.term_chars)
            else:
                raise(IOError("There is no physical connexion established \
with the instrument %s" % self.instr_id_name))
        else:
            answer = msg
        return answer

    def ask(self, msg, num_bytes=None):
        """ writes a command to the instrument and reads its reply """

        answer = None

        if not self.mock_mode:
            if self.instr_intf == INTF_VISA:
                answer = self.instr_connexion.ask(msg)
            elif self.instr_intf in (INTF_SERIAL, INTF_PROLOGIX):
                self.write(msg)
                answer = self.read(num_bytes)
        else:
            answer = msg
        return answer

    def connect(self, instr_port_name=None, **kwargs):
        """implements the connexion to the instrument"""

        if instr_port_name is None:
            instr_port_name = self.instr_port_name

        if self.mock_mode:
            print(
                "Connect %s, named %s on port %s, with %s"
                % (
                    self.instr_id_name,
                    self.instr_user_name,
                    self.instr_port_name,
                    self.instr_intf
                )
            )

        else:
            if self.instr_intf == INTF_VISA:
                # make sure the instrument is not already connected
                self.disconnect()
                # connects through the ressource manager (rm)
                self.instr_connexion = self.rm.open_resource(
                    instr_port_name, **kwargs)
            elif self.instr_intf == INTF_SERIAL:
                # make sure the instrument is not already connected
                self.disconnect()

                # kwargs which must be passed as arguments for Serial
                if "term_chars" in kwargs:
                    self.term_chars = kwargs["term_chars"]
                    kwargs.pop("term_chars")

                if "baud_rate" in kwargs:
                    baud_rate = kwargs["baud_rate"]
                    kwargs.pop("baud_rate")

                    self.instr_connexion = serial.Serial(
                        instr_port_name,
                        baud_rate,
                        **kwargs
                    )
                else:
                    self.instr_connexion = serial.Serial(
                        instr_port_name, **kwargs)

            elif self.instr_intf == INTF_PROLOGIX:
                # only keeps the number of the port
                self.instr_port_name = instr_port_name.replace('GPIB0::', '')

                self.instr_connexion.write(
                    ("++addr %s" % self.instr_port_name))

                # the \n termchar is embedded in the PrologixController class
                self.term_chars = ""

            # should raise ErrorIntf
            else:
                pass

        if self.instr_connexion is not None \
                and self.instr_intf != INTF_PROLOGIX:
            self.instr_port_name = instr_port_name

    def disconnect(self):
        """disconnect the instrument"""

        if self.instr_connexion is not None:
            # should write exception handling as we experience it
            # or do that specifically for the children classes
            self.instr_connexion.close()
