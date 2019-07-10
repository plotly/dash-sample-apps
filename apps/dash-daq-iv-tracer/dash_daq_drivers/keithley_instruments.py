# -*- coding: utf-8 -*-
"""
Driver for the Keithley instruments
Manual for the KT2400 found in 'http://research.physics.illinois.edu/bezryadin/
labprotocol/Keithley2400Manual.pdf'
@author: pierre-francois.duc@netplus.ch
"""
import numpy as np

from .generic_instruments import Instrument, INTF_PROLOGIX


def fake_iv_relation(
    src_type,
    src_val,
    v_oc=20.5,
    i_sc=3.45,
    c1=0.000002694,
    c2=0.077976842
):
    """model of solar cell IV curve
    source: https://www.sciencedirect.com/science/article/pii/S1658365512600120

    src_type should be either 'I' or 'V'
    """
    # Make sure the format is a numpy array
    src_val = np.append(np.array([]), src_val) * 2.1
    # Prepare an answer based on the size of the input
    answer = np.zeros(np.size(src_val))

    if src_type == 'I':
        # Values of the input smaller than the short circuit current
        idx_ok = np.where(src_val < i_sc)
        answer[idx_ok] = c2 * v_oc \
            * np.log(1 + (1 - src_val[idx_ok] / i_sc) / c1)
        return answer
    elif src_type == 'V':
        # Values of the input smaller than the open circuit voltage
        idx_ok = np.where(src_val < v_oc)
        answer[idx_ok] =  \
            i_sc \
            * (1 - c1 * (np.exp(src_val[idx_ok] / (c2 * v_oc)) - 1))
        return answer


INTERFACE = INTF_PROLOGIX

SRC_TYPES = [
    'VOLT',     # Voltage
    'CURR'      # Current
]

SRC_MODES = [
    'FIX',      # Fixed output
    'LIST',     # Variable outputs
    'SWE'       # Sweep outputs
]


class KT2400(Instrument):
    """"driver of the Keithley 2400 SourceMeter"""
    def __init__(self,
                 instr_port_name='',
                 mock_mode=False,
                 instr_user_name='KT 2400',
                 **kwargs):

        # manage the presence of the keyword interface which will determine
        # which method of communication protocol this instrument will use
        if 'interface' in kwargs.keys():

            interface = kwargs.pop('interface')

        else:

            interface = INTERFACE

        instr_mesurands = {
            'V': 'V',   # Voltage in Volt
            'I': 'A'    # Current in Ampere
        }

        if interface == INTF_PROLOGIX:
            kwargs['auto'] = 0

        super(KT2400, self).__init__(instr_port_name,
                                     instr_id_name='KT2400',
                                     instr_user_name=instr_user_name,
                                     mock_mode=mock_mode,
                                     instr_intf=interface,
                                     instr_mesurands=instr_mesurands,
                                     **kwargs)

        self.auto_output_off = False
        self.voltage_compliance = 0
        self.current_compliance = 0

        if instr_port_name:
            self.initialize()

    def _check_arg(self, arg, arg_list):
        """check if the argument is in a list"""
        answer = (arg in arg_list)
        if not answer:
            print("'%s' is not a valid argument" % arg)
            print("Valid arguments are : %s" % str(arg_list))
        return answer

    def _check_is_src_mode(self, arg):
        """check if the argument is a valid source mode"""
        return self._check_arg(arg, arg_list=SRC_MODES)

    def _check_is_src_type(self, arg):
        """check if the argument is a valid source type"""
        return self._check_arg(arg, arg_list=SRC_TYPES)

    def _clear_register(self):
        """refer to p 15-4 of the KT2400 manual"""
        if not self.mock_mode:
            self.write(':STAT:PRES')

    def initialize(self):
        """get the compliance and the auto output parameters"""
        if self.instr_connexion is not None:
            self.auto_output_off = self.enquire_auto_output_off()
            self.voltage_compliance = self.get_voltage_compliance()
            self.current_compliance = self.get_current_compliance()

    def connect(self, instr_port_name, **kwargs):
        super(KT2400, self).connect(instr_port_name, **kwargs)
        self.initialize()

    def measure(self, instr_param):
        if instr_param in self.measure_params:

            if instr_param == 'V':
                if not self.mock_mode:
                    # Initiate a voltage measure (turn output ON)
                    self.write('CONF:VOLT')
                    answer = self.ask(':READ?')
                    # Voltage comes in first position by default
                    answer = float(answer.split(',')[0])
                    # Check that the value is not larger than the compliance
                    if answer >= self.voltage_compliance:
                        print("Measured voltage is at compliance level")
                else:
                    answer = np.random.random()

                self.last_measure[instr_param] = answer

            elif instr_param == 'I':
                if not self.mock_mode:
                    # Initiate a current measure (turn output ON)
                    self.write(':CONF:CURR')
                    answer = self.ask(':READ?')
                    # Voltage comes in second position by default
                    answer = float(answer.split(',')[1])
                    # Check that the value is not larger than the compliance
                    if answer >= self.current_compliance:
                        print("Measured current is above compliance level")
                else:
                    answer = np.random.random()

                self.last_measure[instr_param] = answer

            self.measured_data[instr_param].append(
                self.last_measure[instr_param]
            )

        else:
            print(
                "you are trying to measure a non existent instr_param : "
                + instr_param
            )
            print("existing instr_params :", self.measure_params)
            answer = None
        return answer

    def source_and_measure(self, instr_param, src_val):
        """"set the given source and measure the corresponding measurand
            source voltage => measure current
            source current => measure voltage
        """
        if not self.mock_mode:
            if instr_param == 'V':
                self.configure_voltage_source()
                self.set_voltage(src_val)
                answer = self.measure_current()
                self.disable_output()
            elif instr_param == 'I':
                self.configure_current_source()
                self.set_current(src_val)
                answer = self.measure_voltage()
                self.disable_output()
            else:
                print("The source type should be either 'I' or 'V'")
                answer = np.nan
        else:
            answer = np.round(
                np.squeeze(fake_iv_relation(instr_param, src_val)),
                4
            )
        return answer

    def measure_voltage(self):
        return self.measure('V')

    def get_voltage_compliance(self):
        """voltage compliance in Volt"""
        if not self.mock_mode:
            return float(self.ask(':SENS:VOLT:PROT:LEV?'))
        else:
            return 20

    def measure_current(self):
        return self.measure('I')

    def get_current_compliance(self):
        """current compliance in Ampere"""
        if not self.mock_mode:
            return float(self.ask(':SENS:CURR:PROT:LEV?'))
        else:
            return 1

    def configure_source(self, src_type='CURR', src_mode='FIX'):
        """refer to p 18-73 of the KT2400 manual"""
        if not self.mock_mode:
            if self._check_is_src_type(src_type):
                self.write(':SOUR:FUNC:MODE %s' % src_type)
                if self._check_is_src_mode(src_mode):
                    self.write(':SOUR:%s:MODE %s' % (src_type, src_mode))

    def configure_voltage_source(self, src_mode='FIX'):
        """set the source to output voltage"""
        if not self.mock_mode:
            self.configure_source('VOLT', src_mode)
            self.current_compliance = self.get_current_compliance()

    def configure_current_source(self, src_mode='FIX'):
        """set the source to output current"""
        if not self.mock_mode:
            self.configure_source('CURR', src_mode)
            self.voltage_compliance = self.get_voltage_compliance()

    def set_voltage(self, volt_val):
        """set the voltage for the output, does not turn output on"""
        if not self.mock_mode:
            self.write(':SOUR:VOLT %f' % volt_val)

    def set_current(self, curr_val):
        """set the current (in uA) for the output, does not turn output on"""
        if not self.mock_mode:
            self.write(':SOUR:CURR %f' % (curr_val * 1e-6))

    def enable_output(self):
        """turn the output of the KT2400 on"""
        if not self.mock_mode:
            if not self.auto_output_off:
                print('ENABLING')
                self.write(':OUTP ON;')

    def disable_output(self):
        """shut the output of the KT2400 off"""
        if not self.mock_mode:
            if not self.auto_output_off:
                self.write(':OUTP OFF;')

    def enquire_auto_output_off(self):
        """refer to p. 13 -7 of the KT2400 manual"""
        if not self.mock_mode:
            return bool(int(self.ask(':SOUR:CLE:AUTO?')))
        else:
            return False

    def enable_auto_output_off(self):
        """refer to p. 13 -7 of the KT2400 manual"""
        if not self.mock_mode:
            self.write(':SOUR:CLE:AUTO ON')
            self.auto_output_off = self.ask(':SOUR:CLE:AUTO?')
        else:
            self.auto_output_off = True

    def disable_auto_output_off(self):
        """refer to p. 13 -7 of the KT2400 manual"""
        if not self.mock_mode:
            self.write(':SOUR:CLE:AUTO OFF')
            self.auto_output_off = self.ask(':SOUR:CLE:AUTO?')
        else:
            self.auto_output_off = False


def test_manual_source_and_meas():
    """test source-measure scheme using the low level commands
        first source current and measure voltage,
        then source voltage and measure current
    """
    i = KT2400(
        'GPIB0::11',
        mock_mode=False,
        prologix='COM3',
        auto=0
    )

    i.enable_auto_output_off()

    i.configure_voltage_source()
    i.set_voltage(1.0)
    print(i.measure_current())

    i.configure_current_source()
    i.set_current(0.00001)
    print(i.measure_voltage())


def test_auto_source_and_meas():
    """test source-measure scheme using the ready made method
        first source current and measure voltage,
        then source voltage and measure current
    """
    i = KT2400(
        'GPIB0::11',
        mock_mode=False,
        prologix='COM3',
        auto=0
    )

    i.disable_auto_output_off()

    print(i.ask('*IDN?'))

    print(i.source_and_measure('V', 2))
    print(i.source_and_measure('I', 0.00001))
    print(i.source_and_measure('I', 0.000002))
    print(i.source_and_measure('I', 0.000003))


def test_connect_after_initialization():
    i = KT2400(
        mock_mode=False,
        prologix='COM3',
        auto=0
    )

    i.connect('GPIB0::11')

    i.enable_auto_output_off()

    print(i.source_and_measure('V', 2))
    print(i.source_and_measure('I', 0.00001))
    print(i.source_and_measure('I', 0.000002))
    print(i.source_and_measure('I', 0.000003))


def test_connect_without_prologix():
    i = KT2400(
        mock_mode=False
    )

    i.connect('GPIB0::11')

    i.enable_auto_output_off()

    print(i.ask('*IDN?'))

    print(i.source_and_measure('V', 2))
    print(i.source_and_measure('I', 0.00001))
    print(i.source_and_measure('I', 0.000002))
    print(i.source_and_measure('I', 0.000003))
