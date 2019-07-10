import random
import numpy

import dash_daq as daq
import dash_html_components as html
import dash_core_components as dcc

try:
    import seabreeze.spectrometers as sb
    from seabreeze.spectrometers import SeaBreezeError
except Exception as e:
    print(e)


# abstract base class to represent spectrometers
class DashOceanOpticsSpectrometer:

    def __init__(self, specLock, commLock):
        self._spec = None                 # spectrometer
        self._specmodel = ''              # model name for graph title
        self._lightSources = {}           # dict of light sources, if any
        self._spectralData = [[], []]     # wavelengths and intensities
        self._controlFunctions = {}       # behaviour upon changing controls
        self._int_time_max = 650000000    # maximum integration time (ms)
        self._int_time_min = 1000         # minimum integration time (ms)
        self.comm_lock = commLock         # for communicating with spectrometer
        self.spec_lock = specLock         # for editing spectrometer values

    # refreshes/populates spectrometer properties
    def assign_spec(self):
        return

    # get data for graph
    def get_spectrum(self):
        return self._spectralData

    # send each command; return successes and failures
    def send_control_values(self, commands):
        return ({}, {})

    # live-update light intensity
    def send_light_intensity(self, val):
        return
    
    # getter methods
    
    def model(self):
        return self._specmodel

    def light_sources(self):
        return self._lightSources

    def int_time_max(self):
        return self._int_time_max

    def int_time_min(self):
        return self._int_time_min
    
    
# non-demo version
class PhysicalSpectrometer(DashOceanOpticsSpectrometer):
    
    def __init__(self, specLock, commLock):
        super().__init__(specLock, commLock)
        try:
            self.spec_lock.acquire()
            self.assign_spec()
        except SeaBreezeError:  # if the spec has already been connected
            pass
        finally:
            self.spec_lock.release()
        self._controlFunctions = {
            'integration-time-input':
            "self._spec.integration_time_micros",

            'nscans-to-average-input':
            "self._spec.scans_to_average",

            'continuous-strobe-toggle-input':
            "self._spec.continuous_strobe_set_enable",

            'continuous-strobe-period-input':
            "self._spec.continuous_strobe_set_period_micros",

            'light-source-input':
            "self.update_light_source"
        }

    def assign_spec(self):
        try:
            self.comm_lock.acquire()
            devices = sb.list_devices()
            self._spec = sb.Spectrometer(devices[0])
            self._specmodel = self._spec.model
            self._lightSources = [{'label': ls.__repr__(), 'value': ls}
                                  for ls in list(self._spec.light_sources())]
            self._int_time_min = self._spec.minimum_integration_time_micros()
        except Exception:
            pass
        finally:
            self.comm_lock.release()

    def get_spectrum(self):
        if self._spec is None:
            try:
                self.spec_lock.acquire()
                self.assign_spec()
            except Exception:
                pass
            finally:
                self.spec_lock.release()
        try:
            self.comm_lock.acquire()
            self._spectralData = self._spec.spectrum(correct_dark_counts=True,
                                                     correct_nonlinearity=True)
        except Exception:
            pass
        finally:
            self.comm_lock.release()

        return self._spectralData

    def send_control_values(self, commands):
        failed = {}
        succeeded = {}
        
        for ctrl_id in commands:
            try:
                self.comm_lock.acquire()
                eval(self._controlFunctions[ctrl_id])(commands[ctrl_id])
                succeeded[ctrl_id] = str(commands[ctrl_id])
            except Exception as e:
                failed[ctrl_id] = str(e).strip('b')
            finally:
                self.comm_lock.release()
                
        return(failed, succeeded)

    def send_light_intensity(self, lightSource, intensity):
        try:
            self.comm_lock.acquire()
            lightSource.set_intensity(intensity)
        except Exception:
            pass
        finally:
            self.comm_lock.release()
            
    def model(self):
        try:
            self.spec_lock.acquire()
            self.assign_spec()
        except SeaBreezeError:
            pass
        finally:
            self.spec_lock.release()
        return self._specmodel

    def light_sources(self):
        try:
            self.spec_lock.acquire()
            self.assign_spec()
        except SeaBreezeError:
            pass
        finally:
            self.spec_lock.release()
            
        return self._lightSources

    def int_time_max(self):
        try:
            self.spec_lock.acquire()
            self.assign_spec()
        except SeaBreezeError:
            pass
        finally:
            self.spec_lock.release()
            
        return self._int_time_max

    def int_time_min(self):
        try:
            self.spec_lock.acquire()
            self.assign_spec()
        except SeaBreezeError:
            pass
        finally:
            self.spec_lock.release()
            
        return self._int_time_min

    def update_light_source(self, ls):
        if(ls is not None and ls is not ""):
            ls.set_enable(True)

        
class DemoSpectrometer(DashOceanOpticsSpectrometer):

    def __init__(self, specLock, commLock):
        super().__init__(specLock, commLock)
        try:
            self.spec_lock.acquire()
            self.assign_spec()
        except Exception:
            pass
        finally:
            self.spec_lock.release()
        self.controlFunctions = {
            'integration-time-input':
            "self.integration_time_demo",

            'nscans-to-average-input':
            "self.empty_control_demo",

            'continuous-strobe-toggle-input':
            "self.empty_control_demo",

            'continuous-strobe-period-input':
            "self.empty_control_demo",

            'light-source-input':
            "self.exception_demo"
        }
        self._sample_data_scale = self._int_time_min
        self._sample_data_add = 0

    def assign_spec(self):
        self._specmodel = "USB2000+"
        self._lightSources = [{'label': 'Lamp 1 at 127.0.0.1', 'value': 'l1'},
                              {'label': 'Lamp 2 at 127.0.0.1', 'value': 'l2'}]

    def get_spectrum(self, int_time_demo_val=1000):
        self._spectralData[0] = numpy.linspace(400, 900, 5000)
        self._spectralData[1] = [self.sample_spectrum(wl)
                                 for wl in self._spectralData[0]]

        return self._spectralData

    def send_control_values(self, commands):
        failed = {}
        succeeded = {}

        for ctrl_id in commands:
            try:
                eval(self.controlFunctions[ctrl_id])(commands[ctrl_id])
                succeeded[ctrl_id] = str(commands[ctrl_id])
            except Exception as e:
                failed[ctrl_id] = str(e)

        return(failed, succeeded)

    def send_light_intensity(self, lightSource, intensity):
        if(lightSource == 'l1'):
            return
        elif(lightSource == 'l2'):
            self._sample_data_add = intensity
        else:
            self._sample_data_add = 0
            
    def model(self):
        return self._specmodel
    
    def light_sources(self):
        return self._lightSources

    def int_time_max(self):
        return self._int_time_max

    def int_time_min(self):
        return self._int_time_min

    # demo-specific methods
    
    # generates a sample spectrum that's normally distributed about 500 nm
    def sample_spectrum(self, x):
        return (self._sample_data_scale * (numpy.e**(-1 * ((x-500) / 5)**2) +
                                           0.01 * random.random()) +
                self._sample_data_add * 10)

    def integration_time_demo(self, x):
        self._sample_data_scale = x

    def empty_control_demo(self, _):
        return

    def exception_demo(self, x):
        if(x == "l1"):
            raise Exception("Lamp not found.")
        else:
            return


# class to represent all controls
class Control:
    def __init__(self, new_ctrl_id, new_ctrl_name,
                 new_component_type, new_component_attr):
        self.ctrl_id = new_ctrl_id                # id for callbacks
        self.ctrl_name = new_ctrl_name            # name for label
        self.component_type = new_component_type  # dash-daq component type
        self.component_attr = new_component_attr  # component attributes

    # creates a new control box with defined component, id, and name
    def create_ctrl_div(self, pwrOff):
        # create dash-daq components
        try:
            component_obj = getattr(daq, self.component_type)
        except AttributeError:
            component_obj = getattr(dcc, self.component_type)

        # disable if power is off
        self.component_attr['disabled'] = pwrOff
            
        component = component_obj(**self.component_attr)

        # generate html code
        new_control = html.Div(
            id=self.ctrl_id,
            children=[
                html.Div(
                    className='option-name',
                    children=[
                        self.ctrl_name
                    ]
                ),
                component
            ]
        )
        return new_control

    # gets whether we look for "value", "on", etc.
    def val_string(self):
        if('value' in self.component_attr):
            return 'value'
        elif('on' in self.component_attr):
            return 'on'

    # changes value ('on' or 'value', etc.)
    def update_value(self, new_value):
        self.component_attr[self.val_string()] = new_value

