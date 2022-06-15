import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pyvisa as visa
from keithley2600 import Keithley2600

# We will check if resource manager is open, opens if not, and extracts Keithly resource string:
try:
    rm.list_resources()
except NameError:
    rm = visa.ResourceManager(r'C:\Windows\System32\visa64.dll')
    resources = rm.list_resources()
for resource in resources:
    if resource.startswith("GPIB0"): 
        keithly_string = resource


print(keithly_string)
class Keithley():
    def __init__(self):
            print("Attempting to connect to the  Keithley sourcemeter...")
            try:
                # Opens device:
                self.sm = Keithley2600(keithly_string, raise_keithley_errors=True, visa_library=r'C:\Windows\System32\visa64.dll')
                # Reset both channels for consistancy:
                self.sm.smua.reset()
                self.sm.smub.reset()
                print("Done!")
            except Exception as ex:
                self.sm = None
                print('Could not initilise Keithley sourcemeter! ID given:', keithly_string)
                print(ex)
            self.voltages = []
            self.currents = []

    def measure(self, channel='b'):
        """
        Single measurement: just returns voltage and current measured at a channel
        """
        if channel.lower() == 'a':
            # self.sm.smua.reset()
            # self.sm.smua.sense = self.sm.smua.SENSE_REMOTE
            v = self.sm.smua.measure.v()
            i = self.sm.smua.measure.i()
            return v,i
        elif channel.lower() == 'b':
            # self.sm.smub.reset()
            # self.sm.smub.sense = self.sm.smub.SENSE_REMOTE
            v = self.sm.smub.measure.v()
            i = self.sm.smub.measure.i()
            return v,i
        else: # Default is b 
            print("The provided channel:", channel,"was invalid. Defaulting to b...")
            # self.sm.smub.reset()
            # self.sm.smub.sense = self.sm.smub.SENSE_REMOTE
            v = self.sm.smub.measure.v()
            i = self.sm.smub.measure.i()
            return v,i
        
    def set_oc(self, channel='b'):
        if channel.lower() == 'a':
            self.sm.smua.reset()
        elif channel.lower() == 'b':
            self.sm.smub.reset()
        
        self.set_current_level(0,channel)
 
    def set_current_level(self, level, channel='b'):
        if channel.lower() == 'a':
            self.sm.smua.reset()
            self.sm.smua.source.output = self.sm.smua.OUTPUT_ON
            self.sm.smua.source.func = self.sm.smua.OUTPUT_DCAMPS
            self.sm.smua.source.leveli = level 
        elif channel.lower() == 'b':
            self.sm.smub.reset()

            self.sm.smub.source.output = self.sm.smub.OUTPUT_ON
            self.sm.smub.source.func = self.sm.smub.OUTPUT_DCAMPS
            self.sm.smub.source.leveli = level 
        else:
            print("The provided channel:", channel,"was invalid. Defaulting to b...")
            self.sm.smub.reset()
            self.sm.smub.source.output = self.sm.smub.OUTPUT_ON
            self.sm.smub.source.func = self.sm.smub.OUTPUT_DCAMPS
            self.sm.smub.source.leveli = level    
    
    def set_voltage_level(self, level, channel='b'):  
        if channel.lower() == 'a':
            self.sm.smua.source.output = self.sm.smua.OUTPUT_ON
            self.sm.smua.source.func = self.sm.smua.OUTPUT_DCVOLTS
            self.sm.smua.source.levelv = level 
        elif channel.lower() == 'b':
            self.sm.smub.source.output = self.sm.smub.OUTPUT_ON
            self.sm.smub.source.func = self.sm.smub.OUTPUT_DCVOLTS
            self.sm.smub.source.levelv = level 
        else:
            print("The provided channel:", channel,"was invalid. Defaulting to b...")
            self.sm.smub.source.output = self.sm.OUTPUT_ON
            self.sm.smub.source.func = self.sm.OUTPUT_DCVOLTS
            self.sm.smub.source.levelv = level   
    
    def off(self, channel='b'):
        if channel.lower() == 'a':
            self.sm.smua.reset()
            self.sm.smua.source.output = self.sm.smua.OUTPUT_OFF
        elif channel.lower() == 'b':
            self.sm.smub.reset()
            self.sm.smub.source.output = self.sm.smub.OUTPUT_OFF
        else:
            print("The provided channel:", channel,"was invalid. Defaulting to b...")
            self.sm.smub.reset()
            self.sm.smub.source.output = self.sm.smub.OUTPUT_OFF






