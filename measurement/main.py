import time
import serial
from NP_843R import *
from Ardunio_wheel import *
from Oceanoptics import *
from Shutter import *
import csv
import os

detector_effective_bits = 17.65 # max here

shutter = Shutter()
pm = PowerMeter()
wheel = Wheel(pm, shutter)

def save_spec(savepath, wavel, counts, power='UNKNOWN',integration='UNKNOWN'):
    with open(f"{savepath}/{power}_{integration}.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(zip(wavel,counts))


spec = QEPro()
wavel = spec.update_wls()

##############

num_points = 19
division = 30 # in moter steps 
timewait = 1
savepath = './test4_2'

###############

def integration_from_f(f):
    return 50

for i in ['short', 'long', 'dark_s', 'dark_l']:
    for j in ['in', 'out', 'empty']:
        if not os.path.isdir(f"{savepath}/{j}/{i}"):
            os.makedirs(f"{savepath}/{j}/{i}")



for j in ['in', 'out', 'empty']:
    int_time_short = 50
    int_time_long = 5000

    input(f'Place the sample in the {j} position')
    for i in range(num_points):    
        for _ in range(division):
            wheel.f()
        
        shutter.on()
        power = pm.read()[0]
        print(i,power)
        shutter.off()

        shutter.on()
        ### short


        spec.integration_time_ms = int_time_short
        time.sleep(timewait)
        counts = spec.get_counts()

        while np.max(counts) >= (2**detector_effective_bits) * 0.8:
            int_time_short = int(0.8 * int_time_short)
            spec.integration_time_ms = int_time_short
            time.sleep(timewait)
            counts = spec.get_counts()

        save_spec(f"{savepath}/{j}/short", wavel, counts, power, int_time_short)

        ### long 
        spec.integration_time_ms = int_time_long
        time.sleep(timewait)
        counts = spec.get_counts()
        save_spec(f"{savepath}/{j}/long", wavel, counts, power, int_time_long)

        shutter.off()
        
        ### short
        spec.integration_time_ms = int_time_short
        time.sleep(timewait)
        counts = spec.get_counts()
        save_spec(f"{savepath}/{j}/dark_s", wavel, counts, power, int_time_short)

        ### long 
        spec.integration_time_ms = int_time_long
        time.sleep(timewait)
        counts = spec.get_counts()
        save_spec(f"{savepath}/{j}/dark_l", wavel, counts, power, int_time_long)
    
    
    for _ in range(num_points):
        for _ in range(division):
            wheel.r()
    
    pm.read()[0]




shutter.off()



