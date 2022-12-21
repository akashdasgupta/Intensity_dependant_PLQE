import time
import serial
from NP_843R import *
from Ardunio_wheel import *
from Oceanoptics import *
from Shutter import *
import csv
import os
import atexit

shutter = Shutter()
atexit.register(shutter.close)
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
savepath = r'\\cmfs1.physics.ox.ac.uk\cm\akashdasgupta\int_plqe\Flo\21_12_2022_Krad\1A-aged'

PL_ref_position = 1000
###############


PL_ref_position_index = np.argmin([abs(i-PL_ref_position) for i in wavel])

for i in ['short', 'long', 'dark_s', 'dark_l']:
    for j in ['in', 'out', 'empty']:
        if not os.path.isdir(f"{savepath}/{j}/{i}"):
            os.makedirs(f"{savepath}/{j}/{i}")

min_int_time, max_int_time = spec.integration_time_limits

for j in ['in', 'out', 'empty']:
    int_time_short = 50
    int_time_long = 10000

    input(f'Place the sample in the {j} position')
    for i in range(num_points):    
        for _ in range(division):
            wheel.f()
        
        shutter.on()
        time.sleep(1) # gives time for shutter to open before power read

        power = pm.read()[0]
        print(i,power)
        shutter.off()

        shutter.on()
        ### short


        spec.integration_time_ms = int_time_short
        time.sleep(timewait)
        counts = spec.get_counts()

        while np.max(counts) >= spec.max_intensity * 0.8:
            int_time_short = int(0.8 * int_time_short)

            if int_time_short < min_int_time:
                int_time_short = min_int_time
                spec.integration_time_ms = int_time_short
                time.sleep(timewait)
                counts = spec.get_counts()
                break

            spec.integration_time_ms = int_time_short
            time.sleep(timewait)
            counts = spec.get_counts()

        save_spec(f"{savepath}/{j}/short", wavel, counts, power, int_time_short)

        ### long  
        
        spec.integration_time_ms = int_time_long
        time.sleep(timewait)
        counts = spec.get_counts()

        # checking for PL saturation 
        while np.max(counts[PL_ref_position_index:]) >= spec.max_intensity * 0.8:
            int_time_long = int(0.8 * int_time_long)

            if int_time_long < min_int_time:
                int_time_long = min_int_time
                spec.integration_time_ms = int_time_long
                time.sleep(timewait)
                counts = spec.get_counts()
                break

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
    
    
    for _ in range(num_points): # 
        for _ in range(division):
            wheel.r()
    
    pm.read()[0]


# Pushes wheel into transparent position
for _ in range(4):
    for _ in range(division):
        wheel.r()


shutter.off()



