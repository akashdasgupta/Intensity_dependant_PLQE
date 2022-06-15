import time
import serial
from NP_843R import *
from Ardunio_wheel import *
from Oceanoptics import *
from Shutter import *
import csv
import os

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

num_points = 20
division = 30 # in moter steps 
timewait = 1
savepath = './test'

###############

def integration_from_f(f):
    return 100

for i in ['short', 'long', 'dark_s', 'dark_l']:
    for j in ['in', 'out', 'empty']:
        if not os.path.isdir(f"{savepath}/{j}/{i}"):
            os.makedirs(f"{savepath}/{j}/{i}")

'''
index = []
powers = []

for i in range(num_points):    

    for _ in range(division):
        wheel.f()

    shutter.on()
    power = pm.read()[0]
    powers.append(power)
    print(power)
    index.append(i)
    time.sleep(1)

for i in range(num_points):    

    for _ in range(division):
        wheel.r()

    shutter.on()
    power = pm.read()[0]
    powers.append(power)
    index.append(num_points-i)
    time.sleep(1)

shutter.off()


with open('newtest2.csv','w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(zip(index, powers))


'''
for j in ['in', 'out', 'empty']:
    input(f'Place the sample in the {j} position')
    for i in range(num_points):    
        for _ in range(division):
            wheel.f()
        int_time = integration_from_f(i*division)
        shutter.on()
        power = pm.read()[0]
        print(i,power)
        shutter.off()

        shutter.on()
        ### short
        spec.integration_time_ms = int_time
        time.sleep(timewait)
        counts = spec.get_counts()
        save_spec(f"{savepath}/{j}/short", wavel, counts, power, int_time)

        ### long 
        spec.integration_time_ms = int_time*50
        time.sleep(timewait)
        counts = spec.get_counts()
        save_spec(f"{savepath}/{j}/long", wavel, counts, power, int_time*50)

        shutter.off()
        
        ### short
        spec.integration_time_ms = int_time
        time.sleep(timewait)
        counts = spec.get_counts()
        save_spec(f"{savepath}/{j}/dark_s", wavel, counts, power, int_time)

        ### long 
        spec.integration_time_ms = int_time*50
        time.sleep(timewait)
        counts = spec.get_counts()
        save_spec(f"{savepath}/{j}/dark_l", wavel, counts, power, int_time*50)
    
    
    for _ in range(num_points):
        for _ in range(division):
            wheel.r()




shutter.off()



