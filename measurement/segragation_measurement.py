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
wheel = Wheel(pm, shutter, rezero=False)

spec = QEPro()
wavel = spec.update_wls()

num_points = 60
timestep = 1
integration_time = 2
savepath = r"\\cmfs1.physics.ox.ac.uk\cm\seos\PLQE\segregation\20240321 test\FACsCl\1"

spec.integration_time_ms(integration_time)

if not os.path.isdir(savepath):
    os.makedirs(savepath)


starttime = time.time()
shutter.on()
power = pm.read()[0]

all_counts = []
times = []
for _ in range(num_points):
    counts = spec.get_counts()
    curr_time = time.time() - starttime
    with open(f"{savepath}/{curr_time}_{power}_{integration_time}.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(zip(wavel,counts))
    time.sleep(timestep)


shutter.off()   
time.sleep(1)
input("remove sample")  

shutter.on()   
time.sleep(5)
counts = spec.get_counts()
curr_time = time.time() - starttime
with open(f"{savepath}/EMPTY.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(zip(wavel,counts))
time.sleep(timestep)


shutter.off()   
time.sleep(1)
counts = spec.get_counts()
curr_time = time.time() - starttime
with open(f"{savepath}/DARK.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(zip(wavel,counts))
time.sleep(timestep)




