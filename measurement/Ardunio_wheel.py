import serial
from NP_843R import *
from Shutter import *
import csv
class Wheel():
    def __init__(self, pm, shutter):
        self.serial_device = serial.Serial("COM7",  baudrate=115200)
        self.fcom = b'\x66'
        self.rcom = b'\x72'
        self.initial_step_range = 30

        self.serial_device.baudrate = 115200

        initial_power = pm.read()[0]

        shutter.on()
        while True:
            for _ in range(self.initial_step_range):
                self.f()
            current_power = pm.read()[0]
            print(current_power)
            if initial_power > current_power and abs(initial_power-current_power)/initial_power > 0.2:
                break
            initial_power = current_power
        for _ in range(self.initial_step_range):
            self.f()
        time.sleep(3)
        self.zero = pm.read()[0]
        shutter.off()

        print('DONE CALIBRATION!!!!', self.zero)

    def f(self):
        self.serial_device.write(b'\x66')
        time.sleep(0.1)
        self.serial_device.flush()
        time.sleep(0.1)

    def r(self):
        self.serial_device.write(b'\x72')
        time.sleep(0.1)
        self.serial_device.flush()
        time.sleep(0.1)

    def reset(self, pm, shutter):
        shutter.on()
        while True:
            for _ in range(self.initial_step_range):
                self.f()
            current_power = pm.read()[0]
            if  abs(current_power-self.zero)/self.zero < 0.1:
                break

        shutter.off()


if __name__ == "__main__":
    shutter = Shutter()
    pm = PowerMeter()
    wheel = Wheel(pm, shutter)
    
    shutter.on()
    powers = []

    with open('cal_570.csv','w', newline='') as file:
        writer = csv.writer(file)
        for _ in range(570):
            wheel.f()
            power = pm.read()[0]
            writer.writerow([power])
        for _ in range(570):
            wheel.r()
            power = pm.read()[0]
            writer.writerow([power])
    shutter.off()

