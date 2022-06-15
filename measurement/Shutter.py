import serial
import codecs
# Stops people from freaking out, but also, plz dont conda update!!:
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import time


class Shutter():
    def __init__(self):
        self.serial_device = serial.Serial("COM4")
        self.opencommand = [255, 1,1]
        self.closecommand = [255, 1,0]
    def on(self):
        self.serial_device.write(bytearray(self.opencommand))
        time.sleep(.5)
    def off(self):
        self.serial_device.write(bytearray(self.closecommand))
        time.sleep(.5)
    def close(self):
        self.off()
        self.serial_device.close()


if __name__ == '__main__':
    shutter = Shutter()

    while True:
        a = input()
        if a == 'on':
            shutter.on()
        else:
            shutter.off()
