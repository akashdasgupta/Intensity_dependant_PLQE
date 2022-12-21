import win32com.client
import time
import traceback

class PowerMeter():
    def __init__(self):
        self.COM = win32com.client.Dispatch("OphirLMMeasurement.CoLMMeasurement")
        self.COM.StopAllStreams() 
        self.COM.CloseAll()

        DeviceList = self.COM.ScanUSB()
        for Device in DeviceList:   	# if any device is connected
            DeviceHandle = self.COM.OpenUSBDevice(Device)	# open first device
            exists = self.COM.IsSensorExists(DeviceHandle, 0)
            if exists:
                self.DeviceHandle = DeviceHandle
                break
            else:
                print(f'No Sensor attached to {Device} !!!')
                self.DeviceHandle = None
                break
        
    
    def read(self):
        self.COM.StartStream(self.DeviceHandle, 0)
        time.sleep(.2)				# wait a little for data
        data = self.COM.GetData(self.DeviceHandle, 0)
        data = self.COM.GetData(self.DeviceHandle, 0)
        data = self.COM.GetData(self.DeviceHandle, 0)
        while len(data[0]) == 0:
            data = self.COM.GetData(self.DeviceHandle, 0)
            time.sleep(.2)
            data = self.COM.GetData(self.DeviceHandle, 0)
        reading = float(data[0][0])
        timestap = float(data[1][0])
        status = data[2][0]

        self.COM.StopAllStreams()
        # self.COM.CloseAll()
        return reading, timestap, status

    def close():
        # Stop & Close all devices
        self.COM.StopAllStreams()
        self.COM.CloseAll()
        # Release the object
        self.COM = None




