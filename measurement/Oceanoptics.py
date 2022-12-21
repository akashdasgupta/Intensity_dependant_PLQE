import numpy as np
from seabreeze import spectrometers
import time
import csv

class QEPro(object):
    correct_nonlinearity = True
    correct_dark_counts = True
    _integration_time_ms = 30
    plot_update_ms = 200
    spec = None
    wls = np.array([])
    ax = None
    ani = None
    l = None
    fig = None

    def __init__(self, spec_num=0):
        devs = spectrometers.list_devices()
        try:
            self.spec = spectrometers.Spectrometer(devs[spec_num])
        except Exception as e:
            raise ValueError(f"Can't find a spectrometer, is it powered on and connected? Error message={e}")
        self.integration_time_ms = self._integration_time_ms
        self.update_wls()

    def update_wls(self):
        self.wls = self.spec.wavelengths()
        return self.wls

    @property
    def integration_time_ms(self):
        return self._integration_time_ms

    @property
    def integration_time_limits(self):
        return [i /1000 for i in self.spec.integration_time_micros_limits]
    
    @property
    def max_intensity(self):
        return self.spec.max_intensity

    @integration_time_ms.setter
    def integration_time_ms(self, value):
        try:
            self.spec.integration_time_micros(value * 1000)
            self._integration_time_ms = value
            self.get_counts()  # flush buffer once to forget old data
            self.get_counts()  # flush bufer again to forget in-measurement data
            # For good luck:
            time.sleep(2)
            self.get_counts()

        except Exception as e:
            raise ValueError(f"Error setting integration time to {value}ms: {e}")

    def get_counts(self):
        counts = self.spec.intensities(correct_dark_counts=self.correct_dark_counts, correct_nonlinearity=self.correct_nonlinearity)
        return counts

def noise_calibrate(spec):
    wavel = spec.update_wls()

    for int_time in np.logspace(3,4,3):
        spec.integration_time_ms = int(int_time)
        print(int(int_time))
        holder = []
        for i in range(10):
            counts = spec.get_counts()
            holder.append(counts)
        stds = [np.std([count[i] for count in holder]) for i in range(len(wavel))]

        with open(f"{int(int_time)}.csv",'w',newline='') as file:
            writer = csv.writer(file)
            writer.writerows(zip(wavel, stds))

if __name__ == '__main__':
    spec = QEPro()
    noise_calibrate(spec)

        

