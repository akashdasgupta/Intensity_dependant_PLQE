import numpy as np
from seabreeze import spectrometers

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

    @integration_time_ms.setter
    def integration_time_ms(self, value):
        try:
            self.spec.integration_time_micros(value * 1000)
            self._integration_time_ms = value
            self.get_counts()  # flush buffer once to forget old data
            self.get_counts()  # flush bufer again to forget in-measurement data
        except Exception as e:
            raise ValueError(f"Error setting integration time to {value}ms: {e}")

    def get_counts(self):
        counts = self.spec.intensities(correct_dark_counts=self.correct_dark_counts, correct_nonlinearity=self.correct_nonlinearity)
        return counts