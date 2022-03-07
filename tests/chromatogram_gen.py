import numpy as np
import pandas as pd
import random
import math
import matplotlib.pyplot as plt
from sim_test_data import test_spectra, wavelength


def gaussian(mean, variance, left, right, num_points):
    # returns a Gaussian distribution over the specified range from
    # left to right, with the specified mean and variance
    # the output is a length num_points array
    return np.array([1/math.sqrt(2 * math.pi * variance) *
                    math.exp(-1 * ((x - mean)/variance)**2 / 2)
                    for x in np.linspace(left, right, num_points)])


def generate_elution_profile(elution_time, elution_width, length):
    # generates gaussian elution profile based on elution time,
    # elution width, and spectral length
    return gaussian(elution_time, elution_width, 0, length, length)


def generate_chromatogram(elution_times, elution_widths, concs, length):
    # generates a len(wavelength) x length array of simulated data
    # concs, elution_times, and elution_widths should be
    # 10-element lists corresponding to concentrations, elution times,
    # and elution widths of each of the 10 test compounds
    output_spec = np.zeros((len(wavelength), length))
    for spectra, conc, time, width in \
            zip(test_spectra, concs, elution_times, elution_widths):
        output_spec += conc * np.array(spectra).reshape((-1, 1)) * \
                       generate_elution_profile(time, width, length).reshape((1, -1))
    return output_spec


def generate_random_chromatogram_param(length):
    # generates random parameters for a chromatogram
    times = []
    widths = []
    concs = []
    for i in range(10):
        times.append(math.floor(random.random() * length * 3/4 + length * 1/8))
        widths.append(math.floor(random.random() * length * 1/200 + length * 1/500))
        concs.append(random.random() * 10)
    return times, widths, concs

def generate_test_chromatograms():
    # TEST DATA GENERATION
    test_data = []
    test_data.append(TestChromatogram(1000, concs=[1] * 10,
                                      elution_widths=[5, 10, 5, 5, 3, 10, 3, 6, 15, 15],
                                      elution_times=[100, 150, 200, 250, 300,
                                                     400, 500, 700, 800, 850]))
    test_data.append(TestChromatogram(1000))  # random chromatogram
    test_data.append(TestChromatogram(1000, concs=[100000] + [0] * 9,
                                      elution_times=[200] + [0] * 9))
    test_data.append(TestChromatogram(1000, concs=[1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                                      elution_widths=[10]*10,
                                      elution_times=[50, 0, 0, 0, 55, 0, 0, 0, 0, 0]))  # overlap
    test_data.append(TestChromatogram(1000, concs=[1, 0, 0, 0, 0.2, 0, 0.2, 0, 1, 0],
                                      elution_widths=[15]*10,
                                      elution_times=[50, 0, 0, 0, 51, 0, 53, 0, 60, 0]))  # overlap
    test_data.append(TestChromatogram(1000, concs=[10000, 0, 0, 0, 2000, 0, 2000, 0, 10000, 0],
                                      elution_widths=[15]*10,
                                      elution_times=[50, 0, 0, 0, 51, 0, 53, 0, 60, 0]))
    test_data.append(TestChromatogram(1000, concs=[1] * 10,
                                      elution_widths=[5, 10, 5, 5, 3, 10, 3, 6, 15, 15],
                                      elution_times=[100, 150, 200, 250, 300,
                                              400, 500, 700, 800, 850], add_noise=False))
    test_data.append(TestChromatogram(1000, concs=[1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      elution_widths=[10]*10, elution_times=[150]*10))
    return test_data

def plot_test_data(data):
    """
    data is an input of class Chromatogram

    For use in showing test data plots only, not for actual testing.
    """
    plt.plot(np.sum(data.data, axis=0))
    plt.show()
    plt.plot(data.data)
    plt.show()

class TestChromatogram():
    def __init__(self, length, concs=None, elution_times=None,
                 elution_widths=None, add_noise=True):
        self.times, self.widths, self.concs = generate_random_chromatogram_param(length)
        if elution_times:
            self.times = elution_times
        if elution_widths:
            self.widths = elution_widths
        if concs:
            self.concs = concs
        self.data = generate_chromatogram(self.times, self.widths, self.concs, length)
        if add_noise:
            self.data += np.random.normal(0, 0.000003, (len(wavelength), length))
        self.time = list(range(length))
        self.wavelength = wavelength.copy()

    #def contour_plot(self):
    #    plot = contour_map(data=self.data, time=self.time, wavelength=self.wavelength)
    #    return plot

    #def summed_plot(self):
    #    timepoints = self.time
    #    absorbances = np.sum(self.data, axis=0)
    #    df = pd.DataFrame()
    #    df['Timepoint'] = timepoints
    #    df['Absorbance'] = absorbances
    #    plot = plot_1D_data(df=df, xlabel='Timepoint', ylabel='Absorbance',
    #                        title='Summed Data Plot')
    #    return plot

    def write_csv(self, csv_path, time_interval=0.001):
        df = pd.DataFrame(self.data.T, columns=np.array(wavelength),
                          index=np.array(self.time) * time_interval)
        df.to_csv(csv_path)
