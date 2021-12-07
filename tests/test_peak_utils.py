import numpy as np
import matplotlib.pyplot as plt
from chromatogram_gen import generate_test_chromatograms, plot_test_data
# TODO: make dataset actually point to an object of DADData class,
# rather than Chromatogram class

# MODIFY AS NEEDED FOR TESTING
print_graphs = False

# GET TEST DATA
test_data = generate_test_chromatograms()

# TEST DATA VISUALIZATION
"""
for data in test_data:
    plot_test_data(data)
"""

# ACTUAL TESTS
from mocca.peak.models import PickedPeak
from mocca.peak.utils import get_peak_data, average_peak_spectrum, is_unimodal

def test_get_peak_data():
    peak_1 = PickedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1)
    peak_data = get_peak_data(peak_1)
    assert isinstance(peak_data, np.ndarray)
    assert peak_data.shape[0] == 342
    assert peak_data.shape[1] == 101

def test_average_peak_spectrum():
    peak_1 = PickedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1)
    peak_spectrum = average_peak_spectrum(peak_1)
    if print_graphs:
        plt.plot(peak_spectrum)
        plt.show()
    assert isinstance(peak_spectrum, list)
    assert len(peak_spectrum) == 342

def test_is_unimodal():
    peak_1 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
    peak_data_1 = get_peak_data(peak_1)
    peak_2 = PickedPeak(left=120, right=180, maximum=150, dataset=test_data[0], idx=1)
    peak_data_2 = get_peak_data(peak_2)
    if print_graphs:
        plt.plot(peak_data_1[0, :])
        plt.plot(peak_data_2[0, :])
        plt.show()
    assert not is_unimodal(peak_data_1[0, :])
    assert is_unimodal(peak_data_2[0, :])
