import pytest
import numpy as np
from chromatogram_gen import generate_test_chromatograms, plot_test_data
import matplotlib.pyplot as plt
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


from mocca.peak.utils import check_same_dataset, check_overlap, get_distance_between

def test_check_same_dataset():
    peak_1 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
    peak_2 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[1], idx=1)
    with pytest.raises(Exception):
        check_same_dataset(peak_1, peak_2)

def test_peak_overlap_1():
    # check that non-overlapping peaks are detected as non-overlapping
    peak_1 = PickedPeak(left=100, right=200, maximum=100, dataset=test_data[0], idx=1)
    peak_2 = PickedPeak(left=250, right=350, maximum=300, dataset=test_data[0], idx=1)
    assert not check_overlap(peak_1, peak_2)
    assert not check_overlap(peak_2, peak_1)

def test_peak_overlap_2():
    # check that overlapping peaks are detected as overlapping
    peak_1 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
    peak_2 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
    assert check_overlap(peak_1, peak_2)
    assert check_overlap(peak_2, peak_1)

def test_peak_overlap_3():
    # check that different dataset peaks do throw an exception
    with pytest.raises(Exception):
        peak_1 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
        peak_2 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[1], idx=1)
        check_overlap(peak_1, peak_2)

def test_peak_distance_1():
    # check that peak distances are correct
    peak_1 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
    peak_2 = PickedPeak(left=250, right=350, maximum=300, dataset=test_data[0], idx=1)
    assert get_distance_between(peak_1, peak_2) == 150
    assert get_distance_between(peak_2, peak_1) == 150

def test_peak_distance_2():
    # check that peak distances are correct
    with pytest.raises(Exception):
        peak_1 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
        peak_2 = PickedPeak(left=250, right=350, maximum=300, dataset=test_data[1], idx=1)
        get_distance_between(peak_1, peak_2)

