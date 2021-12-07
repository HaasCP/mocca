#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 11:14:01 2021

@author: haascp
"""

import logging
import matplotlib.pyplot as plt

from mocca.peak.models import PickedPeak, ProcessedPeak
from mocca.peak.utils import average_peak_spectrum
from mocca.databases.peak_db import PeakDatabase

from chromatogram_gen import generate_test_chromatograms, plot_test_data

from mocca.databases.utils import (get_valid_peaks, filter_peaks, 
                                   get_filtered_peaks, sort_peaks_by_compound,
                                   get_filtered_peaks_by_compound,
                                   average_spectra_over_peaks, 
                                   average_ret_times_over_peaks)

# GET TEST DATA
test_data = generate_test_chromatograms()

# TEST DATA VISUALIZATION
"""
for data in test_data:
    plot_test_data(data)
"""

# MODIFY AS NEEDED FOR TESTING
print_graphs = False

# TEST PEAKS

peaks = [] 
peaks.append(ProcessedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1, 
                          saturation=True, pure=True, compound_id="1", integral=123, concentration=12.3))
peaks.append(ProcessedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1, 
                          saturation=False, pure=True, compound_id="2", integral=123, concentration=12.3))
peaks.append(ProcessedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1, 
                          saturation=False, pure=False, compound_id="3", integral=123, concentration=12.3))
peaks.append(ProcessedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1, 
                          saturation=False, pure=True, compound_id=None, integral=123, concentration=None))
peaks.append(ProcessedPeak(left=180, right=220, maximum=200, dataset=test_data[0], idx=1, 
                          saturation=True, pure=True, compound_id="1", integral=123, concentration=12.3))
peaks.append(ProcessedPeak(left=180, right=220, maximum=200, dataset=test_data[0], idx=1, 
                          saturation=False, pure=True, compound_id="2", integral=123, concentration=12.3))
peaks.append(ProcessedPeak(left=180, right=220, maximum=200, dataset=test_data[0], idx=1, 
                          saturation=False, pure=False, compound_id="3", integral=123, concentration=12.3))
peaks.append(ProcessedPeak(left=180, right=220, maximum=200, dataset=test_data[0], idx=1, 
                          saturation=False, pure=True, compound_id=None, integral=123, concentration=None))

# TESTS

def test_get_valid_peaks():
    valid_peaks = get_valid_peaks(peaks)
    assert len(valid_peaks) == 2
    assert valid_peaks[0] == peaks[1]

def filter_func_1(peaks):
    return [peak for peak in peaks if peak.left > 150]

def filter_func_2(peaks):
    return [peak for peak in peaks if peak.compound_id is not None and peak.pure is True]

def test_filter_peaks_1():
    filtered_peaks = filter_peaks(peaks, filter_func_1)
    assert len(filtered_peaks) == 4

def test_filter_peaks_2():
    filtered_peaks = filter_peaks(peaks, filter_func_2)
    assert len(filtered_peaks) == 4

def test_get_filtered_peaks():
    test_db = PeakDatabase()
    test_db.peaks = peaks
    result = get_filtered_peaks(test_db, filter_func_1)
    assert len(result) == 1
    assert result[0] == peaks[5]

def test_sort_peaks_by_compound():
    peak_dict = sort_peaks_by_compound(peaks)
    assert len(peak_dict["1"]) == 2 and len(peak_dict["2"]) == 2 and len(peak_dict["3"]) == 2
    assert peak_dict["1"][0] == peaks[0]

def test_get_filtered_peaks_by_compound():
    test_db = PeakDatabase()
    test_db.peaks = peaks
    peak_dict = get_filtered_peaks_by_compound(test_db, filter_func_1)
    assert len(peak_dict["2"]) == 1

def test_average_spectra_over_peaks():
    peak_1 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
    peak_2 = PickedPeak(left=110, right=190, maximum=150, dataset=test_data[0], idx=1)
    peak_3 = PickedPeak(left=120, right=180, maximum=150, dataset=test_data[0], idx=1)
    peak_spectrum_1 = average_peak_spectrum(peak_1)
    peak_spectrum_2 = average_peak_spectrum(peak_2)
    peak_spectrum_3 = average_peak_spectrum(peak_3)
    spectrum = average_spectra_over_peaks([peak_1, peak_2, peak_3])
    if print_graphs:
        plt.plot(peak_spectrum_1, linestyle='dashed')
        plt.plot(peak_spectrum_2, linestyle='dashed')
        plt.plot(peak_spectrum_3, linestyle='dashed')
        plt.plot(spectrum)
        plt.show()
    assert isinstance(spectrum, list)
    assert len(spectrum) == 342

def test_average_ret_times_over_peaks():
    peak_1 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
    peak_2 = PickedPeak(left=110, right=190, maximum=140, dataset=test_data[0], idx=1)
    peak_3 = PickedPeak(left=120, right=180, maximum=160, dataset=test_data[0], idx=1)
    left, right, maximum = average_ret_times_over_peaks([peak_1, peak_2, peak_3])
    assert left == 110 and right == 190 and maximum == 150

#logging.warning("{}".format(peak_dict))