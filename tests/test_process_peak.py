#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 10:15:37 2021

@author: haascp
"""
from chromatogram_gen import generate_test_chromatograms, plot_test_data

# MODIFY AS NEEDED FOR TESTING
show_peak_purity_analytics = False

# GET TEST DATA
test_data = generate_test_chromatograms()

# TEST DATA VISUALIZATION
"""
for data in test_data:
    plot_test_data(data)
"""

# ACTUAL TESTS
from mocca.peak.models import PickedPeak


# TEST PEAK PROCESS FUNCTIONS
from mocca.peak.expand import expand_peak

def test_peak_expand_1():
    # expand big peak
    peak = PickedPeak(left=848, right=852, maximum=850, dataset=test_data[0], idx=1)
    expanded_peak = expand_peak(peak, absorbance_threshold=0.1)
    left, right = (expanded_peak.left, expanded_peak.right)
    assert left < 835
    assert right > 865

def test_peak_expand_2():
    # expand small peak
    peak = PickedPeak(left=298, right=302, maximum=300, dataset=test_data[0], idx=1)
    expanded_peak = expand_peak(peak, absorbance_threshold=0.1)
    left, right = (expanded_peak.left, expanded_peak.right)
    assert left < 295
    assert right > 305

def test_peak_expand_3():
    # does not expand correct peak
    peak = PickedPeak(left=285, right=315, maximum=300, dataset=test_data[0], idx=1)
    expanded_peak = expand_peak(peak, absorbance_threshold=0.1)
    left, right = (expanded_peak.left, expanded_peak.right)
    assert abs(left - 285) <= 1  # +/- 1 error for tolerance due to noise
    assert abs(right - 315) <= 1

def test_peak_expand_4():
    # does not expand correct peak, no noise dataset
    peak = PickedPeak(left=285, right=315, maximum=300, dataset=test_data[6], idx=1)
    expanded_peak = expand_peak(peak, absorbance_threshold=0.1)
    left, right = (expanded_peak.left, expanded_peak.right)
    assert left == 285
    assert right == 315


from mocca.peak.check import check_peak_saturation, check_peak_purity, check_peak

def test_peak_saturation_1():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1)
    assert not check_peak_saturation(peak, detector_limit=10)


def test_peak_saturation_2():
    peak = PickedPeak(left=150, right=250, maximum=200, dataset=test_data[2], idx=1)
    assert check_peak_saturation(peak, detector_limit=10)

def test_peak_purity_1():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1) 
    assert check_peak_purity(peak, show_analytics=show_peak_purity_analytics)

def test_peak_purity_2():
    # check for impure peak
    peak = PickedPeak(left=40, right=65, maximum=53, dataset=test_data[3], idx=1)
    assert not check_peak_purity(peak, show_analytics=show_peak_purity_analytics)

def test_peak_purity_3():
    # check for impure peak
    peak = PickedPeak(left=40, right=70, maximum=55, dataset=test_data[4], idx=1)
    assert not check_peak_purity(peak, show_analytics=show_peak_purity_analytics)

def test_peak_purity_4():
    # check for impure and saturated peak
    peak = PickedPeak(left=40, right=70, maximum=55, dataset=test_data[5], idx=1)
    assert not check_peak_purity(peak, show_analytics=show_peak_purity_analytics)

def test_check_peak():
    peak = PickedPeak(left=80, right=120, maximum=100, dataset=test_data[0], idx=1)
    checked_peak = check_peak(peak, detector_limit=10, show_analytics=show_peak_purity_analytics)
    assert not checked_peak.saturation
    assert checked_peak.pure

"""
# TEST PICKED PEAK PROCESSOR
from mocca.peak.processor import process_peak

def test_peak_process():
    peak = PickedPeak(left=298, right=302, maximum=300, dataset=test_data[0], idx=1)
    processed_peak = process_peak(peak, detector_limit=10, absorbance_threshold=0.1)
    assert processed_peak.left < 295
    assert processed_peak.right > 305
    assert processed_peak.saturation is not None
    assert processed_peak.pure is not None
    assert type(processed_peak) == ProcessedPeak
"""
