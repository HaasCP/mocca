#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 13:26:48 2021

@author: haascp
"""
import pytest
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

from mocca.peak.models import PickedPeak
from mocca.chromatogram.utils import check_same_dataset, check_overlap, get_distance_between

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

def test_peak_distance_1():
    # check that peak distances are correct
    peak_1 = PickedPeak(left=100, right=200, maximum=150, dataset=test_data[0], idx=1)
    peak_2 = PickedPeak(left=250, right=350, maximum=300, dataset=test_data[0], idx=1)
    assert get_distance_between(peak_1, peak_2) == 150
    assert get_distance_between(peak_2, peak_1) == 150
