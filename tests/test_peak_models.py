#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 11:21:21 2021

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

from mocca.peak.models import BasePeak, PickedPeak, PreprocessedPeak, ProcessedPeak

def test_equals_1():
    peak_1 = BasePeak(left=1, right=3, maximum=2)
    peak_2 = BasePeak(left=1, right=3, maximum=2)
    assert peak_1 == peak_2

def test_equals_2():
    peak_1 = BasePeak(left=1, right=4, maximum=2)
    peak_2 = BasePeak(left=1, right=3, maximum=2)
    assert not peak_1 == peak_2

def test_equals_3():
    peak_1 = PickedPeak(left=1, right=3, maximum=2, dataset=test_data[0], idx=1)
    peak_2 = PickedPeak(left=1, right=3, maximum=2, dataset=test_data[0], idx=2)
    assert peak_1 == peak_2

def test_equals_4():
    peak_1 = PickedPeak(left=1, right=3, maximum=2, dataset=test_data[1], idx=1)
    peak_2 = PickedPeak(left=1, right=3, maximum=2, dataset=test_data[0], idx=1)
    assert not peak_1 == peak_2

def test_equals_5():
    peak_1 = PreprocessedPeak(left=1, right=3, maximum=2, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, matches=[], integral=12)
    peak_2 = PreprocessedPeak(left=1, right=3, maximum=2, dataset=test_data[0], idx=1,
                            saturation=False, pure=True, matches=["compound"], integral=12)
    assert peak_1 == peak_2

def test_equals_6():
    peak_1 = PreprocessedPeak(left=1, right=3, maximum=2, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, matches=None, integral=12)
    peak_2 = PreprocessedPeak(left=1, right=3, maximum=2, dataset=test_data[1], idx=1,
                            saturation=False, pure=False, matches=[], integral=12)
    assert not peak_1 == peak_2

def test_equals_7():
    peak_1 = ProcessedPeak(left=1, right=3, maximum=2, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    peak_2 = ProcessedPeak(left=1, right=3, maximum=2, dataset=test_data[0], idx=1,
                            saturation=False, pure=True, compound_id="compound", integral=12,
                            concentration=12.3)
    assert peak_1 == peak_2

def test_equals_8():
    peak_1 = ProcessedPeak(left=1, right=3, maximum=2, dataset=test_data[0], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    peak_2 = ProcessedPeak(left=1, right=3, maximum=2, dataset=test_data[1], idx=1,
                            saturation=False, pure=False, compound_id=None, integral=12,
                            concentration=None)
    assert not peak_1 == peak_2