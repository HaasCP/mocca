#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 15:34:43 2021

@author: haascp
"""
import logging
import os

from mocca.peak.models import BasePeak

from mocca.dad_data.models import GradientData, CompoundData
from mocca.dad_data.process_funcs import get_peak_locs, merge_peaks, pick_peaks
from mocca.dad_data.utils import sum_absorbance_by_time

file_path = os.path.dirname(os.path.abspath(__file__))
gradient_path = os.path.join(file_path,
                             "example_data/labsolutions_test_gradient.txt")
gradient_data = GradientData("labsolutions", gradient_path)
data_path = os.path.join(file_path,
                         "example_data/labsolutions_api/labsolutions_test_data.txt")
test_input = {'a': 0, 'b': 0.1, 'c': 0, 'd': 1}
compound_data = CompoundData("labsolutions", data_path, test_input, gradient_data)


def test_sum_absorbance_by_time():
    summed_data = sum_absorbance_by_time(compound_data.data)
    assert len(summed_data) == 3001


def test_get_peak_locs():
    absorbance_threshold = 10000
    
    summed_data = sum_absorbance_by_time(compound_data.data)
    new_data_thresh = summed_data.copy()
    new_data_thresh[new_data_thresh < absorbance_threshold] = 0
    
    peak_locs = get_peak_locs(new_data_thresh)
    assert len(peak_locs) == 1

def test_merge_peaks():
    absorbance_threshold = 10000
    
    summed_data = sum_absorbance_by_time(compound_data.data)
    new_data_thresh = summed_data.copy()
    new_data_thresh[new_data_thresh < absorbance_threshold] = 0
    
    peak_locs = get_peak_locs(new_data_thresh)
    # first peak: BasePeak(left=1398, right=1418, maximum=1406)
    peak_locs.append(BasePeak(left=1350, right=1400, maximum=1375))
    
    merged_peaks = merge_peaks(summed_data, peak_locs)
    
    assert len(merged_peaks) == 1
    assert merged_peaks[0].left == 1350
    assert merged_peaks[0].right == 1418
    assert merged_peaks[0].maximum == 1406

def test_pick_peaks():
    absorbance_threshold = 1000
    peaks_high_pass = 1
    peaks_low_pass = 5
    
    chromatogram = pick_peaks(compound_data, absorbance_threshold, peaks_high_pass, peaks_low_pass)
    
    logging.warning("{}".format(len(chromatogram.peaks)))
    assert len(chromatogram.peaks) == 2
    assert any(peak.maximum == 1406 for peak in chromatogram.peaks)
    
    
    