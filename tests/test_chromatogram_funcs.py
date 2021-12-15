#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 17:55:32 2021

@author: haascp
"""
import logging
import copy

from mocca.peak.database import PeakDatabase
from mocca.peak.models import IntegratedPeak, ProcessedPeak
from mocca.components.databases import QualiComponentDatabase
from mocca.dad_data.process_funcs import pick_peaks

from mocca.peak.expand import expand_peak
from mocca.peak.check import check_peak
from mocca.peak.integrate import integrate_peak

from chromatogram_gen import generate_test_chromatograms, plot_test_data

# GET TEST DATA
test_data = generate_test_chromatograms()

# TEST DATA VISUALIZATION
"""
for data in test_data:
    plot_test_data(data)
"""

from mocca.chromatogram.funcs import get_pure_istd_peak, get_istd_peak, get_istd_offset, correct_istd_offset

def create_test_peak_db(test_data, absorbance_threshold=0.1, detector_limit=10,
                        print_purity_check=False):
    chromatogram = pick_peaks(test_data, absorbance_threshold, None, None)
    integrated_peaks = []
    for picked_peak in chromatogram.peaks:
        expanded_peak = expand_peak(picked_peak, absorbance_threshold)
        checked_peak = check_peak(expanded_peak, detector_limit,
                                  show_analytics=print_purity_check)
        integrated_peak = integrate_peak(checked_peak)
        integrated_peaks.append(integrated_peak)
    
    processed_peaks = []
    for i, peak in enumerate(integrated_peaks):
        new_peak = ProcessedPeak(left=peak.left,
                                 right=peak.right,
                                 maximum=peak.maximum,
                                 dataset=peak.dataset,
                                 idx=peak.idx,
                                 saturation=peak.saturation,
                                 pure=peak.pure,
                                 integral=peak.integral,
                                 offset=0,
                                 compound_id=str(i+1),
                                 concentration=peak.integral/(i+1))
        processed_peaks.append(new_peak)

    test_db = PeakDatabase()
    for peak in processed_peaks:
        test_db.insert_peak(peak)
    return test_db

def create_shifted_test_chrom(test_data, absorbance_threshold=0.1, detector_limit=10,
                              print_purity_check=False):
    # Shortens time by 10 --> left shift
    new_test_data = copy.deepcopy(test_data)
    new_test_data.data = new_test_data.data[:, 12:] # will be found even with rel_threshold 0.1 since we double here
    new_test_data.time = new_test_data.time[12:]
    chromatogram = pick_peaks(new_test_data, absorbance_threshold, None, None)
    integrated_peaks = []
    for picked_peak in chromatogram.peaks:
        expanded_peak = expand_peak(picked_peak, absorbance_threshold)
        checked_peak = check_peak(expanded_peak, detector_limit,
                                  show_analytics=print_purity_check)
        integrated_peak = integrate_peak(checked_peak)
        integrated_peaks.append(integrated_peak)
    
    chromatogram.peaks = integrated_peaks
    return chromatogram

def test_get_pure_istd_peak():
    peak_db = create_test_peak_db(test_data[0])
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    
    chromatogram = create_shifted_test_chrom(test_data[0])
    
    istd_peak = get_pure_istd_peak(chromatogram, '3', quali_db, 
                           spectrum_correl_coef_thresh=0.99, relative_distance_thresh=0.01)

def test_get_istd_peak_1():
    peak_db = create_test_peak_db(test_data[0])
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    
    chromatogram = create_shifted_test_chrom(test_data[0])
    
    istd_peak = get_istd_peak(chromatogram, '11', quali_db, 
                           spectrum_correl_coef_thresh=0.99, relative_distance_thresh=0.01)
    assert istd_peak is None  # since istd_key not in quali_db

def test_get_istd_peak_2():
    peak_db = create_test_peak_db(test_data[0])
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    
    chromatogram = create_shifted_test_chrom(test_data[0])
    
    istd_peak = get_istd_peak(chromatogram, '3', quali_db, 
                           spectrum_correl_coef_thresh=0.99, relative_distance_thresh=0.01)
    assert istd_peak.maximum == 188

def test_get_istd_offset_1():
    peak_db = create_test_peak_db(test_data[0])
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    
    chromatogram = create_shifted_test_chrom(test_data[0])
    
    offset = get_istd_offset(chromatogram, '3', quali_db, 
                             spectrum_correl_coef_thresh=0.99, relative_distance_thresh=0.01)
    assert offset == -12

def test_get_istd_offset_2():
    peak_db = create_test_peak_db(test_data[0])
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    
    chromatogram = create_shifted_test_chrom(test_data[0])
    
    offset = get_istd_offset(chromatogram, '11', quali_db, 
                             spectrum_correl_coef_thresh=0.99, relative_distance_thresh=0.01)
    assert offset == 0

def test_correct_istd_offset():
    peak_db = create_test_peak_db(test_data[0])
    quali_db = QualiComponentDatabase()
    quali_db.update(peak_db, peak_filter_function=None)
    
    chromatogram = create_shifted_test_chrom(test_data[0])
    logging.warning("{}".format(chromatogram.peaks[2]))
    
    chromatogram = correct_istd_offset(chromatogram, '3', quali_db, 
                             spectrum_correl_coef_thresh=0.99, relative_distance_thresh=0.01)
    logging.warning("{}".format(chromatogram.peaks[2]))
    assert chromatogram.peaks[2].maximum == quali_db['3'].maximum

#TODO add impure peak tests
    


    
