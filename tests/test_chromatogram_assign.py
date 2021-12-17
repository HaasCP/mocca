#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 14:03:26 2021

@author: haascp
"""
import logging

from mocca.peak.models import PreprocessedPeak, ProcessedPeak
from mocca.peak.database import PeakDatabase
from mocca.chromatogram.model import Chromatogram

from mocca.chromatogram.assign import (sort_peaks_by_best_match, 
                                  get_best_match_compound_id, update_peaks_and_matches,
                                  assign_best_match_peak, assign_matched_peaks,
                                  get_next_unknown_id, assign_unmatched_peaks_react,
                                  assign_peaks_react)

def create_similarity_dicts(compound_id_num, n_dicts):
    simil_dicts = []
    for i in range(n_dicts):
        dic = {}
        dic['compound_id'] = str(compound_id_num+i)
        dic['spectrum_correl_coef'] = 1 - 0.01 * (i + (compound_id_num*0.5))
        dic['distance'] = 0
        dic['relative_distance'] = 0
        simil_dicts.append(dic)
    return simil_dicts

def create_test_peaks():
    simil_dicts_1 = create_similarity_dicts(1, 3)
    simil_dicts_2 = create_similarity_dicts(2, 2)
    simil_dicts_3 = create_similarity_dicts(3, 1)
    
    peak_1 = PreprocessedPeak(left=90, right=110, maximum=100, dataset="a", idx=1,
                              pure=False, saturation=False, integral=0.99, offset=0, matches=simil_dicts_3)
    peak_2 = PreprocessedPeak(left=190, right=210, maximum=200, dataset="b", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=simil_dicts_1)
    peak_3 = PreprocessedPeak(left=290, right=310, maximum=300, dataset="c", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=simil_dicts_2)
    return [peak_1, peak_2, peak_3]

def test_sort_peaks_by_best_match():
    peaks = create_test_peaks()
    sorted_peaks = sort_peaks_by_best_match(peaks)
    assert sorted_peaks[0] == peaks[1] and sorted_peaks[1] == peaks[2] and sorted_peaks[2] == peaks[0]

def test_get_best_match_compound_id():
    peaks = create_test_peaks()
    assert get_best_match_compound_id(peaks[0]) == "3"
    assert get_best_match_compound_id(peaks[1]) == "1"
    assert get_best_match_compound_id(peaks[2]) == "2"


def create_similarity_dicts_2(compound_id_num, n_dicts):
    simil_dicts = []
    for i in range(n_dicts):
        dic = {}
        dic['compound_id'] = str(compound_id_num+i)
        dic['spectrum_correl_coef'] = 1 - 0.01 * (i + 1 / compound_id_num)
        dic['distance'] = 0
        dic['relative_distance'] = 0
        simil_dicts.append(dic)
    return simil_dicts

def create_test_peaks_2():
    simil_dicts_1 = create_similarity_dicts_2(1, 3)
    simil_dicts_2 = create_similarity_dicts_2(2, 2)
    simil_dicts_3 = create_similarity_dicts_2(3, 1)
    
    peak_1 = PreprocessedPeak(left=90, right=110, maximum=100, dataset="a", idx=1,
                              pure=False, saturation=False, integral=0.99, offset=0, matches=simil_dicts_3)
    peak_2 = PreprocessedPeak(left=190, right=210, maximum=200, dataset="b", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=simil_dicts_1)
    peak_3 = PreprocessedPeak(left=290, right=310, maximum=300, dataset="c", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=simil_dicts_2)
    return [peak_1, peak_2, peak_3]  
    
def test_update_peaks_and_matches():
    peaks = create_test_peaks_2()
    sorted_peaks = sort_peaks_by_best_match(peaks)
    
    new_peaks = update_peaks_and_matches(sorted_peaks)
    
    assert len(new_peaks) == 2
    assert peaks[1] in new_peaks and peaks[2] in new_peaks
    for peak in new_peaks:
        assert not any(match['compound_id'] == "3" for match in peak.matches)
        assert any("2" in match['compound_id'] for match in peak.matches)
    
def test_assign_best_match_peak_1():
    peaks = create_test_peaks_2()
    assigned_peak, new_peaks = assign_best_match_peak(peaks)
    assert len(new_peaks) == 2
    assert peaks[1] in new_peaks and peaks[2] in new_peaks
    for peak in new_peaks:
        assert not any(match['compound_id'] == "3" for match in peak.matches)
        assert any("2" in match['compound_id'] for match in peak.matches)
        
    assert type(assigned_peak) == ProcessedPeak
    assert assigned_peak.compound_id == "3" and assigned_peak.maximum == 100

def test_assign_best_match_peak_2():
    peaks = create_test_peaks_2()
    assigned_peak_1, new_peaks = assign_best_match_peak(peaks)
    assigned_peak_2, new_peaks = assign_best_match_peak(new_peaks)
    assigned_peak_3, new_peaks = assign_best_match_peak(new_peaks)
    assert assigned_peak_1.compound_id == "3" and assigned_peak_1.maximum == 100
    assert assigned_peak_2.compound_id == "2" and assigned_peak_2.maximum == 300
    assert assigned_peak_3.compound_id == "1" and assigned_peak_3.maximum == 200
    assert type(new_peaks) == list and len(new_peaks) == 0
    
def test_assign_matched_peaks_1():
    peaks = create_test_peaks_2()
    assigned_peaks, rest_peaks = assign_matched_peaks(peaks)
    assert assigned_peaks[0].compound_id == "3" and assigned_peaks[0].maximum == 100
    assert assigned_peaks[1].compound_id == "2" and assigned_peaks[1].maximum == 300
    assert assigned_peaks[2].compound_id == "1" and assigned_peaks[2].maximum == 200
    assert type(rest_peaks) == list and len(rest_peaks) == 0
    
def test_assign_matched_peaks_2():
    peaks = create_test_peaks_2()
    simil_dicts = create_similarity_dicts_2(1, 3)
    for simil_dict in simil_dicts:
        simil_dict['spectrum_correl_coef'] -= 0.01
    peak_4 = PreprocessedPeak(left=390, right=410, maximum=400, dataset="d", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=simil_dicts)
    peaks.append(peak_4)
    assigned_peaks, rest_peaks = assign_matched_peaks(peaks)
    #logging.warning(f"{rest_peaks} \n {assigned_peaks}")
    assert assigned_peaks[0].compound_id == "3" and assigned_peaks[0].maximum == 100
    assert assigned_peaks[1].compound_id == "2" and assigned_peaks[1].maximum == 300
    assert assigned_peaks[2].compound_id == "1" and assigned_peaks[2].maximum == 200
    assert rest_peaks[0] == peak_4
    
def test_get_next_unknown_id():
    peak_db = PeakDatabase()
    assert get_next_unknown_id(peak_db) == "unknown_1"

def test_assign_unmatched_peaks():
    peak_db = PeakDatabase()
    peaks = []
    peak_5 = PreprocessedPeak(left=380, right=400, maximum=390, dataset="d", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=[])
    peaks.append(peak_5)
    peak_6 = PreprocessedPeak(left=370, right=390, maximum=380, dataset="d", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=[])
    peaks.append(peak_6)
    assigned_peaks = assign_unmatched_peaks_react(peaks, peak_db)
    assert assigned_peaks[0].left == 370 and assigned_peaks[0].compound_id == 'unknown_1'
    assert assigned_peaks[1].left == 380 and assigned_peaks[1].compound_id == 'unknown_2'

def test_assign_peaks_react():
    peak_db = PeakDatabase()
    peaks = create_test_peaks_2()
    simil_dicts = create_similarity_dicts_2(1, 3)
    for simil_dict in simil_dicts:
        simil_dict['spectrum_correl_coef'] -= 0.01
    peak_4 = PreprocessedPeak(left=390, right=410, maximum=400, dataset="d", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=simil_dicts)
    peaks.append(peak_4)
    peak_5 = PreprocessedPeak(left=380, right=400, maximum=390, dataset="d", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=[])
    peaks.append(peak_5)
    peak_6 = PreprocessedPeak(left=370, right=390, maximum=380, dataset="d", idx=1,
                              pure=False, saturation=False, integral=1, offset=0, matches=[])
    peaks.append(peak_6)
    chromatogram = Chromatogram()
    chromatogram.peaks = peaks
    
    chrom = assign_peaks_react(chromatogram, peak_db)
    assert chrom.peaks[0].compound_id == "3" and chrom.peaks[0].maximum == 100
    assert chrom.peaks[1].compound_id == "1" and chrom.peaks[1].maximum == 200
    assert chrom.peaks[2].compound_id == "2" and chrom.peaks[2].maximum == 300
    assert chrom.peaks[3].compound_id == "unknown_1" and chrom.peaks[3].maximum == 380
    assert chrom.peaks[4].compound_id == "unknown_2" and chrom.peaks[4].maximum == 390
    assert chrom.peaks[5].compound_id == "unknown_3" and chrom.peaks[5].maximum == 400
    assert peak_db.unknown_counter == 3
    

    
    









