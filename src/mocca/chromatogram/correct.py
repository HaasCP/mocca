#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 10:54:24 2021

@author: haascp
"""
import copy

from mocca.peak.match import get_filtered_similarity_dicts, get_relative_distance
from mocca.peak.correct import correct_offset


def get_pure_istd_peak(chromatogram, istd_key, quali_component_db, 
                       spectrum_correl_coef_thresh, relative_distance_thresh):
    """
    Doubled relative distance threshold
    """
    istd_peak = None
    best_correl_coef = 0
    pure_peaks = [peak for peak in chromatogram if peak.pure]
    for peak in pure_peaks:
        matches = get_filtered_similarity_dicts(peak, quali_component_db, 
                                                spectrum_correl_coef_thresh,
                                                relative_distance_thresh * 2)
        if matches:
            if any(match['compound_id'] == istd_key for match in matches):
                for match in matches:
                    if match['compound_id'] == istd_key:
                        if match['spectrum_correl_coef'] > best_correl_coef:
                            best_correl_coef = match['spectrum_correl_coef']
                            istd_peak = peak
    return istd_peak


def get_impure_istd_peak(chromatogram, istd_key, quali_component_db, 
                       spectrum_correl_coef_thresh, relative_distance_thresh):
    """
    Doubled relative distance threshold
    """
    istd_component = quali_component_db[istd_key]
    
    impure_peak_targets = [peak for peak in chromatogram if 
                           (not peak.pure and 
                            get_relative_distance(peak, istd_component) <=
                            relative_distance_thresh * 2)]
     # TODO PARAFAC routine on impure peak
     #return istd_peak # None if not found via iterative PARAFAC


def get_istd_peak(chromatogram, istd_key, quali_component_db, 
                  spectrum_correl_coef_thresh, relative_distance_thresh):
    """
    Tries to find an istd peak in the chromatogram from both pure or impure peaks.
    """
    if not istd_key in quali_component_db:
        return None
    else:
        new_quali_component_db = copy.deepcopy(quali_component_db)
        new_quali_component_db.items = [new_quali_component_db[istd_key]]
        istd_peak = get_pure_istd_peak(chromatogram, istd_key, new_quali_component_db, 
                                       spectrum_correl_coef_thresh,
                                       relative_distance_thresh)
        if istd_peak is None:
            istd_peak = get_impure_istd_peak(chromatogram, istd_key,
                                             quali_component_db,
                                             spectrum_correl_coef_thresh,
                                             relative_distance_thresh)
        return istd_peak


def get_istd_offset(chromatogram, istd_key, quali_component_db, 
                    spectrum_correl_coef_thresh, relative_distance_thresh):
    """
    Finds possible istd peak in the chromatogram and calculates the retention
    time offset of the peak compared to its qualitative component in the database.
    """
    istd_offset = 0
    
    if istd_key in quali_component_db:
        istd_peak = get_istd_peak(chromatogram, istd_key, quali_component_db, 
                                  spectrum_correl_coef_thresh,
                                  relative_distance_thresh)
        if istd_peak is not None:
            istd_component = quali_component_db[istd_key]
            istd_offset = istd_peak.maximum - istd_component.maximum
    return istd_offset


def correct_istd_offset(chromatogram, istd_key, quali_component_db, 
                    spectrum_correl_coef_thresh, relative_distance_thresh):
    istd_offset = get_istd_offset(chromatogram, istd_key, quali_component_db, 
                        spectrum_correl_coef_thresh, relative_distance_thresh)
    corrected_peaks = []
    for peak in chromatogram.peaks:
        new_peak = correct_offset(peak, istd_offset)
        corrected_peaks.append(new_peak)
    chromatogram.peaks = corrected_peaks
    return chromatogram
