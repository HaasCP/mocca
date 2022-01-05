#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 16:08:40 2021

@author: haascp
"""

from mocca.peak.expand import expand_peak
from mocca.peak.check import check_peak
from mocca.peak.integrate import integrate_peak
from mocca.chromatogram.correct import correct_istd_offset
from mocca.peak.resolve_impure import get_parafac_peaks
from mocca.peak.match import match_peak


def preprocess_chromatogram(chromatogram, istds, quali_comp_db, 
                            absorbance_threshold, detector_limit, 
                            spectrum_correl_thresh, relative_distance_thresh,
                            print_purity_check = False,
                            print_compound_prediction = False):
    # 1. expand, 2. check, 3. integrate
    integrated_peaks = []
    for picked_peak in chromatogram.peaks:
        expanded_peak = expand_peak(picked_peak, absorbance_threshold)
        checked_peak = check_peak(expanded_peak, detector_limit,
                                  show_analytics=print_purity_check)
        integrated_peak = integrate_peak(checked_peak)
        integrated_peaks.append(integrated_peak)
    
    chromatogram.peaks = integrated_peaks
    
    # 4. correct
    chromatogram = correct_istd_offset(chromatogram, istds, quali_comp_db, 
                                       spectrum_correl_thresh, 
                                       relative_distance_thresh)
    
    # 5. resolve impure
    impure_peaks = [peak for peak in chromatogram if not peak.pure]

    for impure_peak in impure_peaks:
        parafac_peaks = get_parafac_peaks(impure_peak, quali_comp_db,
                                          show_parafac_analytics=False)
        chromatogram.peaks.extend(parafac_peaks)

    # 6. match
    
    matched_peaks = []
    for resolved_peak in chromatogram:
        new_peak = match_peak(resolved_peak, quali_comp_db,
                              spectrum_correl_thresh,
                              relative_distance_thresh,
                              print_compound_prediction)
        matched_peaks.append(new_peak)
    chromatogram.peaks = matched_peaks
    return chromatogram

