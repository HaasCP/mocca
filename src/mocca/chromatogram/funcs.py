#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 10:54:24 2021

@author: haascp
"""

from mocca.peak.preprocessor import preprocess_peak

def preprocess_chromatogram(chromatogram, component_db, absorbance_threshold,
                            detector_limit, spectrum_correl_thresh,
                            relative_distance_thresh):
    preprocessed_peaks = []
    for peak in chromatogram.peaks:
        new_peak = preprocess_peak(peak, component_db, absorbance_threshold,
                                   detector_limit, spectrum_correl_thresh,
                                   relative_distance_thresh)
        preprocessed_peaks.append(new_peak)
    chromatogram.peaks = preprocessed_peaks
    return chromatogram