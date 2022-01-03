#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 10:43:13 2021

@author: haascp
"""

from mocca.dad_data.models import CompoundData

from mocca.dad_data.process_funcs import pick_peaks
from mocca.chromatogram.preprocessor import preprocess_chromatogram

def process_compound_exp(exp, hplc_system_tag, gradient, quali_comp_db, settings):
    compound_data = CompoundData(hplc_system_tag, exp, gradient,
                                 settings.wl_high_pass,
                                 settings.wl_low_pass,
                                 settings.detector_limit)
    chromatogram = pick_peaks(compound_data, settings.absorbance_threshold, 
                              settings.peaks_high_pass,
                              settings.peaks_low_pass)
    chromatogram = preprocess_chromatogram(chromatogram, exp.istd,
                                           quali_comp_db, 
                                           settings.absorbance_threshold, 
                                           settings.detector_limit, 
                                           settings.spectrum_correl_thresh,
                                           settings.relative_distance_thresh)
    
    # assign matched peaks, assign compound peaks, add them to the peak db
    #update component_db after all solvent runs