#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 09:21:06 2021

@author: haascp
"""
from mocca.peak.models import PickedPeak, ProcessedPeak
from mocca.databases.models import ComponentDatabase

from mocca.peak.expand.funcs import expand_peak
from mocca.peak.check.funcs import check_peak
from mocca.peak.assign.funcs import assign_peak

def process_peak_init(picked_peak: PickedPeak, component_db: ComponentDatabase, 
                      absorbance_threshold: int, detector_limit: int,
                      assign_peak_thresh_1: float = 0.9, # spectrum_correl_coef_thresh
                      assign_peak_thresh_2: float = 0.01, # relative_distance_thresh
                      print_purity_check: bool = False, 
                      print_compound_prediction: bool = False) -> ProcessedPeak:
    """
    Peak processing routine
    """
    expanded_peak = expand_peak(picked_peak, absorbance_threshold)
    checked_peak = check_peak(expanded_peak, detector_limit, 
                              show_analytics = print_purity_check)
    assigned_peak = assign_peak(checked_peak, component_db, assign_peak_thresh_1,
                        assign_peak_thresh_2, print_compound_prediction)
    
    
    #return processed_peak
