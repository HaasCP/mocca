#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 09:21:06 2021

@author: haascp
"""
from mocca.peak.models import PickedPeak, ProcessedPeak
from mocca.peak.process.funcs import expand_peak, check_peak_saturation
from mocca.peak.process.purity_predictor import predict_peak_purity

def process_peak(picked_peak: PickedPeak, absorbance_threshold: int, 
                 detector_limit: int, print_purity_check: bool = False) \
    -> ProcessedPeak:
    """
    Peak processing routine
    """
    new_left, new_right = expand_peak(picked_peak, absorbance_threshold)
    new_saturation = check_peak_saturation(picked_peak, detector_limit)
    new_pure = predict_peak_purity(picked_peak, print_purity_check)
    
    return ProcessedPeak(left = new_left,
                         right = new_right,
                         maximum = picked_peak.maximum,
                         dataset = picked_peak.dataset,
                         idx = picked_peak.idx,
                         saturation = new_saturation,
                         pure = new_pure)
