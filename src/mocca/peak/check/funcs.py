#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 09:16:47 2021

@author: haascp
"""
from mocca.peak.models import CheckedPeak
from mocca.peak.check.purity_predictor import predict_peak_purity

def check_peak_saturation(picked_peak, detector_limit):
    """
    Integrates the peak and sets picked_peak.integral to be that value.
    Parameters
    ----------
    detector_limit : int
        Absorbance values above which detector saturation is expected

    Modifies
    --------
    picked_peak.saturation : Sets peak attribute to either True or False
        based on if the peak absorbance exceeds detector_limit.
    """
    max_absorbance = picked_peak.dataset.data[:, picked_peak.maximum].max()
    return max_absorbance > detector_limit

def check_peak(expanded_peak, detector_limit, param=2.5, show_analytics=False):
    new_saturation = check_peak_saturation(expanded_peak, detector_limit)
    new_pure = predict_peak_purity(expanded_peak, show_analytics)
    
    return CheckedPeak(left = expanded_peak.left,
                       right = expanded_peak.right,
                       maximum = expanded_peak.maximum,
                       dataset = expanded_peak.dataset,
                       idx = expanded_peak.idx,
                       saturation = new_saturation,
                       pure = new_pure)