#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 09:16:47 2021

@author: haascp
"""
import numpy as np

from mocca.peak.models import PickedPeak

# Peak processing functions
def expand_peak(picked_peak, absorbance_threshold):
    """
    Expands peak boundaries to those actually in the data. It keeps expanding
    them until the absorbance falls below either expand_threshold or
    rel_threshold times the maximum overall absorbance sum.

    Parameters
    ----------
    absorbance_threshold : int
        Minimum absorbance value when picking peaks

    Modifies
    --------
    picked_peak.left, picked_peak.right : The left and right boundaries of the
    peak are updated
    """
    expand_threshold = absorbance_threshold / 20
    # sum absorbances over all wavelengths
    data = np.sum(picked_peak.dataset.data, axis=0)
    # smoothing filter on time axis with a window length of 5
    data = np.convolve(data, np.ones(5), 'same') / 5

    left = picked_peak.left
    right = picked_peak.right
    
    prev_val = np.inf
    while data[left] > expand_threshold and \
            prev_val > data[left] and left >= 0:
        prev_val = data[left]
        left -= 1

    if prev_val != np.inf:  # if peak was expanded, fix boundary, else don't change
        left += 1

    prev_val = np.inf
    while data[right] > expand_threshold and \
            prev_val > data[right] and right <= len(data) - 1:
        prev_val = data[right]
        right += 1

    if prev_val != np.inf:  # if peak was expanded, fix boundary, else don't change
        right -= 1
    
    return PickedPeak(left = left,
                      right = right,
                      maximum = picked_peak.maximum,
                      dataset = picked_peak.dataset,
                      idx = picked_peak.idx)
