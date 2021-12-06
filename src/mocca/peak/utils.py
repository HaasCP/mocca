#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 10:00:18 2021

@author: haascp
"""
import math
import numpy as np

def get_peak_data(peak):
    """
    Returns absorbance data from the left to the right border of the peak
    for all wavelengths.
    """
    return peak.dataset.data[:, peak.left:(peak.right + 1)]

def average_peak_spectrum(peak):
    """
    Calculates mean spectrum over peak from left to right border.
    """
    return np.average(get_peak_data(peak), axis=1).tolist()

def is_unimodal(L, high_val_threshold=math.inf):
    """
    Checks if a list is unimodal (for use in peak purity).

    Parameters
    ----------
    L : list
        A list to test unimodality for

    high_val_threshold : numeric, optional
        If set, then values above high_val_threshold will not be counted in
        unimodality testing. Default is np.inf (i.e. this threshold is not used).

    Returns
    -------
    TYPE boolean
        True if the list is unimodal ignoring high values; False otherwise.
    """

    passed_turning_point = False
    for idx in range(len(L) - 1):
        if not passed_turning_point:
            if L[idx] <= L[idx + 1] or L[idx] > high_val_threshold:
                continue
            else:
                passed_turning_point = True
        else:
            if L[idx] >= L[idx + 1] or L[idx] > high_val_threshold:
                continue
            else:
                return False
    return True

# Compare two peaks
def check_same_dataset(peak, other):
    """
    Raises Exception if the two peaks are not from the same dataset.
    """
    if peak.dataset != other.dataset:
        raise Exception("Peaks are not from the same dataset, \
                        when comparing peak {} and {}!".format(peak.idx,
                        other.idx))

def check_overlap(peak, other):
    """
    Returns True if peak overlaps with the peak 'other', and False otherwise.
    Raises Exception if the two peaks are not from the same dataset.
    """
    check_same_dataset(peak, other)
    return peak.left <= other.left <= peak.right \
        or other.left <= peak.left <= other.right

def get_distance_between(peak, other):
    """
    Returns the distance from the maxima of peak and the other peak.
    Raises Exception if the two peaks are not from the same dataset.
    """
    check_same_dataset(peak, other)
    return abs(peak.maximum - other.maximum)
  