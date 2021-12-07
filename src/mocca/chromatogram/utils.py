#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 13:26:09 2021

@author: haascp
"""


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
