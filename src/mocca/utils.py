#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 14:58:54 2021

@author: haascp
"""

import math

def is_unimodal(L, high_val_threshold = math.inf): 
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