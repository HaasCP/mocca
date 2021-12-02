#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:05:19 2021

@author: haascp
"""
import numpy as np

from mocca.peak_utils import get_peak_data
def integrate_peak(assigned_peak):
    """
    Integrates the peak.

    Modifies
    --------
    picked_peak.integral : Integral of the peak
    """
    peak_data = get_peak_data(assigned_peak)
    return np.sum(peak_data).tolist()
