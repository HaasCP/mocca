#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 15:30:29 2021

@author: haascp
"""

from mocca.peak.models import CorrectedPeak

def correct_offset(integrated_peak, istd_peaks, offset):
    return CorrectedPeak(left=integrated_peak.left,
                         right=integrated_peak.right,
                         maximum=integrated_peak.maximum,
                         offset=offset,
                         dataset=integrated_peak.dataset,
                         idx=integrated_peak.idx,
                         saturation=integrated_peak.saturation,
                         pure=integrated_peak.pure,
                         integral=integrated_peak.integral,
                         istd=istd_peaks)
