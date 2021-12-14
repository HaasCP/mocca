#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 15:30:29 2021

@author: haascp
"""

from mocca.peak.models import CorrectedPeak

def correct_offset(integrated_peak, offset):
    return CorrectedPeak(left=integrated_peak.left - offset,
                         right=integrated_peak.right - offset,
                         maximum=integrated_peak.maximum - offset,
                         dataset=integrated_peak.dataset,
                         idx=integrated_peak.idx,
                         saturation=integrated_peak.saturation,
                         pure=integrated_peak.pure,
                         integral=integrated_peak.integral,
                         offset=offset)