#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 11:03:08 2021

@author: haascp
"""
from mocca.peak.models import ProcessedPeak


def process_peak(peak, compound_id, conc):
    return ProcessedPeak(left=peak.left,
                         right=peak.right,
                         maximum=peak.maximum,
                         dataset=peak.dataset,
                         idx=peak.idx,
                         saturation=peak.saturation,
                         pure=peak.pure,
                         integral=peak.integral,
                         offset=peak.offset,
                         compound_id=compound_id,
                         concentration=conc)
