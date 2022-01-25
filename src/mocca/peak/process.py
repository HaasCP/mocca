#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 11:03:08 2021

@author: haascp
"""
from mocca.peak.models import ProcessedPeak


def process_peak(peak, compound, is_compound=False):
    return ProcessedPeak(left=peak.left,
                         right=peak.right,
                         maximum=peak.maximum,
                         offset=peak.offset,
                         dataset=peak.dataset,
                         idx=peak.idx,
                         saturation=peak.saturation,
                         pure=peak.pure,
                         integral=peak.integral,
                         istd=peak.istd,
                         compound_id=compound.key,
                         concentration=compound.conc,
                         is_compound=is_compound)
