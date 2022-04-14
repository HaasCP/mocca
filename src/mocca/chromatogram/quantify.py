#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 16:41:38 2022

@author: haascp
"""
from mocca.peak.quantify import quantify_peak


def quantify_peaks(chrom, quant_comp_db, quali_comp_db):
    """
    Quantifies all peaks in the chromatogram.
    """
    new_peaks = []
    for peak in chrom:
        new_peak = quantify_peak(peak, quant_comp_db, quali_comp_db)
        new_peaks.append(new_peak)
    chrom.peaks = new_peaks
    return chrom
