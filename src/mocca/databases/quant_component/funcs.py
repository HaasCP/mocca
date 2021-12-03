#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:54:04 2021

@author: haascp
"""
from mocca.component.utils import average_retention_times_over_peaks, average_spectra_over_peaks
from mocca.component.model import Component

def create_quant_component_from_peaks(peaks):
    """
    Creates a quantitative component object based on the given peaks
    """
    if all(peak.compound_id == peaks[0].compound_id for peak in peaks):
        return QuantComponent(compound_id = compound_id,
                              spectrum: list,
                              calibration_factor: float,
                              compound_peaks: List[ProcessedPeak],
                              istd_peaks: List[ProcessedPeak] = []
    else:
        raise AttributeError("All peaks have to have the same compound_id to "
                             "create a quantitative component")