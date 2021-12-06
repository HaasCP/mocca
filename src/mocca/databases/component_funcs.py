#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:54:04 2021

@author: haascp
"""
from mocca.databases.utils import average_ret_times_over_peaks, average_spectra_over_peaks
from mocca.databases.component_models import QualiComponent

def create_quali_component(peaks):
    """
    Creates a component object based on the given peaks
    """
    # check if all peaks have same compound_id
    if not peaks:
        return
    if all(peak.compound_id == peaks[0].compound_id for peak in peaks):
        mean_left, mean_right, mean_maximum = average_ret_times_over_peaks(peaks)
        mean_spectrum = average_spectra_over_peaks(peaks)
        return QualiComponent(compound_id=peaks[0].compound_id,
                         left=mean_left,
                         right=mean_right,
                         maximum=mean_maximum,
                         spectrum=mean_spectrum,
                         created_from=peaks)
    else:
        raise AttributeError("All peaks have to have the same compound_id to "
                             "create a component")


def create_quant_component(peaks):
    """
    Creates a quantitative component object based on the given peaks
    """
    pass
#    if all(peak.compound_id == peaks[0].compound_id for peak in peaks):
#        return QuantComponent(compound_id = compound_id,
#                              spectrum: list,
#                              calibration_factor: float,
#                              compound_peaks: List[ProcessedPeak],
#                              istd_peaks: List[ProcessedPeak] = []
#    else:
#        raise AttributeError("All peaks have to have the same compound_id to "
#                             "create a quantitative component")
