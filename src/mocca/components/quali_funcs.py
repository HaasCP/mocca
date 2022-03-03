#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:54:04 2021

@author: haascp
"""
from mocca.components.utils import (check_peaks_compound_id,
                                    average_ret_times_over_peaks,
                                    average_spectra_over_peaks)
from mocca.components.models import QualiComponent


def create_quali_component(peaks):
    """
    Creates a qualitative component object based on the given peaks.
    """
    # check if all peaks have same compound_id
    if not peaks:
        return
    compound_id = check_peaks_compound_id(peaks)

    if "unknown" in compound_id or "impurity" in compound_id:
        mean_left, mean_right, mean_maximum, mean_offset =\
            average_ret_times_over_peaks(peaks)
        mean_spectrum = average_spectra_over_peaks(peaks)
        return QualiComponent(compound_id=peaks[0].compound_id,
                              left=mean_left,
                              right=mean_right,
                              maximum=mean_maximum,
                              offset=mean_offset,
                              spectrum=mean_spectrum,
                              created_from=peaks)
    else:
        mean_left, mean_right, mean_maximum, mean_offset =\
            average_ret_times_over_peaks(peaks)
        compound_peaks = [peak for peak in peaks if peak.is_compound]
        mean_spectrum = average_spectra_over_peaks(compound_peaks)
        return QualiComponent(compound_id=compound_id,
                              left=mean_left,
                              right=mean_right,
                              maximum=mean_maximum,
                              offset=mean_offset,
                              spectrum=mean_spectrum,
                              created_from=compound_peaks)
