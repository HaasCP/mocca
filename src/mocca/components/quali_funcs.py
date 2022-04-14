#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:54:04 2021

@author: haascp
"""
from scipy.signal import find_peaks

from mocca.components.utils import (check_peaks_compound_id,
                                    average_ret_times_over_peaks,
                                    average_spectra_over_peaks)
from mocca.components.models import QualiComponent


def get_absorbance_maxima(spectrum):
    """
    Returns absorbance maxima of given spectrum. Maximum must be at least 5%
    intensity of the overall maximum intensity.
    """
    spectrum_maxima, _ = find_peaks(spectrum)
    spectrum_maxima = [m for m in spectrum_maxima if spectrum[m] >
                       0.05 * max(spectrum)]
    return sorted(spectrum_maxima, reverse=True, key=lambda m: spectrum[m])


def create_quali_component(peaks):
    """
    Creates a qualitative component object based on the given peaks.
    """
    # check if all peaks have same compound_id
    if not peaks:
        return
    compound_id = check_peaks_compound_id(peaks)
    mean_left, mean_right, mean_maximum, mean_offset =\
        average_ret_times_over_peaks(peaks)
    mean_spectrum = average_spectra_over_peaks(peaks)
    spectrum_max = get_absorbance_maxima(mean_spectrum)

    if "unknown" in compound_id or "impurity" in compound_id:
        compound_peaks = peaks
    else:
        compound_peaks = [peak for peak in peaks if peak.is_compound]

    return QualiComponent(compound_id=peaks[0].compound_id,
                          left=mean_left,
                          right=mean_right,
                          maximum=mean_maximum,
                          offset=mean_offset,
                          spectrum=mean_spectrum,
                          spectrum_max=spectrum_max,
                          created_from=compound_peaks)
