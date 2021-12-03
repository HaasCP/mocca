#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 13:23:45 2021

@author: haascp
"""
import numpy as np

from mocca.peak.utils import average_peak_spectrum

def get_valid_peaks(peaks):
    """
    Filters list of peaks for pure and unsaturated peaks with a compound_id.
    """
    return [peak for peak in peaks if (peak.pure is True and
                                       peak.saturation is False and
                                       peak.compound_id is not None)]

# TODO: implement different filter options for relevant_peaks (e.g. by date)
def filter_peaks(peaks, filter_function):
    if filter_function is None:
        return peaks
    else:
        if callable(filter_function):
            return filter_function(peaks)
        else:
            raise ValueError("Given parameter {} not callable/not a "
                             "function".format(filter_function))

def get_filtered_peaks(peak_database, filter_function):
    """
    """
    valid_peaks = get_valid_peaks(peak_database.peaks)
    filtered_peaks = filter_peaks(valid_peaks, filter_function)
    return filtered_peaks

def sort_peaks_by_compound(peaks):
    """
    Creates dict with unique compound_id as keys and a list of corresponding
    peaks as values.
    """
    compound_dict = {}
    for peak in peaks:
        if peak.compound_id not in compound_dict:
            compound_dict[peak.compound_id] = []
        compound_dict[peak.compound_id].append(peak)
    return compound_dict

def get_filtered_peaks_by_compound(peak_database, filter_function):
    """
    """
    filtered_peaks = get_filtered_peaks(peak_database.peaks, filter_function)
    compound_dict = sort_peaks_by_compound(filtered_peaks)
    return compound_dict

def average_spectra_over_peaks(peaks):
    """
    Calculates mean spectrum of a list of peaks with averaged spectrum.
    """
    spectra_list = []
    for peak in peaks:
        peak_spectrum = average_peak_spectrum(peak)
        spectra_list.append(peak_spectrum)
    return np.average(np.array(spectra_list), axis=0)

def average_ret_times_over_peaks(peaks):
    """
    Calculates mean retention indices of a list of peaks.
    """
    num_peaks = len(peaks)
    left = int(round(sum([peak.left for peak in peaks]) / num_peaks)),
    right = int(round(sum([peak.right for peak in peaks]) / num_peaks)),
    maximum = int(round(sum([peak.maximum for peak in peaks]) / num_peaks))
    return left, right, maximum