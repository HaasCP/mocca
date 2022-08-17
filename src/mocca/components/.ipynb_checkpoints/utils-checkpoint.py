#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 13:23:45 2021

@author: haascp
"""
import numpy as np

from mocca.peak.utils import average_peak_spectrum


def check_peaks_compound_id(peaks):
    """
    Checks if all given peaks have the same compound_id and, if so, returns
    this compound_id.
    """
    if all(peak.compound_id == peaks[0].compound_id for peak in peaks):
        compound_id = peaks[0].compound_id
        return compound_id
    else:
        raise AttributeError("All peaks have to have the same compound_id to "
                             "create a component")


def get_valid_peaks(peaks):
    """
    Returns a list of peaks from the database which are pure and unsaturated
    and have a compound_id.
    """
    return [peak for peak in peaks if (peak.pure and
                                       not peak.saturation and
                                       peak.compound_id is not None and
                                       peak.idx > 0)]


# TODO: implement different filter options for relevant_peaks (e.g. by date)
def filter_peaks(peaks, filter_function):
    """
    Filters given peaks with regard to the given filter function (which takes
    a list of peaks and returns a filtered list of peaks).
    """
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
    Returns a filtered (by the given filter function) list of peaks from the
    database which are pure and unsaturated and have a compound_id.
    """
    valid_peaks = get_valid_peaks(peak_database.peaks)
    filtered_peaks = filter_peaks(valid_peaks, filter_function)
    return filtered_peaks


def sort_peaks_by_compound(peaks):
    """
    Returns dict with unique compound_id as keys and a list of corresponding
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
    Creates a filtered (by the given filter function) list of peaks from the
    database which are pure and unsaturated and have a compound_id. From this list,
    it returns a dict with unique compound_id as keys and a list of corresponding
    peaks as values.
    """
    filtered_peaks = get_filtered_peaks(peak_database, filter_function)
    compound_dict = sort_peaks_by_compound(filtered_peaks)
    return compound_dict


def average_spectra_over_peaks(peaks):
    """
    Calculates mean spectrum of a list of peaks with averaged spectrum.
    """
    if peaks:
        spectra_list = []
        for peak in peaks:
            peak_spectrum = average_peak_spectrum(peak)
            spectra_list.append(peak_spectrum)
        return np.average(np.array(spectra_list), axis=0).tolist()
    else:
        return []


def average_ret_times_over_peaks(peaks):
    """
    Calculates mean retention indices of a list of peaks.
    """
    if peaks:
        num_peaks = len(peaks)
        left = int(round(sum([peak.left - peak.offset for peak in peaks]) /
                         num_peaks))
        right = int(round(sum([peak.right - peak.offset for peak in peaks]) /
                          num_peaks))
        maximum = int(round(sum([peak.maximum - peak.offset for peak in peaks]) /
                            num_peaks))
        offset = int(round(sum([peak.offset for peak in peaks]) / num_peaks))
        return left, right, maximum, offset


def get_quant_peaks_by_compound(peak_database, filter_function):
    """
    Returns a dict with compound_ids as keys and lists of peaks as values,
    where only peaks are includud which have is_compound True and which
    have a given concentration.
    """
    filtered_peaks = get_filtered_peaks(peak_database, filter_function)
    quant_peaks = [peak for peak in filtered_peaks if (peak.is_compound and
                                                       peak.concentration)]
    compound_dict = sort_peaks_by_compound(quant_peaks)
    return compound_dict
