#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 13:23:45 2021

@author: haascp
"""
import numpy as np

from mocca.peak.utils import average_peak_spectrum
from mocca.components.models import QualiComponent
from mocca.peak.models import ProcessedPeak

def get_closest_peak(component: QualiComponent) -> ProcessedPeak:
    """
    Takes in a QualiComponent and returns the peak it was made from that is the
    closest to the average left, right, and maximum.
    """
    if len(component.created_from) == 0:
        raise ValueError("This QualiComponent has no peaks in created_from attribute!")

    def get_peak_similarity(peak: ProcessedPeak):
        """
        Returns a numeric value computing how similar the input peak is in retention
        time to the input QualiComponent.

        Calculates score by multiplying together the relative offsets of left, right,
        and maximum. The best score is 0, and worse scores are higher in value.
        """
        relative_left_offset = abs((component.left - peak.left) / component.left)
        relative_right_offset = abs((component.right - peak.right) / component.right)
        relative_maximum_offset = abs((component.maximum - peak.maximum) / component.maximum)
        peak_score = relative_left_offset * relative_right_offset * relative_maximum_offset
        return peak_score

    #sort list by peak score and return best score.
    return sorted(component.created_from, key=get_peak_similarity)[0]


def get_valid_peaks(peaks):
    """
    Returns a list of peaks from the database which are pure and unsaturated
    and have a compound_id.
    """
    return [peak for peak in peaks if (peak.pure and
                                       not peak.saturation and
                                       peak.compound_id is not None)]


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
        left = int(round(sum([peak.left for peak in peaks]) / num_peaks))
        right = int(round(sum([peak.right for peak in peaks]) / num_peaks))
        maximum = int(round(sum([peak.maximum for peak in peaks]) / num_peaks))
        return left, right, maximum