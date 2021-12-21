#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 18:55:41 2021

@author: haascp
"""
from operator import attrgetter

from mocca.peak.match import update_matches
from mocca.peak.process import process_peak
from mocca.peak.utils import get_retention_time


def sort_peaks_by_best_match(peaks):
    """
    Sorts peaks by descending spectrum correlation coefficient in their matches.
    """
    matched_peaks = [peak for peak in peaks if peak.matches]
    unmatched_peaks = [peak for peak in peaks if not peak.matches]
    sorted_peaks = sorted(matched_peaks, reverse=True,
                          key=lambda peak: peak.matches[0]['spectrum_correl_coef'])
    return sorted_peaks + unmatched_peaks


def get_best_match_compound_id(peak):
    """
    Returns the compound id of the best match of the given peak.
    """
    return peak.matches[0]['compound_id']


def update_peaks_and_matches(sorted_peaks):
    """
    Triggered after peak assignment. Deletes the peak which was assigned and
    removes the consumed compound id from the matches of all remaining peaks.
    """
    compound_id = get_best_match_compound_id(sorted_peaks[0])
    peaks = sorted_peaks
    peaks.pop(0)
    
    new_peaks = []
    if peaks:
        for peak in peaks:
            new_matches = [match for match in peak.matches if
                           match['compound_id'] != compound_id]
            new_peak = update_matches(peak, new_matches)
            new_peaks.append(new_peak)
    return new_peaks


def assign_best_match_peak(peaks):
    """
    Assigns the peak with the best correlation coefficient with the compound id
    and updates all remaining peaks accordingly.
    """
    sorted_peaks = sort_peaks_by_best_match(peaks)

    compound_id = get_best_match_compound_id(sorted_peaks[0])
    assigned_peak = process_peak(sorted_peaks[0], compound_id, None)
    
    new_peaks = update_peaks_and_matches(sorted_peaks)
    return assigned_peak, new_peaks


def assign_matched_peaks(peaks, assigned_peaks=[]):
    """
    Assigns peaks containing matches with compound ids. In the rare case, that
    some peaks will not contain matches anymore, these are given back as
    unassigned and unmatched peaks.
    """
    assigned_peaks = []
    residual_peaks = peaks
    while any(len(peak.matches) > 0 for peak in residual_peaks):
        assigned_peak, residual_peaks = assign_best_match_peak(residual_peaks)
        assigned_peaks.append(assigned_peak)
        if not residual_peaks:
            break
    return assigned_peaks, residual_peaks


def get_next_unknown_id(peak_db):
    """
    Returns the next unknown compound_id.
    """
    peak_db.increment_unknown_counter()
    return "unknown_" + str(peak_db.unknown_counter)


def assign_unmatched_peaks_react(peaks, peak_db):
    """
    Assigns peaks which do not contain matches with unknown compound ids.
    """
    peaks = sorted(peaks, key=lambda peak: peak.maximum)
    peak_db.update_unknown_counter()
    assigned_peaks = []
    for peak in peaks:
        new_peak = process_peak(peak, get_next_unknown_id(peak_db), None)
        assigned_peaks.append(new_peak)
    return assigned_peaks


def get_matched_peaks(chromatogram):
    return [peak for peak in chromatogram if peak.matches]


def get_unmatched_peaks(chromatogram):
    return [peak for peak in chromatogram if not peak.matches]


def assign_peaks_react(chromatogram, peak_db):
    """
    Assigns peaks of reaction runs with compound ids using unknown compound ids
    for unmatched peaks.
    """
    matched_peaks = get_matched_peaks(chromatogram)
    unmatched_peaks = get_unmatched_peaks(chromatogram)
    
    assigned_peaks, unassigned_peaks = assign_matched_peaks(matched_peaks)
    
    unknown_peaks = assign_unmatched_peaks_react(unmatched_peaks + unassigned_peaks,
                                                 peak_db)
    chromatogram.peaks = sorted(assigned_peaks + unknown_peaks,
                                key=lambda peak: peak.maximum)
    return chromatogram


def get_unknown_impurity_peaks(chromatogram):
    matched_peaks = [peak for peak in chromatogram if peak.matches]
    return [peak for peak in matched_peaks 
            if all(("unknown" in match['compound_id'] or
                    "impurity" in match['compound_id'])
                   for match in  peak.matches)]


def get_max_integral_peak(peaks):
    """
    Returns the peak with the maximum integral value in the given list of peaks.
    """
    if not all(hasattr(peak, 'integral') for peak in peaks):
        raise AttributeError("All given peaks must have integral attribute.")
    return max(peaks, key=attrgetter('integral'))


def get_compound_peak(chromatogram):
    """
    Returns unmatched peak with the biggest integral in the chromatogram. In
    this context, 'impurity' and 'unknown' matched peaks are regarded as unmatched.
    """
    unmatched_peaks = get_unmatched_peaks(chromatogram)
    relevant_peaks = unmatched_peaks + get_unknown_impurity_peaks(chromatogram)
    target_peak = get_max_integral_peak(relevant_peaks)
    if not target_peak.pure:
        raise ValueError("Mocca identified peak at {} min as compound peak."
                         "This peak is impure! Make sure, compound peaks are"
                         "pure.".format(get_retention_time(target_peak)))


def add_quali_component(chromatogram, compound_id, peak_db, quali_comp_db):
    """
    Assignes max integral unmatched peak of chromatogram with compound_id
    """
    component_peak = get_compound_peak(chromatogram)
    processed_peak = process_peak(component_peak, None)
    
    peak_db.insert_peak(processed_peak)
    quali_comp_db.insert_by_compound_id(peak_db, compound_id)
