#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 18:55:41 2021

@author: haascp
"""
from operator import attrgetter

from mocca.peak.match import update_matches, match_peak
from mocca.peak.process import process_peak

from mocca.campaign.user_objects import Compound


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
    assigned_peak = process_peak(sorted_peaks[0], Compound(compound_id),
                                 is_compound=False)
    
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
        if peak.matches is None:
            new_peak = process_peak(peak, Compound(None), is_compound=False)
        else:
            new_peak = process_peak(peak,
                                    Compound(get_next_unknown_id(peak_db)),
                                    False)
        assigned_peaks.append(new_peak)
    return assigned_peaks


def get_matched_peaks(peaks):
    return [peak for peak in peaks if peak.matches]


def get_unmatched_peaks(peaks):
    return [peak for peak in peaks if not peak.matches]


def assign_peaks_react(chromatogram, peak_db):
    """
    Assigns peaks of reaction runs with compound ids using unknown compound ids
    for unmatched peaks.
    """
    matched_peaks = get_matched_peaks(chromatogram.peaks)
    unmatched_peaks = get_unmatched_peaks(chromatogram.peaks)
    
    assigned_peaks, unassigned_peaks = assign_matched_peaks(matched_peaks)
    
    unknown_peaks = assign_unmatched_peaks_react(unmatched_peaks + unassigned_peaks,
                                                 peak_db)
    chromatogram.peaks = sorted(assigned_peaks + unknown_peaks,
                                key=lambda peak: peak.maximum)
    return chromatogram


def get_unknown_impurity_peaks(assigned_peaks):
    return [peak for peak in assigned_peaks if
            "unknown" in peak.compound_id or "impurity" in peak.compound_id]


def get_max_integral_peak(peaks):
    """
    Returns the peak with the maximum integral value in the given list of peaks.
    """
    if not peaks:
        return None
    if not all(hasattr(peak, 'integral') for peak in peaks):
        raise AttributeError("All given peaks must have integral attribute.")
    return max(peaks, key=attrgetter('integral'))


def assign_unmatched_peaks_compound(peaks, compound_id, impurity_counter=0):
    """
    Assigns peaks which do not contain matches with unknown compound ids.
    """
    peaks = sorted(peaks, key=lambda peak: peak.maximum)
    assigned_peaks = []
    for peak in peaks:
        impurity_counter += 1
        compound = Compound(compound_id + "_impurity_" + str(impurity_counter))
        new_peak = process_peak(peak, compound, is_compound=False)
        assigned_peaks.append(new_peak)
    return assigned_peaks


def assign_peaks_compound(chromatogram, compound):
    """
    Assigns all matched peaks with compound_ids. If any compound_id of assigned
    peak matches the compound key, compound attributes are assigned to the peak.
    Else, the highest 
    """
    matched_peaks = get_matched_peaks(chromatogram.peaks)
    assigned_peaks, unassigned_peaks = assign_matched_peaks(matched_peaks)
    unmatched_peaks = get_unmatched_peaks(chromatogram.peaks)
    unmatched_peaks = unmatched_peaks + unassigned_peaks
    
    if any([peak.compound_id == compound.key for peak in assigned_peaks]):
        for peak in assigned_peaks:
            if peak.compound_id == compound.key:
                processed_peak = process_peak(peak, compound, is_compound=True)
                assigned_peaks = [p for p in assigned_peaks if p != peak]
    else:
        max_peak = get_max_integral_peak(unmatched_peaks)
        unmatched_peaks = [p for p in unmatched_peaks if p != max_peak]
        if max_peak:
            if not max_peak.pure:
                chromatogram.bad_data = True
                chromatogram.warnings.append("An impure peak was found to be assigned "
                                             "in a pure compound experiment. Run is "
                                             "therefore dismissed.")
                chromatogram.peaks = sorted(assigned_peaks + unmatched_peaks,
                                            key=lambda peak: peak.maximum)
                return chromatogram
            else:
                processed_peak = process_peak(max_peak, compound, is_compound=True)
        else:
            chromatogram.bad_data = True
            chromatogram.warnings.append("No new peak could be found in the data "
                                         "to which the given compound could be "
                                         "assigned.")
            chromatogram.peaks = sorted(assigned_peaks + unmatched_peaks,
                                        key=lambda peak: peak.maximum)
            return chromatogram

    impurity_peaks = assign_unmatched_peaks_compound(unmatched_peaks,
                                                     compound.key)

    
    if processed_peak.saturation:
        chromatogram.warnings.append("Compound was assigned to a peak possibly "
                                     "affected from saturation effect. User "
                                     "check required!")

    chromatogram.peaks = sorted([processed_peak] + assigned_peaks + impurity_peaks,
                                key=lambda peak: peak.maximum)
    return chromatogram


def reassign_impurities(chromatogram, peak_db, quali_comp_db, spectrum_correl_coef_thresh,
                        relative_distance_thresh, print_similarity_dicts=False):
    """
    This function is only allowed to be run in the process_all_experiments function
    which has to be run everytime a new compound should be added to quali_comp_db.
    """
    impurity_peaks = [peak for peak in chromatogram if 'impurity' in peak.compound_id]
    compound_peaks = [peak for peak in chromatogram if peak not in impurity_peaks]
    compound_id = [peak for peak in compound_peaks if peak.is_compound][0].compound_id
    
    # get impurity counter from peak db
    cur_count = 0
    for peak in peak_db:
        if peak.compound_id.startswith(compound_id + "_impurity_"):
            num = int(peak.compound_id[len(compound_id + "_impurity_"):])
            if num > cur_count:
                cur_count = num

    new_peaks = []
    for peak in impurity_peaks:
        matched_peak = match_peak(peak, quali_comp_db, spectrum_correl_coef_thresh,
                                  relative_distance_thresh, print_similarity_dicts)
        new_peaks.append(matched_peak)
    matched_peaks = get_matched_peaks(new_peaks)
    assigned_peaks, unassigned_peaks = assign_matched_peaks(matched_peaks)
    
    unmatched_peaks = get_unmatched_peaks(new_peaks)
    unmatched_peaks = unmatched_peaks + unassigned_peaks
    impurity_peaks = assign_unmatched_peaks_compound(unmatched_peaks,
                                                     compound_id,
                                                     impurity_counter=cur_count)

    chromatogram.peaks = sorted(compound_peaks + assigned_peaks + impurity_peaks,
                                key=lambda peak: peak.maximum)
    return chromatogram
