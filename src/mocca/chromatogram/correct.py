#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 10:54:24 2021

@author: haascp
"""
import copy

from mocca.peak.match import get_filtered_similarity_dicts
from mocca.peak.correct import correct_offset
from mocca.peak.models import IstdPeak
from mocca.decomposition.utils import check_comp_overlap
from mocca.decomposition.iterative_parafac import iterative_parafac


def get_pure_istd_peak(chromatogram, istd_key, quali_component_db,
                       spectrum_correl_coef_thresh, relative_distance_thresh):
    """
    not doubled relative distance threshold
    """
    istd_peak = None
    best_correl_coef = 0
    pure_peaks = [peak for peak in chromatogram if peak.pure]
    for peak in pure_peaks:
        matches = get_filtered_similarity_dicts(peak, quali_component_db,
                                                spectrum_correl_coef_thresh,
                                                relative_distance_thresh)
        if matches:
            if any(match['compound_id'] == istd_key for match in matches):
                for match in matches:
                    if match['compound_id'] == istd_key:
                        if match['spectrum_correl_coef'] > best_correl_coef:
                            best_correl_coef = match['spectrum_correl_coef']
                            istd_peak = peak
    return istd_peak


def get_impure_istd_peak(chromatogram, istd_key, quali_comp_db, absorbance_threshold,
                         spectrum_correl_coef_thresh, relative_distance_thresh):
    """
    not doubled relative distance threshold
    """
    istd_component = quali_comp_db[istd_key]

    impure_peak_targets = [peak for peak in chromatogram if
                           (not peak.pure and
                            check_comp_overlap(peak, istd_component))]

    istd_peak = None
    best_correl_coef = 0
    for impure_peak in impure_peak_targets:
        parafac_model = iterative_parafac(impure_peak, quali_comp_db,
                                          absorbance_threshold,
                                          relative_distance_thresh,
                                          spectrum_correl_coef_thresh,
                                          show_parafac_analytics=False)
        for peak in parafac_model.peaks:
            matches = get_filtered_similarity_dicts(peak, quali_comp_db,
                                                    spectrum_correl_coef_thresh,
                                                    relative_distance_thresh)
            if matches:
                if any(match['compound_id'] == istd_key for match in matches):
                    for match in matches:
                        if match['compound_id'] == istd_key:
                            if match['spectrum_correl_coef'] > best_correl_coef:
                                best_correl_coef = match['spectrum_correl_coef']
                                istd_peak = peak
    return istd_peak


def get_istd_peak(chromatogram, istd_key, quali_component_db, absorbance_threshold,
                  spectrum_correl_coef_thresh, relative_distance_thresh):
    """
    Tries to find an istd peak in the chromatogram from both pure or impure peaks.
    """
    if istd_key not in quali_component_db:
        return None
    else:
        new_quali_component_db = copy.deepcopy(quali_component_db)
        new_quali_component_db.items = [new_quali_component_db[istd_key]]
        istd_peak = get_pure_istd_peak(chromatogram, istd_key, new_quali_component_db,
                                       spectrum_correl_coef_thresh,
                                       relative_distance_thresh)
        if istd_peak is None:
            istd_peak = get_impure_istd_peak(chromatogram, istd_key,
                                             quali_component_db,
                                             absorbance_threshold,
                                             spectrum_correl_coef_thresh,
                                             relative_distance_thresh)
        return istd_peak


def get_istd_offset(istd_peak, istd_key, quali_component_db):
    """
    Finds possible istd peak in the chromatogram and calculates the retention
    time offset of the peak compared to its qualitative component in the database.
    """
    istd_offset = 0
    if istd_peak is not None:
        istd_component = quali_component_db[istd_key]
        istd_offset = istd_peak.maximum - istd_component.maximum
    return istd_offset


def correct_istd_offset(chromatogram, quali_component_db, absorbance_threshold,
                        spectrum_correl_coef_thresh, relative_distance_thresh):
    """
    Corrects the peaks of the chromatogram by the average of the internal standard
    offsets. Adds the offset to the peak objects.
    """
    istd_peaks = []
    if chromatogram.experiment.istd:
        for istd in chromatogram.experiment.istd:
            istd_p = get_istd_peak(chromatogram, istd.key, quali_component_db,
                                   absorbance_threshold,
                                   spectrum_correl_coef_thresh,
                                   relative_distance_thresh)
            if istd_p:
                istd_offset = get_istd_offset(istd_p, istd.key, quali_component_db)
                istd_peak = IstdPeak(maximum=istd_p.maximum,
                                     integral=istd_p.integral,
                                     offset=istd_offset,
                                     compound_id=istd.key,
                                     concentration=istd.conc)
                istd_peaks.append(istd_peak)
            else:
                chromatogram.bad_data = True
        if istd_peaks:
            istd_offsets = [peak.offset for peak in istd_peaks]
            istd_offset = int(round(sum(istd_offsets)/len(istd_offsets)))
        else:
            istd_offset = 0
    else:
        istd_offset = 0

    corrected_peaks = []
    for peak in chromatogram:
        new_peak = correct_offset(peak, istd_peaks, istd_offset)
        corrected_peaks.append(new_peak)

    chromatogram.peaks = corrected_peaks
    return chromatogram
